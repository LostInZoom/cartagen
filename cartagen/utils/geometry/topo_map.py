"""
Python port (simplified) of GeOxygene's CarteTopo using Shapely for geometry.
This module implements a lightweight topological map with Nodes, Arcs, and Faces.
It is NOT a line-by-line translation; instead it provides the core behavior:
- store nodes and arcs (directed graph)
- sort arcs around nodes by azimuth
- polygonize arcs to build Faces (including an optional "infinite" face)
- left/right face assignment for each arc when possible

Dependencies:
    shapely >= 2.0

Author: Ported for demonstration from "CarteTopo.java" (COGIT / IGN)
"""

from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Iterable, Set
import math
from itertools import combinations

# Shapely imports (runtime requirement)
from shapely.geometry import Point, LineString, Polygon, LinearRing
from shapely.ops import polygonize_full, unary_union, snap
from shapely.ops import shared_paths, linemerge


def _azimuth(p0: Point, p1: Point) -> float:
    """Return azimuth in radians from p0 to p1 in [0, 2*pi)."""
    dx = p1.x - p0.x
    dy = p1.y - p0.y
    ang = math.atan2(dy, dx)
    if ang < 0:
        ang += 2.0 * math.pi
    return ang


@dataclass
class Node:
    id: int
    point: Point
    incident_arcs: Set[int] = field(default_factory=set)

    def degree(self) -> int:
        return len(self.incident_arcs)


@dataclass
class Arc:
    id: int
    geom: LineString
    start: int  # Node id
    end: int    # Node id
    # Faces on each side (ids). By convention, the face on the LEFT of the arc direction is left_face.
    left_face: Optional[int] = None
    right_face: Optional[int] = None
    # Arbitrary attributes
    attrs: Dict[str, object] = field(default_factory=dict)

    def reversed(self) -> "Arc":
        """Return a shallow reversed copy (geometry reversed, start/end swapped, faces swapped)."""
        rev = LineString(list(self.geom.coords)[::-1])
        return Arc(
            id=self.id,
            geom=rev,
            start=self.end,
            end=self.start,
            left_face=self.right_face,
            right_face=self.left_face,
            attrs=dict(self.attrs),
        )

    def length(self) -> float:
        return float(self.geom.length)


@dataclass
class Face:
    id: int
    polygon: Optional[Polygon]  # None is used for the infinite face
    boundary_arcs: List[Tuple[int, bool]] = field(default_factory=list)
    # boundary_arcs: (arc_id, forward=True if oriented from start->end along the boundary ring)


class TopoMap:
    """
    A lightweight topological map inspired by CarteTopo in GeOxygene.

    Key notes vs the original Java:
      - We keep everything in-memory (no persistence or GUI listeners).
      - Geometry operations are handled via Shapely.
      - Face building uses shapely.ops.polygonize_full over the union of arcs.
      - Left/right assignment: for each polygon ring, directed edges imply the
        polygon is on the left; we match edges back to arcs with a tolerance.
    """

    def __init__(self, snap_tolerance: float = 0.0, build_infinite_face: bool = True):
        self.snap_tolerance = float(snap_tolerance)
        self.build_infinite_face = bool(build_infinite_face)

        self.nodes: Dict[int, Node] = {}
        self.arcs: Dict[int, Arc] = {}
        self.faces: Dict[int, Face] = {}

        self._next_node_id = 1
        self._next_arc_id = 1
        self._next_face_id = 1

    # ------------------------ building blocks ------------------------

    def _find_or_create_node(self, pt: Point) -> int:
        """Return an existing node id within snap_tolerance of pt, else create a new one."""
        # Brute force search; for large datasets you should index nodes spatially.
        for nid, node in self.nodes.items():
            if node.point.distance(pt) <= self.snap_tolerance:
                return nid
        nid = self._next_node_id
        self._next_node_id += 1
        self.nodes[nid] = Node(id=nid, point=Point(pt.x, pt.y))
        return nid

    def add_arc_geometry(self, geom: LineString, attrs: Optional[Dict[str, object]] = None) -> int:
        """
        Add an arc given a linestring geometry. The endpoints will be snapped to existing nodes
        within tolerance, or new nodes will be created.
        Returns the arc id.
        """
        if geom.is_empty:
            raise ValueError("Cannot add empty geometry as an arc")
        if len(geom.coords) < 2:
            raise ValueError("Arc geometry must have at least 2 coordinates")
        start_id = self._find_or_create_node(Point(geom.coords[0]))
        end_id = self._find_or_create_node(Point(geom.coords[-1]))

        aid = self._next_arc_id
        self._next_arc_id += 1
        arc = Arc(id=aid, geom=geom, start=start_id, end=end_id, attrs=dict(attrs or {}))
        self.arcs[aid] = arc
        self.nodes[start_id].incident_arcs.add(aid)
        self.nodes[end_id].incident_arcs.add(aid)
        return aid

    def sort_arcs_around_nodes(self) -> Dict[int, List[int]]:
        """
        Sort incident arcs around each node by azimuth of the first segment away from the node.
        Returns a dictionary node_id -> ordered list of arc ids (counter-clockwise from East).
        """
        order: Dict[int, List[int]] = {}
        for nid, node in self.nodes.items():
            pairs = []
            for aid in node.incident_arcs:
                arc = self.arcs[aid]
                if arc.start == nid:
                    # use first segment from start
                    c0, c1 = arc.geom.coords[0], arc.geom.coords[1]
                else:
                    # use last segment toward end (reverse direction away from the node)
                    c0, c1 = arc.geom.coords[-1], arc.geom.coords[-2]
                ang = _azimuth(Point(c0), Point(c1))
                pairs.append((ang, aid))
            pairs.sort(key=lambda t: t[0])  # CCW from East (x positive)
            order[nid] = [aid for _, aid in pairs]
        return order

    # ------------------------ face building ------------------------

    def build_faces(self) -> None:
        """
        Polygonize all arcs and populate self.faces.
        Also attempts to assign left/right faces to arcs using ring orientation.
        """
        if not self.arcs:
            return

        # Union of all arcs (snapped if requested)
        lines = [a.geom for a in self.arcs.values()]
        if self.snap_tolerance > 0:
            # snap lines to themselves to close near-coincident endpoints
            merged = unary_union(lines)
            lines_u = snap(merged, merged, self.snap_tolerance)
        else:
            lines_u = unary_union(lines)

        polygons, dangles, cuts, invalids = polygonize_full(lines_u)

        # Create faces for each polygon
        face_id_map: Dict[Tuple[float, float, float, float], int] = {}
        for poly in polygons.geoms if hasattr(polygons, "geoms") else []:
            fid = self._next_face_id
            self._next_face_id += 1
            self.faces[fid] = Face(id=fid, polygon=poly, boundary_arcs=[])
            # quick bbox key map for later lookup
            face_id_map[tuple(poly.bounds)] = fid

        # Optionally create an infinite face with id 0 (convention)
        if self.build_infinite_face:
            if 0 not in self.faces:
                self.faces[0] = Face(id=0, polygon=None, boundary_arcs=[])

        # Build an index of arcs by normalized segment list for matching
        # (coordinates rounded to reduce floating precision issues)
        def _round_pt(c, nd=9):  # nd=9 ~ mm at continental scales
            return (round(c[0], nd), round(c[1], nd))

        def _norm_seg(a: Tuple[Tuple[float, float], Tuple[float, float]]):
            (x1, y1), (x2, y2) = a
            return (_round_pt((x1, y1)), _round_pt((x2, y2)))

        arc_seg_index: Dict[Tuple[Tuple[float, float], Tuple[float, float]], List[Tuple[int, bool]]] = {}
        for aid, arc in self.arcs.items():
            coords = list(map(_round_pt, arc.geom.coords))
            for i in range(len(coords) - 1):
                seg = (coords[i], coords[i + 1])
                arc_seg_index.setdefault(_norm_seg(seg), []).append((aid, True))
                # Also index reversed for matching reversed orientation
                rseg = (coords[i + 1], coords[i])
                arc_seg_index.setdefault(_norm_seg(rseg), []).append((aid, False))

        # For each polygon ring (exterior and holes), walk edges and map to arcs
        def _register_ring(face_id: int, ring: LinearRing, is_exterior: bool):
            coords = list(map(_round_pt, ring.coords))
            for i in range(len(coords) - 1):  # last == first
                seg = (coords[i], coords[i + 1])
                key = _norm_seg(seg)
                if key not in arc_seg_index:
                    continue  # unmatched tiny segment; can happen after snapping
                # prefer forward matches if any
                candidates = arc_seg_index[key]
                # pick first candidate
                aid, forward = candidates[0]
                # For an exterior ring, polygon interior is on the LEFT of its directed edges if ring is CCW
                # Shapely uses CCW=True for standard exterior rings in valid polygons.
                # If ring is CCW, edges travel CCW; interior is on their LEFT.
                # If ring is CW, interior is on their RIGHT.
                ring_ccw = ring.is_ccw
                if ring_ccw:
                    interior_on_left = True
                else:
                    interior_on_left = False

                # If the segment orientation equals the arc forward orientation, then:
                # - if interior_on_left: the polygon is on LEFT of the arc
                # - else: polygon is on RIGHT of the arc
                # If reversed, flip.
                if forward:
                    if interior_on_left:
                        self._set_arc_face_left(aid, face_id)
                    else:
                        self._set_arc_face_right(aid, face_id)
                    self.faces[face_id].boundary_arcs.append((aid, True))
                else:
                    if interior_on_left:
                        self._set_arc_face_right(aid, face_id)
                    else:
                        self._set_arc_face_left(aid, face_id)
                    self.faces[face_id].boundary_arcs.append((aid, False))

        for fid, face in list(self.faces.items()):
            if face.polygon is None:
                continue
            poly = face.polygon
            _register_ring(fid, LinearRing(poly.exterior.coords), True)
            for interior in poly.interiors:
                _register_ring(fid, LinearRing(interior.coords), False)

        # Assign infinite face to arcs sides that remain unassigned
        if self.build_infinite_face and 0 in self.faces:
            for arc in self.arcs.values():
                if arc.left_face is None:
                    arc.left_face = 0
                    self.faces[0].boundary_arcs.append((arc.id, True))
                if arc.right_face is None:
                    arc.right_face = 0
                    self.faces[0].boundary_arcs.append((arc.id, False))

    def _set_arc_face_left(self, aid: int, face_id: int) -> None:
        arc = self.arcs[aid]
        if arc.left_face is None:
            arc.left_face = face_id

    def _set_arc_face_right(self, aid: int, face_id: int) -> None:
        arc = self.arcs[aid]
        if arc.right_face is None:
            arc.right_face = face_id

    # ------------------------ queries ------------------------

    def adjacent_nodes(self, nid: int) -> Set[int]:
        """Return the set of node ids adjacent to nid via arcs."""
        out: Set[int] = set()
        for aid in self.nodes[nid].incident_arcs:
            arc = self.arcs[aid]
            other = arc.end if arc.start == nid else arc.start
            out.add(other)
        return out

    def arcs_around_node_sorted(self, nid: int) -> List[int]:
        """Convenience to sort arcs around a single node."""
        return self.sort_arcs_around_nodes().get(nid, [])

    def build_from_polygons(self, polygons: Iterable[Polygon]) -> None:
        """
        Construct the topological map from a list of polygons.

        Nodes are created only at endpoints of unique arcs (where boundaries start/stop being shared).
        Shared boundaries between polygons are represented by a single arc.
        """
        lines = []
        for boundary in polygons:
            lines.append(boundary.exterior)
    
        topo_arcs = []
        for line1, line2 in combinations(lines, 2):
            if line1.intersects(line2) != True:
                continue
            
            shared = shared_paths(line1, line2)
            for multi in shared.geoms:
                merged = linemerge(multi)
                if isinstance(merged, LineString):
                    topo_arcs.append(merged)
                else:    
                    for geom in multi.geoms:
                        if isinstance(geom, LineString):
                            topo_arcs.append(geom)

        # now get the remaining arcs that are not shared
        for line in lines:
            current_line = line
            for arc in topo_arcs:
                if current_line.intersects(arc):
                    current_line = current_line.difference(arc)
                    if current_line.is_empty:
                        break
            if current_line.is_empty != True:
                if isinstance(current_line, LineString):
                    topo_arcs.append(current_line)
                else:
                    for geom in current_line.geoms:
                        if isinstance(geom, LineString):
                            topo_arcs.append(geom)
        # add the geometries as arcs
        for arc_geom in topo_arcs:
            self.add_arc_geometry(arc_geom)
        # create the faces from the arcs
        self.sort_arcs_around_nodes()
        self.build_faces()



    def face_of_point(self, pt: Point) -> Optional[int]:
        """Return the id of the face that contains the point (None if no face and infinite face disabled)."""
        for fid, face in self.faces.items():
            if face.polygon is None:
                continue
            if face.polygon.contains(pt):
                return fid
        return 0 if self.build_infinite_face and 0 in self.faces else None

    def clear_faces(self) -> None:
        """
        Remove all faces from the topological map.
        """
        self.faces.clear()
        self._next_face_id = 1


    # ------------------------ export helpers ------------------------

    def to_geojson(self) -> Dict[str, object]:
        """Export nodes, arcs, and (finite) faces to a minimal GeoJSON-like dict."""
        def feat(geom, props):
            return {"type": "Feature", "geometry": geom.__geo_interface__, "properties": props}

        nodes_fc = {
            "type": "FeatureCollection",
            "features": [feat(n.point, {"id": n.id, "degree": n.degree()}) for n in self.nodes.values()],
        }
        arcs_fc = {
            "type": "FeatureCollection",
            "features": [
                feat(a.geom, {"id": a.id, "start": a.start, "end": a.end, "left": a.left_face, "right": a.right_face, **a.attrs})
                for a in self.arcs.values()
            ],
        }
        faces_fc = {
            "type": "FeatureCollection",
            "features": [
                feat(f.polygon, {"id": f.id}) for f in self.faces.values() if f.polygon is not None
            ],
        }
        return {"nodes": nodes_fc, "arcs": arcs_fc, "faces": faces_fc}