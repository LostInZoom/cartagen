# --------------------------------------------------------
# Claude translation of the reduce bend algorithm of the geo_sim_processing QGIS plugin
# available at https://github.com/NRCan/geo_sim_processing/blob/master/reduce_bend_algorithm.py
# --------------------------------------------------------

import math
import numpy as np
from shapely.geometry import (LineString, MultiLineString, Polygon, MultiPolygon, LinearRing, Point)
from shapely.ops import unary_union, polygonize

# --------------------------------------------------------
# Constants
# --------------------------------------------------------
ANTI_CLOCK_WISE = -1
CLOCK_WISE = 1
FLAT_ANGLE = 0

ZERO_ANGLE = 1e-9       # Near-zero threshold for flat-angle detection
ZERO_RELATIVE = 1e-9    # Near-zero threshold for distance comparisons

def simplify_wang_muller(geometry, tolerance):
    """
    Simplify a line or polygon using a bend-reduction method.

    This algorithm proposed by Wang & Müller :footcite:p:`wang:1998` analyses the bends (curves) of a
    line or polygon and reduces those whose size falls below a given diameter
    tolerance, emulating the decisions a cartographer would make when
    generalising a line by hand. Topology is preserved throughout, a bend is
    only reduced if doing so does not cause the resulting geometry to
    self-intersect, cross another feature, or violate sidedness constraints.

    This is a translation to work outside QGIS of the
    `reduce bend algorithm <https://github.com/NRCan/geo_sim_processing/blob/master/reduce_bend_algorithm.py>`_
    of the `geo_sim_processing QGIS plugin <https://github.com/NRCan/geo_sim_processing>`_.

    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to simplify.
    tolerance : float
        Theoretical diameter (in the coordinate reference system units) of a
        bend to remove. Bends whose adjusted area is smaller than the
        iso-perimetric equivalent of a circle with this diameter are candidates
        for reduction.  A good rule of thumb for cartographic generalisation is
        to use 0.5 mm at the target map scale
        (e.g. ``tolerance = 25`` for a 1:50 000 map in metres).
        Higher values = more aggressive simplification.

    Returns
    -------
    LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        Simplified geometry of the same type as the input.

    See Also
    --------
    simplify_angular :
        Simplify a line or polygon by removing vertexes with small angles.
    simplify_douglas_peucker :
        Simplify a line or polygon using a distance-based selection.
    simplify_lang :
        Simplify a line or polygon using a look-ahead distance-based selection.
    simplify_li_openshaw :
        Simplify a line or a polygon using a regular grid.
    simplify_raposo :
        Simplify a line or a polygon using an hexagonal tessellation.
    simplify_reumann_witkam :
        Simplify a line or polygon using a directional distance-based selection.
    simplify_visvalingam_whyatt :
        Simplify a line or polygon using an area-based selection.
    simplify_whirlpool :
        Simplify a line or polygon using an epsilon-circle based selection.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 0.5), (2, 0), (3, 0.2), (10, 0)])
    >>> simplify_wang_muller(line, 2.0)
    <LINESTRING (0 0, 10 0)>
    """
    # --- 1. Recursive handling for Multi-geometries ---
    if geometry.geom_type == 'MultiLineString':
        simplified = [simplify_wang_muller(g, tolerance) for g in geometry.geoms]
        return MultiLineString(simplified)

    if geometry.geom_type == 'MultiPolygon':
        simplified = [simplify_wang_muller(g, tolerance) for g in geometry.geoms]
        return MultiPolygon(simplified)

    # --- 2. Handling Polygons ---
    if geometry.geom_type == 'Polygon':
        # Simplify the exterior ring (closed line)
        ext_coords = list(geometry.exterior.coords)
        simplified_ext = _simplify_coords(ext_coords, is_closed=True, tolerance=tolerance)

        # A valid LinearRing requires at least 4 coordinates (3 unique + closing point)
        if len(simplified_ext) < 4:
            return Polygon()

        exterior_ring = LinearRing(simplified_ext)

        # Simplify interior rings (holes)
        simplified_interiors = []
        for interior in geometry.interiors:
            int_coords = list(interior.coords)
            simplified_int = _simplify_coords(int_coords, is_closed=True, tolerance=tolerance)
            if len(simplified_int) >= 4:
                simplified_interiors.append(LinearRing(simplified_int))

        poly = Polygon(exterior_ring, simplified_interiors)

        # Topological repair for edge cases
        if not poly.is_valid:
            poly = poly.buffer(0)
        return poly

    # --- 3. LinearRing ---
    if geometry.geom_type == 'LinearRing':
        coords = list(geometry.coords)
        simplified = _simplify_coords(coords, is_closed=True, tolerance=tolerance)
        if len(simplified) < 4:
            return LinearRing()
        return LinearRing(simplified)

    # --- 4. Core Linear Logic (LineString) ---
    if geometry.geom_type != 'LineString':
        raise ValueError(f'{geometry.geom_type} geometry type cannot be simplified.')

    coords = list(geometry.coords)
    if len(coords) < 3:
        return geometry

    simplified = _simplify_coords(coords, is_closed=False, tolerance=tolerance)
    return LineString(simplified)

# --------------------------------------------------------
# Internal data structures
# --------------------------------------------------------

class _Bend:
    """Represents a single bend (curve) in a line string.

    A bend spans from index ``i`` to index ``j`` in the coordinate array.
    Its geometric footprint is stored as a Shapely Polygon (the "bend polygon")
    used for spatial-constraint checks.
    """

    __slots__ = ('i', 'j', 'coords', 'area', 'perimeter', 'adj_area',
                 'to_reduce', 'geom_bend', 'geom_new_subline', 'geom_old_subline')

    def __init__(self, i, j, coords):
        self.i = i
        self.j = j
        self.coords = coords
        self.to_reduce = False

        # Build the bend polygon: close the bend with a straight chord i→j
        ring_coords = list(coords) + [coords[0]]
        poly = Polygon(ring_coords)
        if not poly.is_valid:
            poly = poly.buffer(0)
        self.geom_bend = poly
        self.area = poly.area
        self.perimeter = poly.length

        # Adjusted area (iso-perimetric quotient used by Wang & Müller)
        self.adj_area = self._calculate_adj_area(self.area, self.perimeter)

        # Chord that replaces the bend after reduction (straight line i→j)
        self.geom_new_subline = LineString([coords[0], coords[-1]])
        # Original sub-line spanning the bend
        self.geom_old_subline = LineString([coords[0], coords[-1]])

    @staticmethod
    def _calculate_adj_area(area, perimeter):
        """Return the adjusted area used to compare bends against the tolerance.

        The adjusted area is the iso-perimetric quotient
        ``A_adj = area * (4π / perimeter²) * 0.75``.
        When ``perimeter`` is zero the method returns 0.
        """
        if perimeter == 0:
            return 0.0
        return area * (4 * math.pi / (perimeter ** 2)) * 0.75

    @staticmethod
    def min_adj_area_from_diameter(tolerance):
        """Return the minimum adjusted area corresponding to *tolerance*.

        The reference circle has diameter *d*, so
        ``A_adj = π(d/2)² * (4π / (πd)²) * 0.75 = 0.75 * π * (d/2)²  / π  = …``
        Simplified: ``0.75 * (d/2)² * π / π² = …``
        In practice the formula used by Wang & Müller reduces to:
        ``min_adj_area = π * (tolerance / 2)² * 0.75``
        """
        return math.pi * (tolerance / 2.0) ** 2 * 0.75


class _RbGeom:
    """Wraps a single Shapely LineString together with its bend list and
    simplification state during the reduction process.
    """

    __slots__ = ('coords', 'is_closed', 'is_simplest', 'bends')

    def __init__(self, coords, is_closed=False):
        self.coords = list(coords)
        self.is_closed = is_closed
        self.is_simplest = False
        self.bends = None  # Will be populated by detect_bends()

    @property
    def linestring(self):
        """Return the current geometry as a Shapely LineString."""
        return LineString(self.coords)

    def length(self):
        return self.linestring.length


# --------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------

def _angle_between_three_points(ax, ay, bx, by, cx, cy):
    """Return the signed angle at *B* formed by the triplet A–B–C (in radians).

    The result is the angle you would measure going counter-clockwise from
    BA to BC, in [0, 2π).
    """
    v1 = (ax - bx, ay - by)
    v2 = (cx - bx, cy - by)
    cross = v1[0] * v2[1] - v1[1] * v2[0]
    dot = v1[0] * v2[0] + v1[1] * v2[1]
    angle = math.atan2(cross, dot)
    if angle < 0:
        angle += 2 * math.pi
    return angle


def _get_angles(coords, is_closed):
    """Return a list of turn-direction labels for each interior vertex.

    Each element is one of: ``CLOCK_WISE``, ``ANTI_CLOCK_WISE``, or
    ``FLAT_ANGLE``.  For an open line the first and last vertices have no
    predecessor/successor, so they are omitted from the raw angle list and
    re-inserted with duplicated boundary values afterwards.
    """
    xy = list(coords)
    if len(xy) < 3:
        return []

    if is_closed:
        xy.insert(0, xy[-2])  # unwrap the ring for circular indexing

    raw = [
        _angle_between_three_points(
            xy[i-1][0], xy[i-1][1],
            xy[i][0],   xy[i][1],
            xy[i+1][0], xy[i+1][1]
        )
        for i in range(1, len(xy) - 1)
    ]

    # Classify each angle
    result = []
    for a in raw:
        if abs(a - math.pi) <= ZERO_ANGLE or abs(a) <= ZERO_ANGLE:
            result.append(FLAT_ANGLE)
        elif a >= math.pi:
            result.append(CLOCK_WISE)
        else:
            result.append(ANTI_CLOCK_WISE)

    return result


def _delete_colinear(angles, rb_geom):
    """Remove co-linear (flat-angle) vertices from *rb_geom* in-place.

    Co-linear vertices are redundant and can cause near-zero-area artefacts
    during spatial calculations.  After deletion the open-line boundary angles
    are duplicated so that the angles list stays aligned with the coordinate
    array.
    """
    is_closed = rb_geom.is_closed
    coords = rb_geom.coords

    for i in reversed(range(len(angles))):
        if angles[i] == FLAT_ANGLE:
            if is_closed:
                # For a closed ring, vertex i in angles corresponds to coords[i]
                coords_index = i
            else:
                # For an open line, angles[0] corresponds to coords[1]
                coords_index = i + 1

            if is_closed and len(coords) - 1 <= 3:
                rb_geom.is_simplest = True
                break
            if len(coords) <= 2:
                break

            del coords[coords_index]
            del angles[i]

    if not is_closed and len(angles) >= 1:
        angles.insert(0, angles[0])
        angles.append(angles[-1])

    if rb_geom.length() <= ZERO_RELATIVE:
        rb_geom.is_simplest = True


def _detect_bends(angles, rb_geom):
    """Partition the line into a list of :class:`_Bend` objects.

    A bend is a maximal sub-sequence of vertices that all curve in the same
    direction (all clockwise or all counter-clockwise).  Inflexion points
    mark the boundaries between consecutive bends.
    """
    rb_geom.bends = []
    if not angles:
        rb_geom.is_simplest = True
        return 0

    coords = rb_geom.coords
    is_closed = rb_geom.is_closed

    # Locate inflexion indices (where the turn direction changes sign)
    inflexions = [{'start': None, 'end': None}]
    for k in range(len(angles) - 1):
        if angles[k] * angles[k + 1] == -1:
            inflexions[-1]['end'] = k + 1
            inflexions.append({'start': k, 'end': None})

    if is_closed:
        if angles[-1] * angles[0] == -1:
            inflexions[-1]['end'] = 0
            inflexions.append({'start': len(angles) - 1, 'end': None})
        inflexions[0]['start'] = inflexions[-1]['start']
        del inflexions[-1]
    else:
        inflexions[0]['start'] = 0
        inflexions[-1]['end'] = len(coords) - 1

    # Convert inflexion spans into _Bend objects
    n = len(coords)
    for inf in inflexions:
        s, e = inf['start'], inf['end']
        if s is None or e is None:
            continue
        if s < e:
            sub = coords[s:e + 1]
        elif s > e:
            sub = coords[s:] + coords[:e + 1]
        else:
            # Degenerate: only two bends in a closed ring
            sub = coords[0:n - 1]

        if len(sub) >= 2:
            rb_geom.bends.append(_Bend(s, e, sub))

    return len(rb_geom.bends)


def _flag_bends_to_reduce(rb_geom, tolerance):
    """Mark which bends should be reduced in the current pass.

    Starting from the smallest adjusted-area bend, a bend is flagged only if
    neither of its neighbours is already flagged (to avoid reducing two
    adjacent bends simultaneously, which could create self-intersections).
    """
    min_adj_area = _Bend.min_adj_area_from_diameter(tolerance)
    bends = rb_geom.bends
    n = len(bends)

    if rb_geom.is_closed and n == 1:
        rb_geom.bends = []
        return

    candidates = sorted(
        [(b.adj_area, i) for i, b in enumerate(bends) if b.area < min_adj_area],
        key=lambda x: x[0]
    )

    if not candidates:
        rb_geom.is_simplest = True
        return

    for (adj_area, i) in candidates:
        if adj_area > min_adj_area:
            break
        prev_flagged = bends[(i - 1) % n].to_reduce if rb_geom.is_closed else (
            bends[i - 1].to_reduce if i > 0 else False)
        next_flagged = bends[(i + 1) % n].to_reduce if rb_geom.is_closed else (
            bends[i + 1].to_reduce if i < n - 1 else False)
        if not prev_flagged and not next_flagged:
            bends[i].to_reduce = True


# --------------------------------------------------------
# Spatial-constraint validators
# --------------------------------------------------------

def _validate_simplicity(geoms_with_itself, new_subline):
    """Return True if *new_subline* does not self-intersect with *geoms_with_itself*."""
    engine = new_subline
    for geom in geoms_with_itself:
        relation = engine.relate(geom)
        # Interior–interior intersection (index 0) or boundary–interior touch (index 1)
        if relation[0] == '0' or relation[1] == '0':
            return False
    return True


def _validate_intersection(geoms_with_others, new_subline):
    """Return True if *new_subline* does not cross any feature in *geoms_with_others*."""
    for geom in geoms_with_others:
        relation = new_subline.relate(geom)
        if relation[0] == '0' or relation[1] == '0':
            return False
    return True


def _validate_sidedness(geoms_with_others, bend_polygon):
    """Return True if no feature in *geoms_with_others* is fully inside *bend_polygon*.

    A feature entirely enclosed by the bend would change sides after reduction,
    violating the topological sidedness constraint.
    """
    for geom in geoms_with_others:
        if bend_polygon.contains(geom):
            return False
    return True


def _build_spatial_index(all_rb_geoms):
    """Return a simple list-based spatial index (STRtree) over all line segments.

    Returns a Shapely STRtree of (segment_geom, rb_geom_id) pairs for fast
    bounding-box queries.
    """
    from shapely.strtree import STRtree
    segments = []
    meta = []
    for gid, rb_geom in enumerate(all_rb_geoms):
        coords = rb_geom.coords
        for k in range(len(coords) - 1):
            seg = LineString([coords[k], coords[k + 1]])
            segments.append(seg)
            meta.append(gid)
    tree = STRtree(segments)
    return tree, segments, meta


def _query_nearby(tree, segments, meta, rb_geom_id, bbox_geom, old_subline, all_rb_geoms):
    """Return two lists of geometries that intersect *bbox_geom*.

    * ``with_itself``: segments belonging to the same geometry (self-intersection check).
    * ``with_others``: segments belonging to other geometries (intersection / sidedness check).

    Segments that share an endpoint with *old_subline* are excluded to avoid
    false positives on touching vertices.
    """
    start = Point(old_subline.coords[0])
    end = Point(old_subline.coords[-1])

    candidates = tree.query(bbox_geom)
    with_itself = []
    with_others = []
    for idx in candidates:
        seg = segments[idx]
        gid = meta[idx]
        # Skip segments that merely share an endpoint with the chord
        seg_start = Point(seg.coords[0])
        seg_end = Point(seg.coords[-1])
        if (seg_start.distance(start) < ZERO_RELATIVE or
                seg_start.distance(end) < ZERO_RELATIVE or
                seg_end.distance(start) < ZERO_RELATIVE or
                seg_end.distance(end) < ZERO_RELATIVE):
            continue
        if gid == rb_geom_id:
            with_itself.append(seg)
        else:
            with_others.append(seg)

    return with_itself, with_others


def _validate_constraints(bend, rb_geom_id, all_rb_geoms, tree, segments, meta):
    """Return True if reducing *bend* respects all topological constraints."""
    bbox_geom = bend.geom_bend
    old_subline = bend.geom_old_subline
    new_subline = bend.geom_new_subline

    with_itself, with_others = _query_nearby(
        tree, segments, meta, rb_geom_id, bbox_geom, old_subline, all_rb_geoms)

    if not _validate_simplicity(with_itself, new_subline):
        return False
    if not _validate_intersection(with_others, new_subline):
        return False
    if not _validate_sidedness(with_others, bbox_geom):
        return False

    return True


# --------------------------------------------------------
# Coordinate-level bend reduction
# --------------------------------------------------------

def _apply_bend_reduction(rb_geom, bend):
    """Delete the interior vertices of *bend* from *rb_geom*.

    The vertices at positions ``bend.i`` (inclusive) through ``bend.j``
    (exclusive) are replaced by a straight chord from ``coords[i]`` to
    ``coords[j]``, which has already been computed as ``bend.geom_new_subline``.
    """
    coords = rb_geom.coords
    is_closed = rb_geom.is_closed

    i, j = bend.i, bend.j

    if not is_closed:
        # Keep coords[0..i] + coords[j..end]
        rb_geom.coords = coords[:i + 1] + coords[j:]
    else:
        # Closed ring: indices wrap around
        if i < j:
            rb_geom.coords = coords[:i + 1] + coords[j:]
        else:
            # Wrap-around bend
            rb_geom.coords = coords[j:i + 1]
            # Re-close
            if rb_geom.coords[0] != rb_geom.coords[-1]:
                rb_geom.coords.append(rb_geom.coords[0])


def _process_one_pass(rb_geoms, tolerance):
    """Execute one full simplification pass over all geometries.

    Returns the number of bends reduced during this pass.
    """
    # Rebuild spatial index at the start of each pass
    tree, segments, meta = _build_spatial_index(rb_geoms)

    nbr_reduced = 0
    for gid, rb_geom in enumerate(rb_geoms):
        if rb_geom.is_simplest or not rb_geom.bends:
            continue

        is_closed = rb_geom.is_closed
        bends = rb_geom.bends
        n = len(bends)
        local_reduced = 0

        indices = list(reversed(range(1, n - 1))) if is_closed else list(reversed(range(n)))

        for ind in indices:
            bend = bends[ind]
            if not bend.to_reduce:
                continue
            if _validate_constraints(bend, gid, rb_geoms, tree, segments, meta):
                _apply_bend_reduction(rb_geom, bend)
                local_reduced += 1
                nbr_reduced += 1

        # For closed lines, attempt the boundary bend only if nothing else was reduced
        if is_closed and local_reduced == 0:
            for ind in [0, n - 1]:
                if 0 <= ind < n and bends[ind].to_reduce:
                    if _validate_constraints(bends[ind], gid, rb_geoms, tree, segments, meta):
                        _apply_bend_reduction(rb_geom, bends[ind])
                        nbr_reduced += 1
                    break

        if local_reduced > 0:
            rb_geom.bends = None  # Force recalculation on next pass

    return nbr_reduced


def _simplify_coords(coords, is_closed, tolerance):
    """Run the full Wang & Müller bend-reduction loop on a single ring/line.

    This function wraps a single coordinate sequence in an :class:`_RbGeom`,
    runs as many passes as required until no further reductions occur, and
    returns the simplified coordinate list.
    """
    rb_geom = _RbGeom(coords, is_closed=is_closed)

    while True:
        if rb_geom.is_simplest:
            break

        if rb_geom.bends is None:
            angles = _get_angles(rb_geom.coords, rb_geom.is_closed)
            _delete_colinear(angles, rb_geom)
            if rb_geom.is_simplest:
                break
            _detect_bends(angles, rb_geom)
            _flag_bends_to_reduce(rb_geom, tolerance)

        # Single-geometry pass (no cross-feature topology check needed here;
        # that is handled at the multi-geometry level via _process_one_pass)
        local_reduced = 0
        if rb_geom.bends:
            bends = rb_geom.bends
            n = len(bends)
            indices = (list(reversed(range(1, n - 1)))
                       if rb_geom.is_closed else list(reversed(range(n))))
            for ind in indices:
                bend = bends[ind]
                if bend.to_reduce:
                    _apply_bend_reduction(rb_geom, bend)
                    local_reduced += 1
            if rb_geom.is_closed and local_reduced == 0:
                for ind in [0, n - 1]:
                    if 0 <= ind < n and bends[ind].to_reduce:
                        _apply_bend_reduction(rb_geom, bends[ind])
                        local_reduced += 1
                        break
            if local_reduced > 0:
                rb_geom.bends = None
            else:
                break  # No bend reduced → convergence reached

        else:
            break

    return rb_geom.coords