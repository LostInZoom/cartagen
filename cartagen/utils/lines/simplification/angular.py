import math
from shapely.geometry import (
    LineString, MultiLineString, LinearRing,
    Polygon, MultiPolygon, GeometryCollection
)
from shapely.geometry.base import BaseGeometry

def simplify_angular(geometry, angle=10.0):
    """
    Simplify a line or polygon by removing vertexes with small angles.
    
    This algorithm, proposed by McMaster :footcite:p:`mcmaster:1987`, eliminates vertices that represent
    very small turning angles (< 10 degrees). This prevents rounding filters from
    destroying sharp curvatures, granting the resulting geometry a 'manually 
    generalised' characteristic where significant bends remain prominent.
    
    Accept Multi geometries. If a polygon is provided, it also applies the 
    thinning to its holes using the same parameters.
    
    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to thin.
        If an open line is provided, the endpoints are preserved.
        If a closed ring or polygon is provided, the angles wrap around.
    angle : float, optional
        Turning-angle threshold in degrees. Vertices creating an exterior 
        angle below this limit will be iteratively removed. Default is 10.0.
    
    Returns
    -------
    LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        Thinned geometry of the same type as input.

    See Also
    --------
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
    simplify_topographic :
        Simplify a line or polygon and mimic hand-made cartographic generalization.
    simplify_visvalingam_whyatt :
        Simplify a line or polygon using an area-based selection.
    simplify_whirlpool :
        Simplify a line or polygon using an epsilon-circle based selection.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 0.1), (2, 0), (2, 5)])
    >>> simplify_angular(line, angle=10.0)
    <LINESTRING (0 0, 2 0, 2 5)>
    """
    if geometry.is_empty:
        return geometry

    gtype = geometry.geom_type

    # ── LineString ─────────────────────────────────────────────────────────
    if gtype == "LineString":
        coords = _thin_coords(list(geometry.coords), angle, is_ring=False)
        return LineString(coords) if len(coords) >= 2 else geometry

    # ── LinearRing ─────────────────────────────────────────────────────────
    if gtype == "LinearRing":
        coords = _thin_coords(list(geometry.coords), angle, is_ring=True)
        return LinearRing(coords) if len(coords) >= 4 else geometry

    # ── MultiLineString ─────────────────────────────────────────────────────
    if gtype == "MultiLineString":
        thinned = []
        for line in geometry.geoms:
            coords = _thin_coords(list(line.coords), angle, is_ring=False)
            if len(coords) >= 2:
                thinned.append(LineString(coords))
        return MultiLineString(thinned) if thinned else geometry

    # ── Polygon ─────────────────────────────────────────────────────────────
    if gtype == "Polygon":
        ext = _thin_coords(list(geometry.exterior.coords), angle, is_ring=True)
        if len(ext) < 4:
            return geometry
        holes = []
        for interior in geometry.interiors:
            h = _thin_coords(list(interior.coords), angle, is_ring=True)
            if len(h) >= 4:
                holes.append(h)
        return Polygon(ext, holes)

    # ── MultiPolygon ────────────────────────────────────────────────────────
    if gtype == "MultiPolygon":
        polys = []
        for poly in geometry.geoms:
            ext = _thin_coords(list(poly.exterior.coords), angle, is_ring=True)
            if len(ext) < 4:
                continue
            holes = []
            for interior in poly.interiors:
                h = _thin_coords(list(interior.coords), angle, is_ring=True)
                if len(h) >= 4:
                    holes.append(h)
            polys.append(Polygon(ext, holes))
        return MultiPolygon(polys) if polys else geometry

    raise TypeError(f"Unsupported geometry type: {gtype}")

def _turning_angle_deg(a, b, c) -> float:
    """
    Compute the turning angle at vertex B in the triad (A, B, C).
    Returns the deviation from 180° (i.e. 0° = perfectly straight,
    90° = right-angle bend, 180° = U-turn).
    """
    ax, ay = a[0] - b[0], a[1] - b[1]
    cx, cy = c[0] - b[0], c[1] - b[1]
    mag_a = math.hypot(ax, ay)
    mag_c = math.hypot(cx, cy)
    if mag_a == 0 or mag_c == 0:
        return 0.0
    cos_angle = (ax * cx + ay * cy) / (mag_a * mag_c)
    cos_angle = max(-1.0, min(1.0, cos_angle))   # clamp for float safety
    angle_between = math.degrees(math.acos(cos_angle))  # 0–180°
    # "turning angle" = deviation from straight (180° = straight line)
    return 180.0 - angle_between


def _thin_coords(coords: list, threshold_deg: float, is_ring: bool) -> list:
    """
    Apply triad-angle thinning to a coordinate sequence.
    Vertices whose turning angle is < threshold_deg are removed.
    Iterates until no more vertices are removed (converges).
    """
    pts = list(coords)

    # For a closed ring: the last coord duplicates the first — strip it
    if is_ring and len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]

    changed = True
    while changed:
        changed = False
        keep = [True] * len(pts)
        n = len(pts)

        for i in range(n):
            if not keep[i]:
                continue
            if is_ring:
                prev_i = (i - 1) % n
                next_i = (i + 1) % n
            else:
                if i == 0 or i == n - 1:
                    continue          # never remove endpoints of an open line
                prev_i = i - 1
                next_i = i + 1

            # Walk backwards/forwards to the nearest kept neighbour
            p = prev_i
            while not keep[p] and p != i:
                p = (p - 1) % n if is_ring else p - 1
            nx = next_i
            while not keep[nx] and nx != i:
                nx = (nx + 1) % n if is_ring else nx + 1

            if p == i or nx == i:
                continue

            angle = _turning_angle_deg(pts[p], pts[i], pts[nx])
            if angle < threshold_deg:
                keep[i] = False
                changed = True

        pts = [p for p, k in zip(pts, keep) if k]

        # Guard: a ring needs at least 3 points, a line at least 2
        min_pts = 3 if is_ring else 2
        if len(pts) < min_pts:
            break

    if is_ring:
        pts = pts + [pts[0]]    # re-close the ring
    return pts