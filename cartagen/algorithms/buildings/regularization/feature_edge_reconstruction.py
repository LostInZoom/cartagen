import math
from shapely.geometry import (
    LinearRing,
    LineString,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import polygonize, unary_union

def regularize_building_fer(polygon, length=6.0, area=6.0):
    """
    Regularize a polygon using feature edge reconstruction.

    This algorithm was proposed by Yang *et al.* :footcite:p:`yang:2024`.
    FER enforces orthogonality and clean edge alignment on raw segmentation polygons.
    It first estimates the dominant orientation of a building, then iteratively snaps
    near-parallel and near-perpendicular edges to that principal axis, using
    Douglas-Peucker simplification, minimum-rotated-rectangle fitting, and a
    local reconstruction loop.

    This is slightly modified version of the `SamGeo Python package <https://github.com/opengeos/segment-geospatial>`_
    available `here <https://github.com/opengeos/segment-geospatial/blob/main/samgeo/fer.py>`_.
    It is basically the same but without the osgeo dependencies.

    Parameters
    ----------
    polygon : Polygon or MultiPolygon
        Input building footprint(s) to regularize.
    length : float
        Minimum edge length used during the reconstruction pass (default 6 m).
    area : float
        Polygons whose area is below this value (m²) are discarded (default 6 m²).

    Returns
    -------
    Polygon or MultiPolygon or None
        The regularised geometry, preserving the input type.
        Returns None if the geometry is entirely filtered out.

    See Also
    --------
    enclosing_rectangle :
        Construct an enclosing rectangle from a polygon.
    regularize_building_rectangle :
        Transform a polygon into a rectangle.
    regularize_building_regression
        Regularize a polygon using recursive linear regression.

    Notes
    -----
    FER preserves overall shape well and handles irregular, non-axis-aligned
    buildings. But it assumes roughly rectilinear architecture and tuning length and
    area may be needed for dense urban scenes or very small footprints.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 2), (1, 2), (1, 1), (2, 1), (2, 0), (0, 0)])
    >>> regularize_building_fer(polygon, 1.0, 1.0)
    <POLYGON ((0 0, 0 2, 1.5 2, 1.5 0, 0 0))>
    """
    if polygon is None or polygon.is_empty:
        return None

    if isinstance(polygon, Polygon):
        polygons = [polygon]
    elif isinstance(polygon, MultiPolygon):
        polygons = list(polygon.geoms)
    else:
        raise TypeError(f"Expected Polygon or MultiPolygon, got {type(polygon).__name__}")

    results = []
    for poly in polygons:
        if poly.is_empty or poly.area < area:
            continue
        regularised = _process_ring(
            list(poly.exterior.coords),
            min_length=length,
            min_area=area,
        )
        if regularised is not None and not regularised.is_empty:
            results.append(regularised)

    if not results:
        return None
    if len(results) == 1:
        return results[0]
    return MultiPolygon(results)

# ---------------------------------------------------------------------------
# Internal data structures
# ---------------------------------------------------------------------------

class _Vec:
    """2-D directed segment used throughout the algorithm."""

    __slots__ = ("x1", "y1", "x2", "y2", "index", "x", "y")

    def __init__(self, x1, y1, x2, y2, index=0):
        self.x1, self.y1 = x1, y1
        self.x2, self.y2 = x2, y2
        self.index = index
        self.x = x2 - x1
        self.y = y2 - y1

    # ---- geometry helpers -------------------------------------------------

    def length(self):
        return math.hypot(self.x, self.y)

    def angle_deg(self, other):
        """Angle in degrees between self and *other* (0-180)."""
        denom = self.length() * other.length()
        if denom == 0:
            return 0.0
        cos = (self.x * other.x + self.y * other.y) / denom
        cos = max(-1.0, min(1.0, cos))
        return math.degrees(math.acos(cos))

    def slope(self):
        if self.x1 == self.x2:
            return float("inf")
        return (self.y2 - self.y1) / (self.x2 - self.x1)

    def valid(self):
        return not (self.x1 == self.x2 and self.y1 == self.y2)


class _AppendV:
    __slots__ = ("v1", "v2", "domain", "priority", "fit")

    def __init__(self, v1, v2, domain, priority, fit):
        self.v1 = v1
        self.v2 = v2
        self.domain = domain
        self.priority = priority
        self.fit = fit


# ---------------------------------------------------------------------------
# Low-level geometry utilities
# ---------------------------------------------------------------------------

def _coords_to_vecs(coords):
    """Convert a coordinate sequence (ring-like) to a list of _Vec."""
    pts = list(coords)
    if pts[0] == pts[-1]:           # drop closing duplicate
        pts = pts[:-1]
    vecs = []
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        if (x1, y1) != (x2, y2):
            vecs.append(_Vec(x1, y1, x2, y2, i))
    return vecs


def _vecs_to_coords(vecs):
    """Convert a list of _Vec back to a closed coordinate list."""
    coords = [(v.x1, v.y1) for v in vecs]
    if coords:
        coords.append(coords[0])    # close the ring
    return coords


def _ring_area(coords):
    """Shoelace formula — signed area (positive = CCW)."""
    n = len(coords)
    area = 0.0
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return area / 2.0


def _ring_centroid(coords):
    """Centroid of a polygon given its exterior coordinate list."""
    n = len(coords)
    cx = cy = a = 0.0
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
        a += cross
    a *= 3.0
    if a == 0:
        return (sum(c[0] for c in coords) / n, sum(c[1] for c in coords) / n)
    return cx / a, cy / a


# ---------------------------------------------------------------------------
# Step 1 – Douglas-Peucker compression
# ---------------------------------------------------------------------------

def _compress(pts, p1_idx, p2_idx, delete_ids, D=0.5):
    """Recursive Douglas-Peucker simplification."""
    m, n = p1_idx, p2_idx
    if n == m + 1:
        return

    p1 = pts[m]
    p2 = pts[n]
    A = p1[1] - p2[1]
    B = p2[0] - p1[0]
    C = p1[0] * p2[1] - p2[0] * p1[1]
    denom = math.sqrt(A * A + B * B)
    if denom == 0:
        for i in range(m + 1, n):
            delete_ids.add(i)
        return

    distances = [
        abs(A * pts[i][0] + B * pts[i][1] + C) / denom
        for i in range(m + 1, n)
    ]
    dmax = max(distances)

    if dmax <= D:
        for i in range(m + 1, n):
            delete_ids.add(i)
    else:
        mid_local = distances.index(dmax)
        mid_global = m + 1 + mid_local
        _compress(pts, m, mid_global, delete_ids, D)
        _compress(pts, mid_global, n, delete_ids, D)


def _simplify_ring(coords, D=0.5):
    """Apply Douglas-Peucker to a closed ring; return simplified coords."""
    pts = list(coords)
    if pts[0] == pts[-1]:
        pts = pts[:-1]
    n = len(pts)
    if n < 4:
        return pts

    delete_ids = set()
    _compress(pts, 0, n - 1, delete_ids, D)
    kept = [pt for i, pt in enumerate(pts) if i not in delete_ids]
    return kept


# ---------------------------------------------------------------------------
# Step 2 – Rectangle similarity check / minimum rotated rectangle
# ---------------------------------------------------------------------------

def _rec_similar(polygon: Polygon):
    """
    Return (is_rect_like, simplified_polygon).

    If the polygon is close enough to its minimum rotated rectangle, return
    that rectangle; otherwise return the original polygon unchanged.
    """
    if not polygon.is_valid or polygon.is_empty:
        return False, polygon

    try:
        min_rect = polygon.minimum_rotated_rectangle
    except Exception:
        return False, polygon

    p_area = min_rect.area
    if p_area == 0:
        return False, polygon

    area_diff = p_area / polygon.area - 1

    fits = (
        (p_area < 100 and area_diff < 0.25)
        or (100 <= p_area < 300 and area_diff < 0.20)
        or (300 <= p_area < 500 and area_diff < 0.15)
        or (500 <= p_area < 1000 and area_diff < 0.10)
        or (1000 <= p_area < 2000 and area_diff < 0.10)
        or (2000 <= p_area < 5000 and area_diff < 0.05)
        or (p_area >= 5000 and area_diff < 0.02)
    )
    if fits:
        return True, min_rect
    return False, polygon


# ---------------------------------------------------------------------------
# Step 3 – Principal direction estimation
# ---------------------------------------------------------------------------

def _principal_direction(vecs):
    """Find the orientation axis that minimises the angular deviation sum."""
    if not vecs:
        return _Vec(0, 0, 1, 0)   # default: horizontal

    x0, y0 = vecs[0].x1, vecs[0].y1
    best_vec = _Vec(x0, y0, x0 + 10, y0)
    min_weight = float("inf")

    for deg in range(90):
        rad = math.radians(deg)
        vx = _Vec(x0, y0, x0 + 10 * math.cos(rad), y0 + 10 * math.sin(rad))
        weight = 0.0
        for v in vecs:
            angle = vx.angle_deg(v)
            weight += min(abs(angle), abs(90 - angle), abs(180 - angle)) / 90.0 * v.length()
        if weight < min_weight:
            min_weight = weight
            best_vec = vx

    return best_vec


# ---------------------------------------------------------------------------
# Step 4 – Smooth / straighten near-collinear and near-orthogonal edges
# ---------------------------------------------------------------------------

def _intersect_pt(v1, v2, vx):
    """Project midpoint of the gap between v1 and v2 along principal direction."""
    k = vx.slope()
    if math.isinf(k):
        return v1.x1, v2.y2

    if k == 0:
        x = v1.x1
        y = v2.y2
        x1 = v2.x2
        y1 = v1.y1
    else:
        x = (k * v1.x1 - v1.y1 + (1 / k) * v2.x2 + v2.y2) / (k + 1 / k)
        y = k * (x - v1.x1) + v1.y1
        x1 = ((-1 / k) * v1.x1 - v1.y1 - k * v2.x2 + v2.y2) / (-1 / k - k)
        y1 = (-1 / k) * (x1 - v1.x1) + v1.y1

    xc = (v1.x2 + v2.x1) / 2
    yc = (v1.y2 + v2.y1) / 2
    if math.hypot(x - xc, y - yc) > math.hypot(x1 - xc, y1 - yc):
        return x1, y1
    return x, y


def _line_intersection(p1, p2, p3, p4):
    """Compute intersection of lines (p1-p2) and (p3-p4). Raises if parallel."""
    xdiff = (p1[0] - p2[0], p3[0] - p4[0])
    ydiff = (p1[1] - p2[1], p3[1] - p4[1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        raise ValueError("Lines are parallel — no intersection")

    d = (det(p1, p2), det(p3, p4))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


def _smooth(vecs, vx, angle_diff=25, length_diff=None):
    """
    Merge near-collinear or near-orthogonal consecutive edges.
    Modifies the list in-place (copy beforehand if you need the original).
    """
    if length_diff is None:
        length_diff = float("inf")

    changed = True
    while changed and len(vecs) >= 4:
        changed = False
        n = len(vecs)
        for i in range(n):
            k = (i + 1) % n
            v1, v2 = vecs[i], vecs[k]
            angle = v1.angle_deg(v2)

            # Near-collinear → merge into one vector
            if min(angle, abs(180 - angle)) < angle_diff:
                merged = _Vec(v1.x1, v1.y1, v2.x2, v2.y2, v1.index)
                vecs.pop(max(i, k))
                vecs.pop(min(i, k))
                if merged.valid():
                    vecs.insert(i if i < k else k, merged)
                changed = True
                break

            # Near-perpendicular small edges → snap to a clean right angle
            if (
                abs(90 - angle) < angle_diff - 5
                and abs(90 - angle) > 1
                and v1.length() < length_diff
                and v2.length() < length_diff
            ):
                v1a = v1.angle_deg(vx)
                v2a = v2.angle_deg(vx)
                if abs(90 - v1a) < angle_diff or abs(90 - v2a) < angle_diff:
                    x, y = _intersect_pt(v1, v2, vx)
                    ov1 = _Vec(v1.x1, v1.y1, x, y, v1.index)
                    ov2 = _Vec(x, y, v2.x2, v2.y2, v2.index)
                    vecs.pop(max(i, k))
                    vecs.pop(min(i, k))
                    pos = min(i, k)
                    if ov2.valid():
                        vecs.insert(pos, ov2)
                    if ov1.valid():
                        vecs.insert(pos, ov1)
                    changed = True
                    break

    return vecs


# ---------------------------------------------------------------------------
# Step 5 – Feature-line extraction
# ---------------------------------------------------------------------------

def _feature_line(vecs, vx, angle_diff, min_len1, min_len2):
    """Keep only edges that are sufficiently aligned or long enough."""
    result = []
    for v in vecs:
        a = v.angle_deg(vx)
        min_a = min(a, abs(a - 90), 180 - a)
        if (min_a < angle_diff and v.length() > min_len1) or v.length() > min_len2:
            result.append(v)
    return result


def _fill_gaps(vecs):
    """Insert gap-bridging vectors between non-consecutive feature edges."""
    if not vecs:
        return []
    result = []
    n = len(vecs)
    for i in range(n):
        k = (i + 1) % n
        result.append(vecs[i])
        if (vecs[i].x2, vecs[i].y2) != (vecs[k].x1, vecs[k].y1):
            bridge = _Vec(vecs[i].x2, vecs[i].y2, vecs[k].x1, vecs[k].y1)
            result.append(bridge)
    return result


# ---------------------------------------------------------------------------
# Step 6 – Domain classification & parallel-distance handling
# ---------------------------------------------------------------------------

def _domain(v1, v2, angle_diff, area, length_diff):
    a1 = v1.angle_deg(v2)
    a2 = abs(a1 - 90)
    a3 = 180 - a1

    if a1 < angle_diff and a1 <= a2 and a1 <= a3:
        return 0, a1
    if (
        a2 < angle_diff - 5
        and a2 < a1
        and a2 < a3
        and ((v1.x2 != v2.x1 or v1.y2 != v2.y1) or a2 >= 0.01)
    ):
        if (area >= 200 and v1.length() < length_diff and v2.length() < length_diff) or area < 200:
            return 90, a2
        else:
            if v1.x2 != v2.x1 or v1.y2 != v2.y1:
                return 901, a2
    if a3 < angle_diff and a3 < a1 and a3 < a2:
        return 180, a3

    return -1, a1


def _para_distance(v1, v2, vx, length_diff, area, len2):
    """
    Decide how to reconcile two nearly-parallel edges with a gap between them.
    Returns (x1, y1, x2, y2, decision) where decision:
      -1 → skip, 0 → merge, 1 → snap endpoints, 3 → merge (large area)
    """
    mx = (v1.x2 + v2.x1) / 2
    my = (v1.y2 + v2.y1) / 2

    k = vx.slope()
    if math.isinf(k) or k == 0:
        vn1 = _Vec(mx, v1.y1, mx, my)
        vn2 = _Vec(mx, my, mx, v2.y1)
    else:
        ik = -1 / k
        b = my - ik * mx
        yn = ik * (mx + 1) + b
        vn1 = _Vec(mx + 1, yn, mx, my)
        vn2 = _Vec(mx, my, mx + 1, yn)

    x1, y1 = _intersect_pt(v1, vn1, vx)
    x2, y2 = _intersect_pt(v2, vn2, vx)
    d1 = math.hypot(mx - x1, my - y1)
    d2 = math.hypot(mx - x2, my - y2)

    vp1 = _Vec(v1.x2, v1.y2, x1, y1)
    vp2 = _Vec(v2.x1, v2.y1, x2, y2)
    v12 = _Vec(v1.x2, v1.y2, v2.x1, v2.y1)

    if (vp1.length() == 0 or vp2.length() == 0) and d1 + d2 > length_diff:
        return 0, 0, 0, 0, -1

    if d1 + d2 <= length_diff or v12.length() < length_diff / 1.5:
        if (area >= 100 and v1.length() < len2 and v2.length() < len2) or area < 100:
            return 0, 0, 0, 0, 0
        return 0, 0, 0, 0, 3

    if d1 + d2 > length_diff:
        if (area >= 200 and v1.length() < len2 and v2.length() < len2) or area < 200:
            return x1, y1, x2, y2, 1
        return 0, 0, 0, 0, -1

    return 0, 0, 0, 0, -1


def _line_relation(v1, v2, vx, angle_diff, area, len2):
    domain, angle = _domain(v1, v2, angle_diff, area, len2)
    _, av1 = _domain(v1, vx, angle_diff, area, len2)
    _, av2 = _domain(v2, vx, angle_diff, area, len2)
    priority = v1 if av1 <= av2 else v2
    return _AppendV(priority, v2, domain, priority, angle)


# ---------------------------------------------------------------------------
# Step 7 – Local reconstruction
# ---------------------------------------------------------------------------

def _local_resc(vecs, vx, angle_diff, min_len1, min_len2, area):
    """
    Iteratively adjust edge pairs to enforce orthogonality / parallelism.
    """
    prev_len = -1
    while True:
        n = len(vecs)
        if n < 2 or n == prev_len:
            break
        prev_len = n
        modified = False
        for i in range(n):
            k = (i + 1) % n
            v1, v2 = vecs[i], vecs[k]
            av = _line_relation(v1, v2, vx, angle_diff, area, min_len2 * 3)

            if av.domain == 90:
                x, y = _intersect_pt(v1, v2, vx)
                ov1 = _Vec(v1.x1, v1.y1, x, y, v1.index)
                ov2 = _Vec(x, y, v2.x2, v2.y2, v2.index + 1)
                vecs.pop(max(i, k))
                vecs.pop(min(i, k))
                pos = min(i, k)
                if ov2.valid():
                    vecs.insert(pos, ov2)
                if ov1.valid():
                    vecs.insert(pos, ov1)
                modified = True
                break

            if av.domain == 901:
                try:
                    x, y = _line_intersection(
                        (v1.x1, v1.y1), (v1.x2, v1.y2),
                        (v2.x1, v2.y1), (v2.x2, v2.y2),
                    )
                except ValueError:
                    continue
                ov1 = _Vec(v1.x1, v1.y1, x, y, v1.index)
                ov2 = _Vec(x, y, v2.x2, v2.y2, v2.index + 1)
                vecs.pop(max(i, k))
                vecs.pop(min(i, k))
                pos = min(i, k)
                if ov2.valid():
                    vecs.insert(pos, ov2)
                if ov1.valid():
                    vecs.insert(pos, ov1)
                modified = True
                break

            if av.domain in (0, 180):
                x1, y1, x2, y2, d = _para_distance(
                    v1, v2, vx, min_len1, area, min_len2 * 3
                )
                if d == -1:
                    continue
                if d == 1:
                    if v1.index == v2.index:
                        continue
                    ov1 = _Vec(v1.x1, v1.y1, x1, y1, v1.x1 * min_len1 + v2.x2)
                    ov2 = _Vec(x2, y2, v2.x2, v2.y2, v1.x1 * min_len1 + v2.x2)
                    vecs.pop(max(i, k))
                    vecs.pop(min(i, k))
                    pos = min(i, k)
                    if ov2.valid():
                        vecs.insert(pos, ov2)
                    if ov1.valid():
                        vecs.insert(pos, ov1)
                    modified = True
                    break
                if d in (0, 3):
                    merged = _Vec(v1.x1, v1.y1, v2.x2, v2.y2, v1.index)
                    vecs.pop(max(i, k))
                    vecs.pop(min(i, k))
                    if merged.valid():
                        vecs.insert(min(i, k), merged)
                    modified = True
                    break

        if not modified:
            break

        # Re-fill gaps and re-filter after each pass
        vecs = _feature_line(
            _fill_gaps(vecs), vx, angle_diff, min_len1 + 0.5, min_len2
        )

    return vecs


# ---------------------------------------------------------------------------
# Step 8 – Self-intersection repair
# ---------------------------------------------------------------------------

def _fix_self_intersection(coords):
    """
    If the ring self-intersects, return the largest simple sub-polygon.
    Returns (coords, was_fixed).
    """
    if len(coords) < 4:
        return coords, False

    ls = LineString(coords[:-1] if coords[0] == coords[-1] else coords)
    lr = LineString(list(ls.coords) + [ls.coords[0]])

    if lr.is_simple:
        return coords, False

    mls = unary_union(lr)
    best = None
    best_area = 0.0
    for poly in polygonize(mls):
        if poly.area > best_area:
            best_area = poly.area
            best = poly

    if best is None:
        return coords, False

    return list(best.exterior.coords), True


# ---------------------------------------------------------------------------
# Step 9 – Area-consistency guard
# ---------------------------------------------------------------------------

def _area_ok(orig_coords, new_coords, min_area=6):
    """
    Reject the new ring if its area or centroid differs too much from the original.
    """
    if not orig_coords or not new_coords:
        return False

    try:
        poly_orig = Polygon(orig_coords)
        poly_new = Polygon(new_coords)
    except Exception:
        return False

    i_area = poly_orig.area
    o_area = poly_new.area
    if i_area == 0:
        return True

    cx_i, cy_i = poly_orig.centroid.x, poly_orig.centroid.y
    cx_o, cy_o = poly_new.centroid.x, poly_new.centroid.y
    centroid_dist = math.hypot(cx_i - cx_o, cy_i - cy_o)

    if i_area > 1000 and abs(i_area - o_area) / i_area < 0.25:
        return True
    if abs(i_area - o_area) / i_area >= 0.25 or centroid_dist > 1:
        return False
    return True


# ---------------------------------------------------------------------------
# Core per-ring processing
# ---------------------------------------------------------------------------

def _process_ring(coords, min_length=6, min_area=6):
    """
    Apply FER to a single exterior (or interior) ring.
    Returns a Polygon (possibly simplified/rectangularised) or None to skip.
    """
    # --- Douglas-Peucker simplification ---
    pts = _simplify_ring(coords)
    if len(pts) < 3:
        return None

    try:
        poly = Polygon(pts)
    except Exception:
        return None

    if not poly.is_valid:
        poly = poly.buffer(0)
    if poly.is_empty or poly.area < min_area:
        return None

    # --- Minimum-rectangle test ---
    is_rect, rect_poly = _rec_similar(poly)
    if is_rect:
        return rect_poly

    if len(pts) == 3:   # triangle after compression
        return Polygon(pts)

    # --- Build vector list ---
    vecs = _coords_to_vecs(poly.exterior.coords)
    if len(vecs) < 3:
        return poly

    # --- Principal direction ---
    try:
        vx = _principal_direction(vecs)
    except Exception:
        return poly

    avg_len = poly.exterior.length / len(vecs) * 1.5
    area = poly.area

    # --- Smooth pass ---
    vecs = _smooth(list(vecs), vx, angle_diff=25, length_diff=avg_len)
    smooth_coords = _vecs_to_coords(vecs)

    if not _area_ok(list(poly.exterior.coords), smooth_coords):
        return poly

    if len(vecs) == 3:
        return Polygon(smooth_coords)

    # --- Iterative feature-line reconstruction ---
    ref_coords = smooth_coords
    m = 0.5
    result_coords = smooth_coords

    while m <= min_length:
        vecs = _feature_line(vecs, vx, 20, m, avg_len)
        if len(vecs) < 3:
            m += 10
            continue

        vecs = _local_resc(vecs, vx, 20, m, avg_len, area)
        candidate = _vecs_to_coords(vecs)
        candidate, _ = _fix_self_intersection(candidate)

        if m == 0.5:
            if not _area_ok(smooth_coords, candidate):
                result_coords = smooth_coords
                break
            ref_coords = candidate
        else:
            if not _area_ok(ref_coords, candidate):
                result_coords = ref_coords
                break
            ref_coords = candidate

        result_coords = candidate
        m += 0.5

    result_coords, _ = _fix_self_intersection(result_coords)
    try:
        result = Polygon(result_coords)
        if not result.is_valid:
            result = result.buffer(0)
        return result if not result.is_empty else poly
    except Exception:
        return poly