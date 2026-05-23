import math
import numpy as np
from shapely.geometry import Polygon

from cartagen.utils.geometry.polygon import polygon_orientation

def square_polygon_orthogonal(polygon, orientation='mbtr', force_all=True, angle_tolerance=8.0, correct_tolerance=0.6, remove_flat=True):
    """
    Squares a building by forcing it into an orthogonal coordinate system.

    This algorithm described by Hakim and Tsai :footcite:p:`hakim:2026`
    is a heuristic approach that leverages the Manhattan World assumption.
    It identifies the orientation of the building and
    performs a global rotation to align the footprint with the X/Y grid.
    Facades are then classified and merged into strictly horizontal or vertical
    lines. Perfect 90-degree corners are reconstructed via direct line
    intersections before an inverse rotation is applied to restore the building
    to its original geospatial location.

    Parameters
    ----------
    polygon : Polygon
        The polygon to square.
    orientation : str, optional
        The method to calculate the orientation of the polygon. This defines
        the dominant axis used for the global rotation.

        - **'primary'** calculates the orientation of the
          longest side of the provided polygon.
        - **'mbr'** calculates the orientation
          of the long side of the minimum rotated bounding rectangle.
        - **'mbtr'** calculates the orientation
          of the long side of the minimum rotated bounding touching rectangle.
        - **'swo'** or statistical weighted orientation calculates
          the orientation using the statistical weighted orientation method.

    force_all : bool, optional
        If set to True (default), every segment is snapped to whichever
        cardinal direction (H or V) it is closest to, regardless of the
        deviation. The ``angle_tolerance`` and ``correct_tolerance``
        parameters are ignored in this mode.
        Set to False to apply tolerance-based squaring: only segments within
        ``angle_tolerance`` of a cardinal direction are snapped, and vertices
        already within ``correct_tolerance`` of right or flat are left alone.
    angle_tolerance : float, optional
        Only used when ``force_all`` is False. Tolerance in degrees within
        which a segment is snapped to the nearest cardinal direction.
        Segments whose deviation exceeds this value are left unchanged ('U').
    correct_tolerance : float, optional
        Only used when ``force_all`` is False. Tolerance in degrees below
        which a vertex is considered already correct and skipped.
    remove_flat : bool, optional
        If set to True, vertices whose interior angle is flat (close to 180°)
        after squaring are removed from the result.

    Returns
    -------
    Polygon

    See Also
    --------
    square_polygon_naive :
        Squares a polygon according to its orientation, vertex by vertex.
    square_polygon_ls :
        Squares a polygon using the method of least squares.
    polygon_orientation :
        Calculate the orientation of a polygon.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 1), (1.1, 1), (1, 0)])
    >>> square_polygon_orthogonal(polygon)
    <POLYGON ((0 1, 1.05 1, 1.05 0, 0 0, 0 1))>
    """
    if force_all:
        flat_tolerance_rad = math.pi  # remove all flat vertices unconditionally
    else:
        angle_tolerance_rad = angle_tolerance * math.pi / 180
        flat_tolerance_rad = correct_tolerance * math.pi / 180

    # --- Extract unique coordinates (drop the closing duplicate point) ---
    coords = np.array(polygon.exterior.coords)
    if np.allclose(coords[0], coords[-1]):
        coords = coords[:-1]

    N = len(coords)
    if N < 3:
        return polygon

    # Delegate to the shared polygon_orientation utility
    dominant_angle = polygon_orientation(polygon, orientation)

    # --- STEP 2: Global rotation to align the building with the X axis ---
    cos_a = math.cos(-dominant_angle)
    sin_a = math.sin(-dominant_angle)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    coords_rot = coords @ R.T

    # --- STEP 3: Classify segments as Horizontal (H) or Vertical (V) ---
    seg_types = []
    for i in range(N):
        p1 = coords_rot[i]
        p2 = coords_rot[(i + 1) % N]
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        # seg_angle in [0, pi/2]: 0 = perfectly horizontal, pi/2 = perfectly vertical
        seg_angle = abs(math.atan2(abs(dy), abs(dx)))

        dist_h = seg_angle                # deviation from horizontal
        dist_v = math.pi / 2 - seg_angle  # deviation from vertical

        if force_all:
            # Always snap to whichever cardinal direction is closest
            seg_types.append('H' if dist_h <= dist_v else 'V')
        else:
            # Only snap if within tolerance, otherwise leave unclassified
            if dist_h <= angle_tolerance_rad:
                seg_types.append('H')
            elif dist_v <= angle_tolerance_rad:
                seg_types.append('V')
            else:
                seg_types.append('U')

    # --- STEP 4: Merge consecutive segments of the same type into facades ---
    # Walk cyclically starting from the first type-change boundary (seam),
    # so no facade is split across the wrap-around.
    seam = 0
    for i in range(N):
        if seg_types[i] != seg_types[(i + 1) % N]:
            seam = (i + 1) % N
            break

    merged_segments = []
    current_type = seg_types[seam]
    current_pts = [coords_rot[seam]]

    for step in range(N):
        i = (seam + step) % N
        i_next = (i + 1) % N
        t = seg_types[i]

        if t != current_type:
            merged_segments.append({'type': current_type, 'pts': current_pts})
            current_type = t
            current_pts = [coords_rot[i]]

        current_pts.append(coords_rot[i_next])

    merged_segments.append({'type': current_type, 'pts': current_pts})

    # A valid orthogonal building needs at least 4 facades
    if len(merged_segments) < 4:
        return polygon

    # --- STEP 5: Compute the representative position of each facade ---
    for seg in merged_segments:
        pts = np.array(seg['pts'])
        if seg['type'] == 'H':
            seg['val'] = np.mean(pts[:, 1])
        elif seg['type'] == 'V':
            seg['val'] = np.mean(pts[:, 0])
        else:
            seg['val'] = None

    # --- STEP 6: Reconstruct corners via alternating intersections ---
    new_coords_rot = []
    M = len(merged_segments)
    for i in range(M):
        seg_curr = merged_segments[i]
        seg_next = merged_segments[(i + 1) % M]

        if seg_curr['type'] == 'H' and seg_next['type'] == 'V':
            new_coords_rot.append([seg_next['val'], seg_curr['val']])
        elif seg_curr['type'] == 'V' and seg_next['type'] == 'H':
            new_coords_rot.append([seg_curr['val'], seg_next['val']])
        else:
            new_coords_rot.append(seg_curr['pts'][-1].tolist())

    if not new_coords_rot:
        return polygon

    new_coords_rot = np.array(new_coords_rot)

    # --- STEP 7: Inverse rotation to restore the building to its original place ---
    cos_inv = math.cos(dominant_angle)
    sin_inv = math.sin(dominant_angle)
    R_inv = np.array([[cos_inv, -sin_inv], [sin_inv, cos_inv]])
    new_coords = new_coords_rot @ R_inv.T

    # --- STEP 8: Remove flat vertices if requested ---
    if remove_flat:
        filtered = []
        n = len(new_coords)
        for i in range(n):
            p_prev = new_coords[(i - 1) % n]
            p_curr = new_coords[i]
            p_next = new_coords[(i + 1) % n]

            v1 = p_prev - p_curr
            v2 = p_next - p_curr
            norm1 = np.linalg.norm(v1)
            norm2 = np.linalg.norm(v2)

            if norm1 < 1e-10 or norm2 < 1e-10:
                continue

            cos_angle = np.clip(np.dot(v1, v2) / (norm1 * norm2), -1.0, 1.0)
            interior_angle = math.acos(cos_angle)

            if abs(math.pi - interior_angle) <= flat_tolerance_rad:
                continue

            filtered.append(p_curr)

        if len(filtered) >= 3:
            new_coords = np.array(filtered)

    # Close the polygon cleanly for Shapely
    new_coords = np.vstack([new_coords, new_coords[0]])
    return Polygon(new_coords)