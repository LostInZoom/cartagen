import math
import numpy as np
from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

def remove_flat_vertices(geom, tolerance=5.0):
    """
    Removes flat vertices from a geometry.

    A vertex is considered flat when the interior angle it forms with its
    two neighbours is close to 180°, meaning the three points are nearly
    collinear and the vertex brings no geometric information. Removing such
    vertices simplifies the geometry without altering its shape.

    Parameters
    ----------
    geom : LineString, MultiLineString, Polygon, MultiPolygon, or LinearRing
        The geometry to simplify.
    tolerance : float, optional
        Tolerance in degrees. A vertex whose interior angle deviates from
        180° by less than or equal to this value is removed.
        Defaults to 0.6°.

    Returns
    -------
    LineString, MultiLineString, Polygon, MultiPolygon, or LinearRing
        A simplified geometry of the same type as the input.

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1)])
    >>> remove_flat_vertices(polygon)
    <POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))>
    """
    if geom.is_empty:
        return geom

    tolerance_rad = tolerance * math.pi / 180

    def _simplify_coords(coords, is_ring):
        coords_arr = np.array(coords)
        if len(coords_arr) == 0:
            return coords_arr

        # Check if coordinates explicitly close back on themselves
        is_closed = len(coords_arr) > 1 and np.allclose(coords_arr[0], coords_arr[-1])

        # If it's a ring, drop the closing duplicate to calculate wrapping angles properly
        if is_ring and is_closed:
            coords_arr = coords_arr[:-1]

        n = len(coords_arr)
        if is_ring and n < 3:
            return coords
        if not is_ring and n < 2:
            return coords

        filtered = []
        for i in range(n):
            # For open linestrings, always keep the first and last endpoints
            if not is_ring and (i == 0 or i == n - 1):
                filtered.append(coords_arr[i])
                continue

            p_prev = coords_arr[(i - 1) % n]
            p_curr = coords_arr[i]
            p_next = coords_arr[(i + 1) % n]

            v1 = p_prev - p_curr
            v2 = p_next - p_curr
            norm1 = np.linalg.norm(v1)
            norm2 = np.linalg.norm(v2)

            # Degenerate vertex (zero-length edge): remove it
            if norm1 < 1e-10 or norm2 < 1e-10:
                continue

            cos_angle = np.clip(np.dot(v1, v2) / (norm1 * norm2), -1.0, 1.0)
            interior_angle = math.acos(cos_angle)

            # Keep the vertex only if its angle is not flat
            if abs(math.pi - interior_angle) > tolerance_rad:
                filtered.append(p_curr)

        # Enforce minimum lengths and closures
        if is_ring:
            if len(filtered) < 3:
                return coords
            filtered.append(filtered[0])  # Re-close the ring
        else:
            if len(filtered) < 2:
                return coords

        return np.array(filtered)

    # Route processing based on geometry type
    gtype = geom.geom_type

    if gtype == 'LineString':
        return LineString(_simplify_coords(geom.coords, is_ring=False))
        
    elif gtype == 'LinearRing':
        return LinearRing(_simplify_coords(geom.coords, is_ring=True))
        
    elif gtype == 'Polygon':
        ext = _simplify_coords(geom.exterior.coords, is_ring=True)
        ints = [_simplify_coords(ring.coords, is_ring=True) for ring in geom.interiors]
        return Polygon(ext, ints)
        
    elif gtype == 'MultiLineString':
        lines = [LineString(_simplify_coords(line.coords, is_ring=False)) for line in geom.geoms]
        return MultiLineString(lines)
        
    elif gtype == 'MultiPolygon':
        polys = []
        for poly in geom.geoms:
            ext = _simplify_coords(poly.exterior.coords, is_ring=True)
            ints = [_simplify_coords(ring.coords, is_ring=True) for ring in poly.interiors]
            polys.append(Polygon(ext, ints))
        return MultiPolygon(polys)
        
    else:
        raise TypeError(f"Unsupported geometry type: {gtype}. Expected LineString, MultiLineString, Polygon, MultiPolygon, or LinearRing.")