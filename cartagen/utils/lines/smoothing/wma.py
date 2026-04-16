import numpy as np
import math
from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, LinearRing, GeometryCollection

def smooth_wma(geometry, iterations=1, weights=[1.0, 2.0, 1.0]):
    """
    Smooth a line or polygon using a low-pass filter.

    This algorithm applies a weighted moving average (WMA) window over the vertices.
    The default weights `[1, 2, 1]` heavily weight the central point to suppress
    unwanted short-wavelength noise while minimizing overall geometric shrinkage.
    
    Accept Multi geometries. If a polygon is provided, it also applies the 
    smoothing to its holes using the same parameters.
    
    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to smooth.
        If an open line is provided, the endpoints are preserved.
        If a closed ring or polygon is provided, the smoothing wraps around.
    iterations : int, optional
        Number of iterations to apply the filter. Default is 1.
    weights : sequence of odd length, optional
        Weights for the moving average window. Default is `[1, 2, 1]`.
    
    Returns
    -------
    LineString, Polygon, MultiLineString, MultiPolygon, LinearRing
        Smoothed geometry of the same type as input.

    See Also
    --------
    smooth_catmull_rom :
        Smooth a line or polygon and preserve vertexes.
    smooth_chaikin :
        Smooth a line or polygon by cutting corners.
    smooth_gaussian :
        Smooth a line or a polygon and attenuate its inflexion points.
    smooth_platre :
        Smooth a line and preserve the integrity of sharp turns.
    smooth_taubin :
        Smooth a line or polygon and prevent shrinkage.

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (3, 1), (4, 0)])
    >>> smooth_wma(line, 1)
    <LINESTRING (0 0, 1 0.5, 2 0.5, 3 0.5, 4 0)>
    """
    if geometry.is_empty or geometry.geom_type in ['Point', 'MultiPoint']:
        return geometry

    # --- Recursive Handling ---
    if geometry.geom_type in ['MultiLineString', 'MultiPolygon', 'GeometryCollection']:
        parts = [smooth_wma(part, iterations, weights) for part in geometry.geoms]
        if geometry.geom_type == 'MultiLineString':
            return MultiLineString(parts)
        elif geometry.geom_type == 'MultiPolygon':
            return MultiPolygon(parts)
        else:
            return GeometryCollection(parts)

    if geometry.geom_type == 'Polygon':
        ext = smooth_wma(geometry.exterior, iterations, weights)
        ints = []
        for ring in geometry.interiors:
            new_ring = smooth_wma(ring, iterations, weights)
            if len(new_ring.coords) >= 4:
                ints.append(new_ring)
        if len(ext.coords) < 4:
            return geometry
        poly = Polygon(ext, ints)
        return poly if poly.is_valid else poly.buffer(0)

    # --- Core Linear Logic ---
    coords = np.array(geometry.coords)
    is_closed = geometry.geom_type == 'LinearRing' or (
        len(coords) > 0 and np.allclose(coords[0], coords[-1])
    )
    
    if len(coords) < 3:
        return geometry

    # Strip closing duplicate for processing
    if is_closed:
        coords = coords[:-1]

    # Save endpoints for open lines
    start_pt, end_pt = coords[0].copy(), coords[-1].copy()

    # WMA Math
    w = np.array(weights, dtype=float)
    w /= w.sum()
    half = len(w) // 2
    result = coords.copy()

    for _ in range(iterations):
        if is_closed:
            padded = np.pad(result, ((half, half), (0, 0)), mode='wrap')
        else:
            padded = np.pad(result, ((half, half), (0, 0)), mode='edge')
            
        smoothed = np.zeros_like(result)
        for i in range(len(result)):
            window = padded[i : i + len(w)]
            smoothed[i] = np.sum(window * w[:, np.newaxis], axis=0)
        
        result = smoothed
        if not is_closed:
            result[0], result[-1] = start_pt, end_pt
    
    # Re-close and Re-construct
    if is_closed:
        result = np.vstack([result, result[0]])
        if geometry.geom_type == 'LinearRing':
            return LinearRing(result) if len(result) >= 4 else geometry
    
    return LineString(result)