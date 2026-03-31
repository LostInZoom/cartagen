import numpy as np
from shapely.geometry import LineString, Polygon, MultiLineString, MultiPolygon, LinearRing

def smooth_chaikin(geometry, iterations=3, keep_ends=True):
    """
    Smooth a line or polygon by cutting corners.
    
    This algorithm was proposed by Chaikin :footcite:p:`chaikin:1974` and this version
    makes use of the equivalent multi-step algorithm introduced by Wu *et al.*
    :footcite:p:`wu:2004`.
    It is a simple subdivision scheme that repeatedly cuts corners of a line or polygon.
    Each iteration doubles the number of vertices and produces a smoother curve that
    converges to a quadratic B-spline.
    
    Accept Multi geometries. If a polygon is provided, it also applies the 
    smoothing to its holes using the same parameters.

    This implementation is a translation of the C++ implementation available
    `here <https://github.com/philipschall/shapelysmooth>`_.
    
    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The line or polygon to smooth.
    iterations : int, optional
        Number of subdivision iterations (k parameter). Default is 3.
        Each iteration doubles the number of vertices (2^k subdivisions per segment).
    keep_ends : bool, optional
        Whether to keep the original endpoints fixed for open lines. Default is True.
        For closed geometries, this parameter is ignored and corners are always cut.
    
    Returns
    -------
    LineString, Polygon, MultiLineString, MultiPolygon, LinearRing
        Smoothed geometry of the same type as input.

    Warning
    -------
    Note that the algorithm (roughly) doubles the amount of nodes at each iteration,
    therefore care should be taken when selecting the number of iterations.

    See Also
    --------
    smooth_catmull_rom :
        Smooth a line or polygon and preserve vertexes.
    smooth_gaussian :
        Smooth a line or a polygon and attenuate its inflexion points.
    smooth_platre :
        Smooth a line and preserve the integrity of sharp turns.
    smooth_taubin :
        Smooth a line or polygon and prevent shrinkage.

    References
    ----------
    .. footbibliography::
    
    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (3, 1)])
    >>> smooth_chaikin(line, iterations=2, keep_ends=True)
    <LINESTRING (0 0, 0.3 0.3, ... 3 1)>
    
    >>> polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    >>> smooth_chaikin(polygon, iterations=2)
    <POLYGON ((0.1 0.2, ... 0.1 0.2))>
    """
    # Handle different geometry types
    if isinstance(geometry, MultiLineString):
        return MultiLineString([
            smooth_chaikin(line, iterations, keep_ends)
            for line in geometry.geoms
        ])
    
    if isinstance(geometry, MultiPolygon):
        return MultiPolygon([
            smooth_chaikin(poly, iterations, keep_ends)
            for poly in geometry.geoms
        ])
    
    if isinstance(geometry, Polygon):
        # Smooth exterior ring
        exterior = smooth_chaikin(geometry.exterior, iterations, keep_ends)
        # Smooth interior rings (holes)
        interiors = [
            smooth_chaikin(interior, iterations, keep_ends)
            for interior in geometry.interiors
        ]
        return Polygon(exterior, interiors)
    
    # Process LineString or LinearRing
    coords = list(geometry.coords)
    is_closed = np.array_equal(coords[0], coords[-1])
    
    # Add extra point for closed geometries
    if is_closed:
        coords.append(coords[1])
    
    size = len(coords)
    exp = 2 ** iterations
    
    # Precompute weight coefficients
    F = np.zeros(exp)
    G = np.zeros(exp)
    H = np.zeros(exp)
    
    for j in range(exp):
        F[j] = 0.5 - 0.5 / exp - j * (1.0 / exp - (j + 1.0) * 0.5 / exp / exp)
        G[j] = 0.5 + 0.5 / exp + j * (1.0 / exp - (j + 1.0) / exp / exp)
        H[j] = j * (j + 1.0) * 0.5 / exp / exp
    
    # Build new linestring
    new_coords = [coords[0]]
    
    for i in range(1, size - 1):
        p_prev = np.array(coords[i - 1])
        p_curr = np.array(coords[i])
        p_next = np.array(coords[i + 1])
        
        for j in range(exp):
            # Weighted sum of three consecutive points
            new_point = F[j] * p_prev + G[j] * p_curr + H[j] * p_next
            new_coords.append(tuple(new_point))
    
    # Handle endpoints
    if keep_ends and not is_closed:
        new_coords.append(coords[-1])
    else:
        # Adjust first point
        p0 = np.array(coords[0])
        p1 = np.array(coords[1])
        new_coords[0] = tuple(
            0.5 * (1.0 + 1.0 / exp) * p0 + 0.5 * (1.0 - 1.0 / exp) * p1
        )
        
        if is_closed:
            # Close the ring
            if isinstance(geometry, LinearRing):
                return LinearRing(new_coords)
            return LineString(new_coords)
        
        # Adjust last point
        p_n2 = np.array(coords[size - 2])
        p_n1 = np.array(coords[size - 1])
        new_coords.append(tuple(
            0.5 * (1.0 - 1.0 / exp) * p_n2 + 0.5 * (1.0 + 1.0 / exp) * p_n1
        ))
    
    # Return appropriate geometry type
    if isinstance(geometry, LinearRing):
        return LinearRing(new_coords)
    return LineString(new_coords)