import numpy as np
from shapely.geometry import LineString, Polygon, MultiLineString, MultiPolygon, LinearRing

def smooth_catmull_rom(geometry, subdivisions=10, alpha=0.5,):
    """
    Smooth a line or polygon and preserve vertexes.
    
    This algorithm was proposed by Catmull and Rom :footcite:p:`catmull:1974`, this is the
    version proposed by Barry and Goldman :footcite:p:`barry:1988` that makes use of the
    de Boor's algorithm for evaluating spline curves in B-spline form.

    Catmull-Rom splines are a family of cubic interpolating splines that pass 
    through all control points. The alpha parameter controls the type of 
    parameterization: 0 for uniform, 0.5 for centripetal (default), and 1 for 
    chordal. Centripetal parameterization avoids cusps and self-intersections.
    
    Accept Multi geometries. If a polygon is provided, it also applies the 
    smoothing to its holes using the same parameters.

    This implementation is a translation of the C++ implementation available
    `here <https://github.com/philipschall/shapelysmooth>`_.
    
    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The line or polygon to smooth.
        The spline passes through all original vertices.
    subdivisions : int, optional
        Number of interpolated points between each pair of control points.
        Default is 10. Higher values produce smoother curves.
    alpha : float, optional
        Parameterization type. Default is 0.5 (centripetal):

        - 0.0: uniform parameterization
        - 0.5: centripetal parameterization (recommended, prevents loops)
        - 1.0: chordal parameterization

    Returns
    -------
    LineString, Polygon, MultiLineString, MultiPolygon, LinearRing
        Smoothed geometry of the same type as input.

    See Also
    --------
    smooth_chaikin :
        Smooth a line or polygon by cutting corners.
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
    >>> smooth_catmull_rom(line, alpha=0.5, subdivisions=10)
    <LINESTRING (0 0, 0.1 0.12, ... 3 1)>
    
    >>> polygon = Polygon([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)])
    >>> smooth_catmull_rom(polygon, alpha=0.5, subdivisions=20)
    <POLYGON ((0 0.1, ... 0 0.1))>
    """
    # Handle different geometry types
    if isinstance(geometry, MultiLineString):
        return MultiLineString([
            smooth_catmull_rom(line, alpha, subdivisions)
            for line in geometry.geoms
        ])
    
    if isinstance(geometry, MultiPolygon):
        return MultiPolygon([
            smooth_catmull_rom(poly, alpha, subdivisions)
            for poly in geometry.geoms
        ])
    
    if isinstance(geometry, Polygon):
        # Smooth exterior ring
        exterior = smooth_catmull_rom(geometry.exterior, alpha, subdivisions)
        # Smooth interior rings (holes)
        interiors = [
            smooth_catmull_rom(interior, alpha, subdivisions)
            for interior in geometry.interiors
        ]
        return Polygon(exterior, interiors)
    
    # Process LineString or LinearRing
    coords = list(geometry.coords)
    is_closed = np.array_equal(coords[0], coords[-1])
    size = len(coords)
    
    # Add control points for boundary conditions
    if is_closed:
        # For closed curves, wrap around
        coords.insert(0, coords[-2])
        coords.append(coords[2])
    else:
        # For open curves, extrapolate endpoints
        p0 = np.array(coords[0])
        p1 = np.array(coords[1])
        p_n1 = np.array(coords[-1])
        p_n2 = np.array(coords[-2])
        
        # Extrapolate first control point
        coords.insert(0, tuple(2 * p0 - p1))
        # Extrapolate last control point
        coords.append(tuple(2 * p_n1 - p_n2))
    
    # Build new linestring
    new_coords = [coords[1]]
    
    # Process each segment
    for k in range(len(coords) - 3):
        # Extract four control points for this segment
        slice_points = [np.array(coords[k + i]) for i in range(4)]
        
        # Compute knot vector (tangent parameters)
        tangents = [0.0]
        for j in range(3):
            dist = np.linalg.norm(slice_points[j + 1] - slice_points[j])
            tangents.append(tangents[-1] + dist ** alpha)
        
        # Generate interpolation parameters
        t_start = tangents[1]
        t_end = tangents[2]
        seg_length = (t_end - t_start) / subdivisions
        t_values = [t_start + i * seg_length for i in range(1, subdivisions)]
        
        # Evaluate spline at each parameter value
        for t in t_values:
            # Recursive evaluation (de Boor algorithm)
            point = _recursive_catmull_rom_eval(slice_points, tangents, t)
            new_coords.append(tuple(point))
        
        # Add the control point
        new_coords.append(tuple(slice_points[2]))
    
    # Return appropriate geometry type
    if isinstance(geometry, LinearRing):
        return LinearRing(new_coords)
    return LineString(new_coords)


def _recursive_catmull_rom_eval(points, tangents, t):
    """
    Evaluate Catmull-Rom spline at parameter t using recursive algorithm.
    
    This implements the de Boor algorithm for B-spline evaluation.
    
    Parameters
    ----------
    points : list of np.ndarray
        Four control points for the spline segment.
    tangents : list of float
        Knot vector with four values.
    t : float
        Parameter value at which to evaluate the spline.
    
    Returns
    -------
    np.ndarray
        Interpolated point coordinates.
    """
    # Initialize pyramid of points
    pyramid = [points.copy()]
    
    # Recursive evaluation
    for r in range(1, 4):
        level = []
        idx = max(r - 2, 0)
        
        for i in range(4 - r):
            # Compute blend weights
            t_left = tangents[i + r - idx]
            t_right = tangents[i + idx]
            denominator = t_left - t_right
            
            if abs(denominator) < 1e-10:
                # Handle degenerate case
                left_weight = 0.5
                right_weight = 0.5
            else:
                left_weight = (t_left - t) / denominator
                right_weight = (t - t_right) / denominator
            
            # Linear interpolation between previous level points
            new_point = (
                left_weight * pyramid[r - 1][i] + 
                right_weight * pyramid[r - 1][i + 1]
            )
            level.append(new_point)
        
        pyramid.append(level)
    
    # Return the final interpolated point
    return pyramid[-1][0]