from cartagen.utils.lines.smoothing.wma import smooth_wma
from cartagen.utils.lines.simplification.angular import simplify_angular

def smooth_topographic(geometry, iterations=2, angle=10.0, weights=[1.0, 2.0, 1.0]):
    """
    Smooth a line or polygon and mimic hand-made cartographic generalization.

    This algorithm described in Müller and Wang :footcite:p:`muller:1992`
    combines both a low-pass filter using weighted moving average (WMA), and
    a high-pass filter using an angular threshold. It is suitable to smooth topographic features
    such as rivers, coastlines, lakes, *etc.* and simplify it to mimic hand-made generalization.
    
    Accept Multi geometries. If a polygon is provided, it also applies the 
    smoothing to its holes using the same parameters.
    
    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to smooth.
        If an open line is provided, the endpoints are preserved.
        If a closed ring or polygon is provided, the smoothing wraps around.
    iterations : int, optional
        Number of low-pass filter passes before thinning. Default is 2.
    angle : float, optional
        Turning-angle threshold in degrees for the high-pass step. Default is 10.0.
    weights : sequence of odd length, optional
        Weights for the moving average window. Default is `[1, 2, 1]`.
    
    Returns
    -------
    LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        Simplified geometry of the same type as input.

    See Also
    --------
    smooth_gaussian :
        Smooth a line or a polygon and attenuate its inflexion points.
    smooth_platre :
        Smooth a line and preserve the integrity of sharp turns.
    smooth_snake :
        Smooth a line or polygon using the snake method.
    smooth_taubin :
        Smooth a line or polygon and prevent shrinkage.
    smooth_wma :
        Smooth a line or polygon using a low-pass filter.

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 0.9), (2, 0.1), (3, 1), (4, 0)])
    >>> simplify_topographic(line, 1, 15.0)
    <LINESTRING (0 0, 2 0.5, 4 0)>
    """
    smoothed = smooth_wma(geometry, iterations=iterations, weights=weights)
    final = simplify_angular(smoothed, angle=angle)
    return final