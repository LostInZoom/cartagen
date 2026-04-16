from cartagen.utils.lines.smoothing.wma import smooth_wma
from cartagen.utils.lines.simplification.angular import simplify_angular

def simplify_topographic(geometry, iterations=2, angle=10.0, weights=[1.0, 2.0, 1.0]):
    """
    Simplify a line or polygon and mimic hand-made cartographic generalization.

    This algorithm combines both a low-pass filter using weighted moving average (WMA), and
    a high-pass filter using an angular threshold. It is suitable to smooth topographic features
    such as rivers, coastlines, lakes, *etc.* and simplify it to mimic hand-made generalization.
    
    Accept Multi geometries. If a polygon is provided, it also applies the 
    smoothing to its holes using the same parameters.
    
    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to simplify.
        If an open line is provided, the endpoints are preserved.
        If a closed ring or polygon is provided, the simplification wraps around.
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

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 0.9), (2, 0.1), (3, 1), (4, 0)])
    >>> simplify_topographic(line, 1, 15.0)
    <LINESTRING (0 0, 2 0.5, 4 0)>
    """
    smoothed = smooth_wma(geometry, iterations=iterations, weights=weights)
    final = simplify_angular(smoothed, angle=angle)
    return final