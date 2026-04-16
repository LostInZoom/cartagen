from shapely import simplify

def simplify_douglas_peucker(geometry, threshold, preserve_topology=True):
    """
    Simplify a line or polygon using a distance-based selection.

    This algorithm was proposed by Ramer :footcite:p:`ramer:1972` and by Douglas and Peucker.
    :footcite:p:`douglas:1973` It is a line filtering algorithm, which means that it
    filters the vertices of the line (or polygon) to only retain the most important ones
    to preserve the shape of the line. The algorithm iteratively searches the most
    characteristics vertices of portions of the line and decides to retain
    or remove them given a distance threshold.

    The algorithm tends to unsmooth geographic lines, and is rarely used to simplify geographic features.
    But it can be very useful to quickly filter the vertices of a line inside another algorithm.

    This is a simple wrapper around :func:`shapely.simplify() <shapely.simplify()>`.

    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to simplify.
    threshold : float
        The distance threshold to remove the vertex from the line.
    preserve_topology : bool, optional
        If set to True, the algorithm will prevent invalid geometries
        from being created (checking for collapses, ring-intersections, etc).
        The trade-off is computational expensivity.

    Returns
    -------
    LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

    See Also
    --------
    simplify_angular :
        Simplify a line or polygon by removing vertexes with small angles.
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
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> simplify_douglas_peucker(line, 1.0)
    <LINESTRING (0 0, 2 0, 5 3)>
    """
    return simplify(geometry, threshold, preserve_topology=preserve_topology)