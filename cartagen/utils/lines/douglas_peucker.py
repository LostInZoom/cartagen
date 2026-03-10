from shapely import simplify

def douglas_peucker(line, threshold, preserve_topology=True):
    """
    Distance-based line simplification.

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
    line : LineString
        The line to simplify.
    threshold : float
        The distance thresholdto remove the vertex from the line.
    preserve_topology : bool, optional
        If set to True, the algorithm will prevent invalid geometries
        from being created (checking for collapses, ring-intersections, etc).
        The trade-off is computational expensivity.

    Returns
    -------
    LineString

    See Also
    --------
    visvalingam_whyatt :
        Area-based line simplification.
    raposo :
        Hexagon-based line simplification.
    li_openshaw :
        Square grid-based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> douglas_peucker(line, 1.0)
    <LINESTRING (0 0, 2 0, 5 3)>
    """
    return simplify(line, threshold, preserve_topology=preserve_topology)