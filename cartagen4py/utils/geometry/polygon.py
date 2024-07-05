import shapely, numpy

def polygon_compactness(polygon):
    """
    Calculate the compactness of a polygon.

    This function calculates the compactness of a polygon
    using the Miller index. This index gives a maximum
    value of 1 for circles.

    Parameters
    ----------
    polygon : Polygon
        The polygon to calculated the compactness from.

    Returns
    -------
    float

    Notes
    -----
    The Miller index is calculated using :math:`(4·pi·area)/(perimeter^2)`

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    >>> polygon_compactness(polygon)
    0.7853981633974483
    """
    return 4 * numpy.pi * polygon.area / (polygon.exterior.length * polygon.exterior.length)

def polygon_elongation(polygon):
    """
    Calculate the elongation of a polygon.

    This function calculates the elongation of a polygon
    using the :func:`minimum_rotated_rectangle <shapely.minimum_rotated_rectangle>`.
    It is the ratio between the length and the width of this rectangle.

    Parameters
    ----------
    polygon : Polygon
        The polygon to calculated the elongation from.

    Returns
    -------
    float

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 1), (2, 1), (2, 0), (0, 0)])
    >>> polygon_elongation(polygon)
    2.0
    """
    points = polygon.exterior.coords
    ssr = shapely.MultiPoint(points).minimum_rotated_rectangle.exterior.coords

    x0, x1, x2 = ssr[0][0], ssr[1][0], ssr[2][0]
    y0, y1, y2 = ssr[0][1], ssr[1][1], ssr[2][1]

    length = numpy.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0))
    width = numpy.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))

    if length < width:
        length, width = width, length

    return length / width

def polygon_concavity(polygon):
    """
    Calculate the concavity of a polygon.

    This function calculates the concavity of a polygon as
    its area divided by the area of its convex hull.

    Parameters
    ----------
    polygon : Polygon
        The polygon to calculated the concavity from.

    Returns
    -------
    float

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 2), (1, 1), (2, 2), (2, 0), (0, 0)])
    >>> polygon_concavity(polygon)
    0.75
    """
    hull = shapely.convex_hull(polygon)
    return polygon.area / hull.area