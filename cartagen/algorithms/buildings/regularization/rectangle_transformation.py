import shapely

from cartagen.utils.geometry.polygon import enclosing_rectangle

def regularize_building_rectangle(polygon, factor=1.0, method='mbr'):
    """
    Transform a polygon into a rectangle.

    This function transforms a polygon to a rectangle using
    the minimum rotated rectangle and scale it up or down.

    Parameters
    ----------
    polygon : Polygon
        The polygon to regularize.
    factor : float, optional
        The scaling factor to apply.
    method : str, optional
        The method to calculate the rectangle:

        - **'mbr'** calculate the minimum rotated bounding rectangle.
        - **'mbtr'** calculate minimum rotated bounding touching rectangle.
          It is the same as the mbr but the rectangle and the polygon
          must have at least one side in common.

    Returns
    -------
    Polygon

    See Also
    --------
    regularize_building_regression :
        Regularize a polygon using recursive linear regression.
    regularize_building_fer :
        Regularize a polygon using feature edge reconstruction.

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 2), (1, 2), (1, 1), (2, 1), (2, 0), (0, 0)])
    >>> regularize_building_rectangle(polygon)
    <POLYGON ((2 0, 2 2, 0 2, 0 0, 2 0))>
    """
    if method == 'mbr':
        mbr = enclosing_rectangle(polygon, mode='hull')
    elif method == 'mbtr':
        mbr = enclosing_rectangle(polygon, mode='input')
    else:
        raise Exception('Selected method does not exists: {0}'.format(method))

    return shapely.affinity.scale(mbr, xfact=factor, yfact=factor, origin=mbr.centroid)