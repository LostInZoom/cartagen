from shapely.ops import transform
from shapely import Point, MultiPoint, LineString, MultiLineString

from cartagen.utils.partitioning.tessellation import HexagonalTessellation

def raposo(line, initial_scale, final_scale, centroid=True, tobler=False):
    """
    Hexagon-based line simplification.
    
    This algorithm proposed by Raposo :footcite:p:`raposo:2013` simplifies lines based on a
    hexagonal tessellation. The algorithm also works for the simplification of the border of a polygon object.
    The idea of the algorithm is to put a hexagonal tessallation on top of the line to simplify,
    the size of the cells depending on the targeted granularity of the line.
    Similarly to the Li-Openshaw algorithm, only one vertex is kept inside each cell.
    This point can be the centroid of the removed vertices, or a projection on the initial line of this centroid.
    The shapes obtained with this algorithm are less sharp than the ones obtained with other algorithms such as Douglas-Peucker.

    The algorithm is dedicated to the smooth simplification of natural features such as rivers, forests, coastlines, lakes.

    Parameters
    ----------
    line : LineString, MultiLineString
        The line to simplify.
    initial_scale : float
        Initial scale of the provided line (25000.0 for 1:25000 scale).
    final_scale : float
        Final scale of the simplified line.
    centroid : bool, optional
        If true, uses the center of the hexagonal cells as the new vertex,
        if false, the center is projected on the nearest point in the initial line.
    tobler : bool, optional
        If True, compute cell resolution based on Tobler’s formula, else uses Raposo's formula

    Returns
    -------
    LineString, MultiLineString

    See Also
    --------
    douglas_peucker :
        Distance-based line simplification.
    visvalingam_whyatt :
        Area-based line simplification.
    li_openshaw :
        Square grid-based line simplification.

    Notes
    -----
    The Tobler based formula to compute hexagonal cell size is :math:`cellsize=5·l·s`
    where :math:`l` is the width of the line in the map in meters
    (*e.g.* 0.0005 for 0.5 mm), and :math:`s` is the target scale denominator.

    Raposo’s formula to compute hexagonal cell size is :math:`cellsize={l/n}·{t/d}`
    where :math:`l` is the length of the line, :math:`n` the number of vertices of the line,
    :math:`t` the denominator of the target scale, and :math:`d` the denominator of the initial scale.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.raposo(line, 5000, 10000)
    <LINESTRING (0 0, 1 0.3333333333333333, 5 3)>
    """

    if line.geom_type not in ['LineString', 'MultiLineString']:
        raise Exception('{0} geometry type cannot be simplified.'.format(line.geom_type))

    if line.geom_type == 'MultiLineString':
        geoms = []
        for geom in list(line.geoms):
            geoms.append(raposo(geom, initial_scale, final_scale, centroid, tobler))
        return MultiLineString(geoms)

    width = 0
    if tobler:
        width = final_scale * 5 / 4000 
    else:
        firstFactor = line.length / (len(line.coords)-1)
        secondFactor = final_scale / initial_scale
        width = firstFactor * secondFactor

    # compute hexagon tessellation
    tessellation = HexagonalTessellation(line.envelope, width)

    current_index = 0
    final_coords = []
    # append the first point of the line, but without the z coordinate
    twod_point = transform(lambda *args: args[:2], Point(line.coords[0]))
    final_coords.append(twod_point.coords[0])
    previous_cell = None
    while (current_index < len(line.coords)-1):
        current_cell = None
        # now loop on the vertices from current index
        # builds a point cloud as a multi-point geometry with all line vertices
        # contained in the current cell.
        point_cloud = []
        for i in range(current_index,len(line.coords)):
            # get the cells containing the point
            point = line.coords[i]
            containing_cells = tessellation.get_containing_cells(point)
            if(current_cell is None):
                if(previous_cell in containing_cells):
                    containing_cells.remove(previous_cell)
                current_cell = containing_cells[0]
                point_cloud.append(point)
                continue

            if (current_cell in containing_cells):
                point_cloud.append(point)
                current_index = i+1
            else:
                current_index = i
                previous_cell = current_cell
                break
        multipoint = MultiPoint(point_cloud)
        if (centroid):
            # replace the points by the centroid of the vertices in the cell
            final_coords.append(multipoint.centroid.coords[0])
        else:
            # find the nearest vertex to the centroid
            nearest = nearest_points(multipoint,multipoint.centroid)
            final_coords.append(nearest[0].coords[0])
        if (current_index == len(line.coords) - 1):
            twod_final_point = transform(lambda *args: args[:2], Point(line.coords[current_index]))
            final_coords.append(twod_final_point.coords[0])
    # add the final point if it is not in the line
    twod_final_point = transform(lambda *args: args[:2], Point(line.coords[len(line.coords) - 1]))
    if(twod_final_point.coords[0] not in final_coords):
        final_coords.append(twod_final_point.coords[0])
    
    # print(final_coords)
    return LineString(final_coords)