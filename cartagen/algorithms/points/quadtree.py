import warnings

import geopandas as gpd

from shapely.geometry import Polygon
from cartagen.utils.partitioning.quadtree import PointSetQuadTree

def quadtree_selection(points, depth, column, quadtree=False):
    """
    Reduce a set of points by selection using a quadtree.

    For each quadtree cell, the point with the largest value is selected.

    This algorithm was proposed by Bereuter & Weibel. :footcite:p:`bereuter:2012`
    The quadtree algorithm iteratively divide the point set into four chunks, creating clusters,
    until the depth parameter is reach or only one point remain per cluster.

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    depth : int
        The maximum depth of the quadtree. This acts as a
        selector for the wanted degree of generalisation.
        The lower the value, the more generalised the point set will be.
    column : str
        Name of the column to use.
    quadtree : bool, optional
        If set to True, returns a tuple with the reduced points and the quadtree.
    
    Returns
    -------
    points : GeoDataFrame of Point
        The provided set of points with two new attributes:

        - 'selected_quadtree': Set to True if the point has been selected.
        - 'cell_count': The number of points in the quadtree cell represented by this point.
    
    quadtree : ``QuadTree``, optional

    See Also
    --------
    quadtree_simplification :
        Reduce a set of points by simplification using a quadtree.
    quadtree_aggregation :
        Reduce a set of points by aggregation using a quadtree.

    References
    ----------
    .. footbibliography::
    """
    return __reduce_quadtree(points, depth, mode='selection', column=column, quadtree=quadtree)

def quadtree_simplification(points, depth, quadtree=False):
    """
    Reduce a set of points by simplification using a quadtree.

    For each quadtree cell, the point closest to the cluster's centroid is selected.

    This algorithm was proposed by Bereuter & Weibel. :footcite:p:`bereuter:2012`
    The quadtree algorithm iteratively divide the point set into four chunks, creating clusters,
    until the depth parameter is reach or only one point remain per cluster.

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    depth : int
        The maximum depth of the quadtree. This acts as a
        selector for the wanted degree of generalisation.
        The lower the value, the more generalised the point set will be.
    quadtree : bool, optional
        If set to True, returns a tuple with the reduced points and the quadtree.
    
    Returns
    -------
    points : GeoDataFrame of Point
        The provided set of points with two new attributes:

        - 'selected_kmeans': Set to True if the point has been selected.
    
    quadtree : ``QuadTree``, optional

    See Also
    --------
    quadtree_selection :
        Reduce a set of points by selection using a quadtree.
    quadtree_aggregation :
        Reduce a set of points by aggregation using a quadtree.

    References
    ----------
    .. footbibliography::
    """
    return __reduce_quadtree(points, depth, mode='simplification', column=None, quadtree=quadtree)

def quadtree_aggregation(points, depth, column=None, quadtree=False):
    """
    Reduce a set of points by aggregation using a quadtree.

    For each quadtree cell, the centroid is calculated.

    This algorithm was proposed by Bereuter & Weibel. :footcite:p:`bereuter:2012`
    The quadtree algorithm iteratively divide the point set into four chunks, creating clusters,
    until the depth parameter is reach or only one point remain per cluster.

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    depth : int
        The maximum depth of the quadtree. This acts as a
        selector for the wanted degree of generalisation.
        The lower the value, the more generalised the point set will be.
    column : str, optional
        Name of a numeric column to calculate the sum and mean.
    quadtree : bool, optional
        If set to True, returns a tuple with the reduced points and the quadtree.
    
    Returns
    -------
    points : GeoDataFrame of Point
        A new set of points representing the centroids with:

        - 'cell_count': The number of points in the quadtree cell represented by this point.
        - 'column_total': The sum of the column values if provided.
        - 'column_mean': The mean of the column values if provided.
    
    quadtree : ``QuadTree``, optional

    See Also
    --------
    quadtree_selection :
        Reduce a set of points by selection using a quadtree.
    quadtree_simplification :
        Reduce a set of points by simplification using a quadtree.

    References
    ----------
    .. footbibliography::
    """
    return __reduce_quadtree(points, depth, mode='aggregation', column=column, quadtree=quadtree)

def __reduce_quadtree(points, depth, mode, column=None, quadtree=False):
    """
    Reduce a set of points using a quadtree.

    This algorithm was proposed by Bereuter & Weibel. :footcite:p:`bereuter:2012`
    The quadtree algorithm iteratively divide the point set into four chunks, creating clusters,
    until the depth parameter is reach or only one point remain per cluster.

    Parameters
    ----------
    points : GeoDataFrame of Point
        The point set to reduce.
    depth : int
        The maximum depth of the quadtree. This acts as a
        selector for the wanted degree of generalisation.
        The lower the value, the more generalised the point set will be.
    mode : str, optional
        There are three available modes:

        - *'selection'*: for one cell, the algorithm retains
          the point with the largest value in the chosen column, weighted by the depth of the point.
          This option requires the column parameter to be provided.
        - *'simplification'*: the point retained in the cell
          is the closest to the center of the cell.
        - *'aggregation'*: the points are all aggregated to
          the centroid of the cell. The count of point is added as a new attribute.
          If a column name is provided, also adds the sum of the attribute.

    column : str, optional
        Name of the column to use.
    quadtree : bool, optional
        If set to True, returns a tuple with the reduced points and the quadtree.
    
    Returns
    -------
    GeoDataFrame of Point or tuple (GeoDataFrame of Point, ``QuadTree``)

    References
    ----------
    .. footbibliography::
    """

    if column is not None and mode == 'simplification':
        warnings.warn("Warning: There is no need to indicate a column name in simplification mode.")

    # First get the extent of the quadtree
    xmin, ymin, xmax, ymax = points.geometry.total_bounds
    xcenter, ycenter = xmin + ((xmax - xmin) / 2), ymin + ((ymax - ymin) / 2)
    length = max(xmax - xmin, ymax - ymin) / 2
    xdmin, ydmin, xdmax, ydmax = xcenter - length, ycenter - length, xcenter + length, ycenter + length
    domain = Polygon([(xdmin, ydmin), (xdmin, ydmax), (xdmax, ydmax), (xdmax, ydmin)])

    # Then create the quadtree and populate it with the points
    qtree = PointSetQuadTree(domain, 1)
    qtree.populate(points)

    cells = []
    # get all the quadtree cells at the good depth
    cells_to_process = []
    cells_to_process.append(qtree.sw)
    cells_to_process.append(qtree.nw)
    cells_to_process.append(qtree.se)
    cells_to_process.append(qtree.ne)
    while len(cells_to_process) > 0:
        current_cell = cells_to_process.pop()
        if(current_cell.depth > depth):
            continue
        elif current_cell.depth < depth and current_cell.divided:
            cells_to_process.append(current_cell.sw)
            cells_to_process.append(current_cell.nw)
            cells_to_process.append(current_cell.se)
            cells_to_process.append(current_cell.ne)
        elif current_cell.depth < depth and not current_cell.divided:
            cells.append(current_cell)
            current_cell
        else:
            cells.append(current_cell)

    results = []
    records = points.to_dict('records')
    # then generalise the points based on the chosen mode
    match mode:
        case 'selection':
            if column is None:
                raise Exception('Provide an attribute name in selection mode.')

            indexes, values = [], []
            # loop on the cells
            for cell in cells:
                # get all the points in this cell
                cell_points = cell.get_all_points()
                if len(cell_points) == 0:
                    continue

                # retain the largest value in each cell for the chosen attribute of the point
                selected = None
                largest = 0
                for point, depth in cell_points:
                    value = point[column]
                    if value*depth > largest:
                        largest = value*depth
                        selected = point

                if selected is not None:
                    indexes.append(selected.name)
                    values.append(len(cell_points))
            
            for i, r in enumerate(records):
                selected, value = False, None
                if i in indexes:
                    selected = True
                    value = values[indexes.index(i)]
                r['selected_quadtree'] = selected
                r['cell_count'] = value
                results.append(r)

        case 'simplification':
            indexes = []
            # loop on the cells
            for cell in cells:
                # get all the points in this cell
                cell_points = cell.get_all_points()
                if len(cell_points) == 0:
                    continue

                # the point retained in the cell is the closest to the center of the cell
                center = cell.envelope.centroid
                mindist = float("inf")
                nearest = None
                for point, depth in cell_points:
                    dist = point['geometry'].distance(center)
                    if dist < mindist:
                        mindist = dist
                        nearest = point

                if nearest is not None:
                    indexes.append(nearest.name)
                    # results.append(nearest)
            
            for i, r in enumerate(records):
                selected = False
                if i in indexes:
                    selected = True
                r['selected_quadtree'] = selected
                results.append(r)

        case 'aggregation':
            # loop on the cells
            for cell in cells:
                # get all the points in this cell
                cell_points = cell.get_all_points()
                if len(cell_points) == 0:
                    continue
                
                # the points are all aggregated to the centroid of the cell.
                center = cell.envelope.centroid
                total = 0
                count = 0
                geoms = []
                for point, depth in cell_points:
                    if column is not None:
                        total += point[column]
                    count += 1
                    geoms.append(point['geometry'])

                entry = { 'cell_count': count }
                if column is not None:
                    entry['{0}_total'.format(column)] = total
                    entry['{0}_mean'.format(column)] = total/count
                entry['geometry'] = center
                results.append(entry)
    
    output = gpd.GeoDataFrame(results, crs=points.crs)
    if quadtree:
        return output, qtree
    else:
        return output