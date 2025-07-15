import warnings
import random
import geopandas as gpd

from shapely import MultiPoint

def kmeans_selection(points, ratio, column):
    """
    Reduce a set of points by selection using K-Means clustering.

    For each K-Means cluster calculated depending on the ratio, the point
    with the largest value is selected.

    For more information about K-Means clustering,
    here is a link to the related `wikipedia article
    <https://en.wikipedia.org/wiki/K-means_clustering>`_

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    ratio : float
        A value between 0 (all points are removed) and 1 (all points are kept).
    column : str
        Name of the column to use.

    Returns
    -------
    points : GeoDataFrame of Point
        The provided set of points with two new attributes:

        - 'selected_kmeans': Set to True if the point has been selected.
        - 'cluster_count': The number of points in the K-Means cluster represented by this point.

    See Also
    --------
    kmeans_simplification :
        Reduce a set of points by simplification using K-Means clustering.
    kmeans_aggregation :
        Reduce a set of points by aggregation using K-Means clustering.
    """
    return __reduce_kmeans(points, ratio, mode='selection', column=column)

def kmeans_simplification(points, ratio):
    """
    Reduce a set of points by simplification using K-Means clustering.

    For each K-Means cluster calculated depending on the ratio, the point
    closest to the cluster's centroid is selected.

    For more information about K-Means clustering,
    here is a link to the related `wikipedia article
    <https://en.wikipedia.org/wiki/K-means_clustering>`_

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    ratio : float
        A value between 0 (all points are removed) and 1 (all points are kept).

    Returns
    -------
    points : GeoDataFrame of Point
        The provided set of points with one new attributes:

        - 'selected_kmeans': Set to True if the point has been selected.

    See Also
    --------
    kmeans_selection :
        Reduce a set of points by selection using K-Means clustering.
    kmeans_aggregation :
        Reduce a set of points by aggregation using K-Means clustering.
    """
    return __reduce_kmeans(points, ratio, mode='simplification')

def kmeans_aggregation(points, ratio, column=None):
    """
    Reduce a set of points by aggregation using K-Means clustering.

    For each K-Means cluster calculated depending on the ratio, the
    centroid is calculated.

    For more information about K-Means clustering,
    here is a link to the related `wikipedia article
    <https://en.wikipedia.org/wiki/K-means_clustering>`_

    Parameters
    -------
    points : GeoDataFrame of Point
        The points to reduce.
    ratio : float
        A value between 0 (all points are removed) and 1 (all points are kept).
    column : str, optional
        Name of a numeric column to calculate the sum and mean.
    
    Returns
    -------
    points : GeoDataFrame of Point
        A new set of points representing the centroids with:

        - 'cluster_count': The number of points in the K-Means cluster represented by this point.
        - 'column_total': The sum of the column values if provided.
        - 'column_mean': The mean of the column values if provided.

    See Also
    --------
    kmeans_selection :
        Reduce a set of points by selection using K-Means clustering.
    kmeans_simplification :
        Reduce a set of points by simplification using K-Means clustering.
    """
    return __reduce_kmeans(points, ratio, mode='aggregation', column=column)
    

# @gtouya
def __reduce_kmeans(points, ratio, mode, column=None):
    """
    Reduce a set of points using K-Means clustering.

    For more information about K-Means clustering,
    here is a link to the related `wikipedia article
    <https://en.wikipedia.org/wiki/K-means_clustering>`_

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    ratio : float
        A value between 0 (all points are removed) and 1 (all points are kept).
    mode : str, optional
        There are three available modes:

        - *'selection'*: inside the cluster, only point with the largest
          value in the chosen column is retained.
          This option requires the column parameter to be provided.
        - *'simplification'*: the point retained in the cluster
          is the closest to the centroid of the cluster.
        - *'aggregation'*: the points are all aggregated to
          the centroid of the cluster. The count of point is added as a new attribute.
          If a column name is provided, also adds the sum of the attribute.
    column : str, optional
        Name of the column to use.

    Returns
    -------
    points : GeoDataFrame of Point
        
        - In *'selection'* and *'simplification'* mode, returns the provided set of points with a
          boolean attribute set to True if the point has been selected.
        - In *'aggregation'* mode, returns a new set of point.
    """
    if column is None and mode == 'selection':
        raise Exception('Provide an attribute name in selection mode.')

    if column is not None and mode == 'simplification':
        warnings.warn("Warning: There is no need to indicate a column name in simplification mode.")

    geometries = list(points.geometry)

    final_pts = []

    # first compute the number of points kept in the generalised point set
    k = int(round(len(geometries) * ratio, 0))

    if k == 0:
        return gpd.GeoDataFrame()
    if k == len(geometries):
        return points

    indexes = [ ip for ip in range(0, len(geometries)) ]
    
    #***********************************************
    # create k cluster with the K-Means algorithm
    clusters = []
    movement_tol = 3
    # initialise the clusters with random centers
    initial_centers = random.sample(geometries, k)
    centers = initial_centers.copy()
    while movement_tol > 2:
        movement_tol = 0
        # build empty clusters
        clusters.clear()
        for i in range(0,k):
            clusters.append([])
        # put the points in the closest cluster
        for ip in indexes:
            mindist = float("inf")
            nearest = -1
            i = 0
            for center in centers:
                dist = center.distance(geometries[ip])
                if(dist < mindist):
                    nearest = i
                    mindist = dist
                i += 1
            clusters[nearest].append(ip)
        # now compute the new centers of the clusters
        for i in range(0, k):
            if len(clusters[i]) == 0:
                continue # do not change the center in this case

            multi_cluster = MultiPoint([ geometries[j] for j in clusters[i] ])
            new_center = multi_cluster.centroid
            # compute the distance between the old center and the new one for this cluster
            dist = new_center.distance(centers[i])
            # check if the distance if bigger than the movement tolerance
            if dist > movement_tol:
                movement_tol = dist
            # apply the new center for this cluster
            centers[i] = new_center

    results = []
    records = points.to_dict('records')

    match mode:
        case 'selection' | 'simplification':
            indexes, values = [], []
            if mode == 'selection':
                for cluster in clusters:
                    # retain the largest value in each cluster for the chosen column
                    selected = None
                    largest = 0
                    for index in cluster:
                        point = records[index]
                        value = point[column]
                        if value > largest:
                            largest = value
                            selected = index

                    if selected is not None:
                        indexes.append(selected)
                        values.append(len(cluster))
                
                for i, r in enumerate(records):
                    selected, value = False, None
                    if i in indexes:
                        selected = True
                        value = values[indexes.index(i)]
                    r['selected_kmeans'] = selected
                    r['cluster_count'] = value
                    results.append(r)

            else:
                for cluster in clusters:
                    # the point retained in the cluster is the closest to the centroid of the cluster
                    # get the centroid of the cluster
                    multi = MultiPoint([ geometries[j] for j in cluster ])
                    center = multi.centroid
                    mindist = float("inf")
                    nearest = None

                    for index in cluster:
                        point = records[index]
                        dist = point['geometry'].distance(center)
                        if dist < mindist:
                            mindist = dist
                            nearest = index

                    if nearest is not None:
                        indexes.append(nearest)

                for i, r in enumerate(records):
                    selected = False
                    if i in indexes:
                        selected = True
                    r['selected_kmeans'] = selected
                    results.append(r)
            
        case 'aggregation':
            for cluster in clusters:
                # the points are all aggregated to the centroid of the cluster.
                multi = MultiPoint([ geometries[j] for j in cluster ])
                center = multi.centroid
                total = 0
                count = 0
                geoms = []
                for index in cluster:
                    point = records[index]
                    if column is not None:
                        total += point[column]
                    count += 1
                    geoms.append(point['geometry'])

                entry = { 'cluster_count': count }
                if column is not None:
                    entry['{0}_total'.format(column)] = total
                    entry['{0}_mean'.format(column)] = total/count
                entry['geometry'] = center
                results.append(entry)

    return gpd.GeoDataFrame(results, crs=points.crs)