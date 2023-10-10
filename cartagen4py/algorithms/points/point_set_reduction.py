# This file contains point set reduction algorithms, following the classification of point set operators from Bereuter & Weibel (2013).
# A reduction means that there is still a set of points as output of the algorithm, but fewer points.

from shapely.geometry import MultiPoint, Point, Polygon
from cartagen4py.algorithms.points.point_set_quadtree import PointSetQuadTree
import random

# A point set reduction algorithm based on a K-Means clustering. The reduction is based an a shrink ratio parameter between 0 (all points are removed) and 1 (all points are kept).
# Two options are possible: either keeping one of the initial points to replace a cluster (default option) or replace the cluster by its centroid.
# @gtouya
def kmeans_point_set_reduction(points, shrink_ratio, centroid_option = False):
    final_pts = []

    # first compute the number of points kept in the generalised point set
    k = int(round(len(points) * shrink_ratio, 0))

    if k == 0:
        return final_pts
    if k == len(points):
        return points
    
    #***********************************************
    # create k cluster with the K-Means algorithm
    clusters = []
    movement_tol = 3
    # initialise the clusters with random centers
    initial_centers = random.sample(points, k)
    centers = initial_centers.copy()
    while(movement_tol > 2):
        movement_tol = 0
        # build empty clusters
        clusters.clear()
        for i in range(0,k):
            clusters.append([])
        # put the points in the closest cluster
        for point in points:
            mindist = float("inf")
            nearest = -1
            i = 0
            for center in centers:
                dist = center.distance(point)
                if(dist < mindist):
                    nearest = i
                    mindist = dist
                i += 1
            clusters[nearest].append(point)
        # now compute the new centers of the clusters
        for i in range(0, k):
            if len(clusters[i]) == 0:
                continue # do not change the center in this case
            multi_cluster = MultiPoint(clusters[i])
            new_center = multi_cluster.centroid
            # compute the distance between the old center and the new one for this cluster
            dist = new_center.distance(centers[i])
            # check if the distance if bigger than the movement tolerance
            if dist > movement_tol:
                movement_tol = dist
            # apply the new center for this cluster
            centers[i] = new_center

    #***********************************************
    # then collapse each cluster into one point
    for cluster in clusters:
        # get the centroid of the cluster
        multi = MultiPoint(cluster)
        center = multi.centroid
        if(centroid_option):
            # append center to the list of final points
            final_pts.append(center)
        else:
            # get the nearest point from the center of the cluster
            nearest = center
            mindist = float("inf")
            for point in cluster:
                dist = center.distance(point)
                if dist < mindist:
                    nearest = point
                    mindist = dist
            final_pts.append(nearest)

    return final_pts

def quadtree_point_set_reduction(points, depth, mode='simplification', attribute = ""):
    """Algorithm to reduce a point set based on a quadtree. The algorithm was proposed by Bereuter & Weibel (2012).
    
    The points should be a GeoDataframe from Geopandas with a Point geometry. The 'depth' parameter is the degree of generalisation as it retains one point per cell when the quadtree is sliced at the given depth.
    
    There are different modes for this algorithm:
    - mode = 'selection' means that for one cell, the algorithm retains the point with the largest value in the chosen attribute, weighted by the depth of the point. 
    - mode = 'simplification' means that the point retained in the cell is the closest to the center of the cell
    - mode = 'aggregation' means that the points are all aggregated to the centroid of the points.

    The algorithm returns the list of tuples (geometry, index, nb_of_points) where index is the index of the point in the initial Geodataframe (-1 if the point was created), and nb_of_points gives the amount of initial points replaced (which can be used to weight the size of the symbol of this point). 
    """

    # First get the extent of the quadtree

    # Then create the quadtree and populate it with the points
    xmin, ymin, xmax, ymax = points.geometry.total_bounds
    length = max(xmax-xmin, ymax-ymin)
    domain = Polygon([(xmin, ymin), (xmin + length, ymin), (xmin + length, ymin + length), (xmin, ymin + length), (xmin,ymin)])
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

    # loop on the cells
    output = []
    for cell in cells:
        # get all the points in this cell
        cell_points = cell.get_all_points()
        if len(cell_points) == 0:
            continue

        # then generalise the points based on the chosen mode
        match mode:
            case 'selection':
                # retain the largest value in each cell for the chosen attribute of the point
                selected = None
                largest = 0
                for point, depth in cell_points:
                    value = point[attribute]
                    if value*depth > largest:
                        largest = value*depth
                        selected = point
                output.append((selected['geometry'],selected.index,len(cell_points)))

            case 'simplification':
                # the point retained in the cell is the closest to the center of the cell
                center = cell.envelope.centroid
                mindist = float("inf")
                nearest = None
                for point in cell_points:
                    dist = point[0]['geometry'].distance(center)
                    if dist < mindist:
                        mindist = dist
                        nearest = point
                output.append((nearest[0]['geometry'], nearest.index, len(cell_points)))

            case 'aggregation':
                # the points are all aggregated to the centroid of the points.
                geoms = []
                for point in cell_points:
                    geoms.append(point[0]['geometry'])
                multi = MultiPoint(geoms)
                centroid = multi.centroid
                output.append((centroid, -1,len(cell_points)))
    
    return output, qtree


