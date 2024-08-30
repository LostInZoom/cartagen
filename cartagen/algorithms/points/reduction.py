# This file contains point set reduction algorithms, following the classification of point set operators from Bereuter & Weibel (2013).
# A reduction means that there is still a set of points as output of the algorithm, but fewer points.

from shapely.geometry import MultiPoint, Point, Polygon
from cartagen.utils.partitioning.quadtree import PointSetQuadTree
import random
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

# @gtouya
def reduce_kmeans(points, shrink_ratio, centroid=False):
    """
    Reduce a set of points using K-Means clustering.

    For more information about K-Means clustering,
    here is a link to the related `wikipedia article
    <https://en.wikipedia.org/wiki/K-means_clustering>`_

    Parameters
    ----------
    points : list of Point
        The points to reduce.
    shrink_ratio : float
        A value between 0 (all points are removed) and 1 (all points are kept).
    centroid : bool, optional
        If set to True, replace the cluster by its centroid.
        If False, replace the cluster by the point closest to the centroid of the cluster.

    Returns
    -------
    list of Point

    See Also
    --------
    reduce_quadtree :
        Reduce a set of points using a quadtree.
    reduce_labelgrid :
        Reduce a set of points using the Label Grid method.
    """

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
        if(centroid):
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

def reduce_quadtree(points, depth, mode='simplification', attribute=None):
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

        - *'selection'* means that for one cell, the algorithm retains
          the point with the largest value in the chosen attribute,
          weighted by the depth of the point. This option requires
          the attribute parameter to be provided.
        - *'simplification'* means that the point retained in the cell
          is the closest to the center of the cell.
        - *'aggregation'* means the points are all aggregated to
          the centroid of the points.
    
    Returns
    -------
    reduced : list of tuple
        The reduced points as tuples composed of three elements:

        #. The geometry of the reduced point.
        #. The index of the point in the initial Geodataframe (-1 if the point was created).
        #. The amount of initial points replaced (which can be used to weight the size of the symbol of this point).
    
    quadtree : QuadTree
        The quadtree object.

    See Also
    --------
    reduce_kmeans :
        Reduce a set of points using K-Means clustering.
    reduce_labelgrid :
        Reduce a set of points using the Label Grid method.

    References
    ----------
    .. footbibliography::
    """

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
                if attribute is None:
                    raise Exception('Provide an attribute name in selection mode.')
                # retain the largest value in each cell for the chosen attribute of the point
                selected = None
                largest = 0
                for point, depth in cell_points:
                    value = point[attribute]
                    if value*depth > largest:
                        largest = value*depth
                        selected = point

                if selected is not None:
                    output.append((selected['geometry'], selected.name, len(cell_points)))

            case 'simplification':
                # the point retained in the cell is the closest to the center of the cell
                center = cell.envelope.centroid
                mindist = float("inf")
                nearest = None
                for point, depth in cell_points:
                    dist = point['geometry'].distance(center)
                    if dist < mindist:
                        mindist = dist
                        nearest = point
                output.append((nearest['geometry'], nearest.name, len(cell_points)))

            case 'aggregation':
                # the points are all aggregated to the centroid of the points.
                geoms = []
                for point in cell_points:
                    geoms.append(point[0]['geometry'])
                multi = MultiPoint(geoms)
                centroid = multi.centroid
                output.append((centroid, -1, len(cell_points)))
    
    return output, qtree

def reduce_labelgrid(points, attribute, width, height, shape='square', mode='selection', grid=False):
    """
    Reduce a set of points using the Label Grid method.

    This algorithm was proposed by Gröbe & Burghardt. :footcite:p:`grobe:2021`
    It assigns each point to a grid cell of a given height, width and shape
    before reducing the points per cells either by selection or aggregation.
    
    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    attribute : str
        Name of the attribute to rank the points.
    width : float
        Width of the grid cells.
    height : float
        Height of the grid cells.
    shape : str, optional
        Shape of the grid cells, can be 'square', 'diamond', 'hexagonal'.
    mode : str, optional
        The reduction method used, can be:

        - **'selection'** keeps the point with the highest attribute value.
        - **'aggregation'** keeps the centroid of the cell and adds the number of points within that cell
          as an attribute.
    grid : bool, optional
        If set to True, returns a tuple with the points and the grid.

    Returns
    -------
    GeoDataFrame of Point or tuple (GeoDataFrame of Point, GeoDataFrame of Polygon)

    See Also
    --------
    reduce_kmeans :
        Reduce a set of points using K-Means clustering.
    reduce_quadtree :
        Reduce a set of points using a quadtree.

    References
    ----------
    .. footbibliography::
    """
    lg = LabelGrid(points, attribute, width, height, shape, mode)
    lg.set_point_label_grid()
    result = lg.getPointResults()
    if grid:
        rgrid = lg.getGrid()
        return result, rgrid
    else:
        return result

"""
Created on Sun Jun 16 22:40:42 2024.

@author: vpech
"""
class LabelGrid():
    """
    A class for label grid algorithm.
    
    The algorithm was proposed by Gröbe, M. and Burghardt, D., 2021. 
    Scale-Dependent Point Selection Methods for Web Maps. KN - Journal of 
    Cartography and Geographic Information.
    """
    
    def __init__(self, points, attribute, width, height, typ='square', mode='selection'):
        """
        Construct all the necessary attributes for the LabelGrid object.

        Parameters
        ----------
        points : GeoDataFrame
            GeoDataFrame of object points we want to apply the label grid algorithm.
        width : int
            Width of a cell.
        height : int
            Height of a cell.
        typ : str
            Form of a cell. The default is 'carre'.
        """
        self.__attribute = attribute
        self.__width = width
        self.__height = height
        self.__typ = typ
        self.__mode = mode
        self.__points = points
        self.__points_res = None
        self.__grid = None

    def __create_hexagonalCell(self, coords_min, coords_max):
        """
        Create a grid of hexagonal cells.

        Parameters
        ----------
        coords_min : tuple
            Coordinates of the bottom left cell first point.
        coords_max : tuple
            Coordinates of the top right cell last point.

        Returns
        -------
        polygons : GeoDataFrame
            Contains hexagonal polygon objects.
        """
        xstep = self.__width*3
        ystep = self.__width*2
        cols = list(np.arange(coords_min[0], coords_max[0]+xstep, xstep))
        rows = list(np.arange(coords_min[1], coords_max[1]+ystep, ystep))
        
        polygons = [Polygon([(x, y), 
                             (x + self.__width, y),
                             (x+self.__width*3/2, y+self.__width),
                             (x + self.__width, y+2*self.__width),
                             (x, y+2*self.__width),
                             (x-self.__width/2, y+self.__width)]) 
                    for x in cols[:-1] for y in rows[:-1]]
        return polygons

    def __createGrid(self):
        """
        Create the grid of the label grid.

        Returns
        -------
        GeoDataFrame
            Contains the grid polygons.

        """
        xmin, ymin, xmax, ymax = self.__points.total_bounds
        
        if self.__typ == 'square':
            cols = list(np.arange(xmin, xmax + self.__width, self.__width))
            rows = list(np.arange(ymin, ymax + self.__height, self.__height))
        
            polygons = [Polygon([(x, y), 
                                 (x + self.__width, y), 
                                 (x + self.__width, y + self.__height), 
                                 (x, y + self.__height)]) 
                        for x in cols[:-1] for y in rows[:-1]]
            
        if self.__typ == 'diamond':
            cols = list(np.arange(xmin, xmax + self.__width*2, self.__width))
            rows = list(np.arange(ymin-self.__width, ymax + self.__width, self.__width))
            
            polygons = [Polygon([(x, y), 
                                 (x + self.__width / 2, y + self.__width / 2), 
                                 (x, y + self.__width), 
                                 (x - self.__width / 2, y + self.__width / 2)]) 
                        for x in cols[:-1] for y in rows[:-1]]
            
        if self.__typ == 'hexagonal':

            odd_coords_min = (xmin-self.__width, ymin-self.__width)
            odd_coords_max = (xmax+self.__width, ymax)
            
            odd_poly = self.__create_hexagonalCell(odd_coords_min, odd_coords_max)
            
            not_odd_coords_min = (odd_coords_min[0]+3/2*self.__width, odd_coords_min[1]-self.__width)
            not_odd_coords_max = (odd_coords_max[0], odd_coords_max[1]+self.__width)
            
            not_odd_poly = self.__create_hexagonalCell(not_odd_coords_min, not_odd_coords_max)
            
            polygons = odd_poly + not_odd_poly
            
        return gpd.GeoDataFrame({'geometry':polygons})

    def set_point_label_grid(self):
        """Set the attributes points_res and grid of the LabelGrid object."""
        # self.__points.reset_index(drop=True, inplace=True)
        # self.__points = self.__points[['geometry', self.__attribute]]
        
        self.__grid = self.__createGrid()
        
        lst_in_value = [self.__points.loc[cell.contains(self.__points['geometry']), self.__attribute] 
                        for cell in self.__grid['geometry']]

        lst_inter_value = [self.__points.loc[cell.touches(self.__points['geometry']), self.__attribute] 
                           for cell in self.__grid['geometry']]
        
        lst_all_value = [lst_inter_value[i] 
                         if lst_in_value[i].empty 
                         else pd.concat([lst_in_value[i],lst_inter_value[i]]) 
                         for i in range(len(lst_in_value))]
        
        if self.__mode == 'selection':
            ind = [e.nlargest(1).index[0] for e in lst_all_value if not e.empty]
            p = [self.__points.iloc[i] for i in ind]
            point_results = gpd.GeoDataFrame(p, geometry="geometry")
            
        if self.__mode == 'aggregation':
            aggreg = [(geom, len(cell)) for geom, cell in zip(self.__grid.centroid, lst_all_value) if not cell.empty]
            df = pd.DataFrame(aggreg, columns=["geometry", "radius"])
            point_results = gpd.GeoDataFrame(df, geometry="geometry")
            
        self.__points_res = point_results

    def getPointResults(self):
        """Get points results."""
        return self.__points_res
    
    def getGrid(self):
        """Get the grid."""
        return self.__grid
    
    def draw(self):
        """Draw the grid with the points given and the grid with point results."""
        fig, ax1 = plt.subplots()
        self.__grid.boundary.plot(ax=ax1, color='goldenrod')
        self.__points.plot(ax=ax1, color='black')
        
        fig, ax2 = plt.subplots()
        self.__grid.boundary.plot(ax=ax2, color='goldenrod')
        
        if self.__mode == 'selection':
            self.__points_res.plot(ax=ax2, color='chocolate')
            
        if self.__mode == 'aggregation':
            self.__points_res.plot(ax=ax2, color='chocolate', markersize=self.__width*self.__points_res["radius"])