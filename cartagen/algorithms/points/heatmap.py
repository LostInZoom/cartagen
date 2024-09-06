import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

def heatmap(points, cell_size, radius, column=None, method='quartic', clip=None):
    """
    Create a heatmap using the kernel density estimation technique (KDE).

    This function performs a spatial smoothing with the kernel density
    estimation technique (KDE), also known as heatmap. More information
    about heatmaps can be found in Wilkinson & Friendly. :footcite:p:`wilkinson:2009`
    
    For more information about KDE,
    here is a link to the related `Wikipedia article.
    <https://en.wikipedia.org/wiki/Kernel_density_estimation>`_
    This code is partially based on `this script.
    <https://niitdigital.medium.com/how-to-create-a-heatmap-from-scratch-in-python-234f602e856e>`_

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points used to calculates density.
    cell_size : int
        The size of the cell of the grid
        containing density values.
        Smaller size means smoother result,
        but also higher computation time.
    radius : int
        The radius used for the density calculation
        in each grid cells. For each centroid
        of grid cell, all the points within the
        radius are taken in account for density calculation.
        Higher radius means more generalized results.
    column : str, optional 
        Name of the column of the
        point to use to weight the density value.
    method : str, optional
        Name of the smoothing method that calculates
        the density value of each point within the radius.
        Each method impacts the way distance is
        important in the density calculation.
        Default to 'quartic'.
    clip : GeoDataFrame of Polygon, optional
        Polygons to clip the resulting heatmap grid.
        Be aware that it can return MultiPolygon.

    Returns
    -------
    grid : GeoDataFrame of Polygon
        The grid containing the values of density

    See Also
    --------
    reduce_quadtree :
        Reduce a set of points using a quadtree.
    reduce_labelgrid :
        Reduce a set of points using the Label Grid method.

    References
    ----------
    .. footbibliography::
    """
    #create square grid 
    xmin = points.total_bounds[0]
    xmax = points.total_bounds[2]
    ymin = points.total_bounds[1]
    ymax = points.total_bounds[3]

    cols = list(np.arange(xmin, xmax + cell_size, cell_size))
    rows = list(np.arange(ymin, ymax + cell_size, cell_size))

    polygons = gpd.GeoSeries([Polygon([(x, y), (x + cell_size, y), (x + cell_size, y + cell_size), 
                        (x, y + cell_size)]) for x in cols[:-1] for y in rows[:-1]])
    
    #retrieving grid cells centroid
    centroids = polygons.centroid
    
    #retrieving value attribute and store in a list
    if column is not None:
        lst_values = list(points[column])

    #storing x and y value of points and centroid into lists
    centroids_x = list(centroids.x)
    centroids_y = list(centroids.y)
    points_x = list(points.geometry.x)
    points_y = list(points.geometry.y)

    #definition de la fonction de lissage par méthode des noyaux (quadratique)
    if method == 'quartic':
        def kernel(d,radius):
            dn=d/radius
            P=(15/16)*(1-dn**2)**2
            return P

    #density calculation
    intensity_list = [] # list to store density value of each cell
    for i in range(len(centroids)):
        kde_value_list = []
        for k in range(len(points_x)):
            d = np.sqrt((centroids_x[i]-points_x[k])**2+(centroids_y[i]-points_y[k])**2)
            if d <= radius:
                p = kernel(d,radius)
                if column is None:
                    kde_value_list.append(p)
                else: 
                    kde_value_list.append(p*lst_values[k])

        total_density = sum(kde_value_list)
        intensity_list.append(total_density)

    final_grid = gpd.GeoDataFrame({'geometry': polygons,'density':intensity_list}, crs=points.crs)

    if clip is not None:
        return final_grid.clip(clip)
    else:
        return final_grid