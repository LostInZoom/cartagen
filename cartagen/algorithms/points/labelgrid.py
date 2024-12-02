from cartagen.utils.partitioning.tessellation import tessellate
import warnings
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

def labelgrid_selection(points, width, height, column, shape='square', grid=False):
    """
    Reduce a set of points by selection using the Label Grid method.

    For each cell of the given shape, the point with the largest value is selected.

    This algorithm was proposed by Gröbe & Burghardt. :footcite:p:`grobe:2021`
    It assigns each point to a grid cell of a given height, width and shape
    before reducing the points for each cell.

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    width : float
        Width of the grid cells.
    height : float
        Height of the grid cells.
    column : str
        Name of the column to use.
    shape : str, optional
        Shape of the grid cells, can be 'square', 'diamond', 'hexagonal'.
    grid : bool, optional
        If set to True, returns a tuple with the points and the grid cell.

    Returns
    -------
    points : GeoDataFrame of Point
        The provided set of points with two new attributes:

        - 'selected_labelgrid': Set to True if the point has been selected.
        - 'cell_count': The number of points in the cell represented by this point.
    
    grid : GeoDataFrame of Polygon, optional

    See Also
    --------
    labelgrid_simplification :
        Reduce a set of points by simplification using the Label Grid method.
    labelgrid_aggregation :
        Reduce a set of points by aggregation using the Label Grid method.

    References
    ----------
    .. footbibliography::
    """
    return __reduce_labelgrid(points, width, height, shape, mode='selection', column=column, grid=grid)

def labelgrid_simplification(points, width, height, shape='square', grid=False):
    """
    Reduce a set of points by simplification using the Label Grid method.

    For each cell of the given shape, the point closest to the grid's centroid is selected.

    This algorithm was proposed by Gröbe & Burghardt. :footcite:p:`grobe:2021`
    It assigns each point to a grid cell of a given height, width and shape
    before reducing the points for each cell.

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    width : float
        Width of the grid cells.
    height : float
        Height of the grid cells.
    shape : str, optional
        Shape of the grid cells, can be 'square', 'diamond', 'hexagonal'.
    grid : bool, optional
        If set to True, returns a tuple with the points and the grid cell.

    Returns
    -------
    points : GeoDataFrame of Point
        The provided set of points with two new attributes:

        - 'selected_labelgrid': Set to True if the point has been selected.
    
    grid : GeoDataFrame of Polygon, optional

    See Also
    --------
    labelgrid_selection :
        Reduce a set of points by selection using the Label Grid method.
    labelgrid_aggregation :
        Reduce a set of points by aggregation using the Label Grid method.

    References
    ----------
    .. footbibliography::
    """
    return __reduce_labelgrid(points, width, height, shape, mode='simplification', column=None, grid=grid)

def labelgrid_aggregation(points, width, height, column=None, shape='square', grid=False):
    """
    Reduce a set of points by aggregation using the Label Grid method.

    For each cell of the given shape, the point closest to the grid's centroid is selected.

    This algorithm was proposed by Gröbe & Burghardt. :footcite:p:`grobe:2021`
    It assigns each point to a grid cell of a given height, width and shape
    before reducing the points for each cell.

    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    width : float
        Width of the grid cells.
    height : float
        Height of the grid cells.
    column : str, optional
        Name of a numeric column to calculate the sum and mean.
    shape : str, optional
        Shape of the grid cells, can be 'square', 'diamond', 'hexagonal'.
    grid : bool, optional
        If set to True, returns a tuple with the points and the grid cell.

    Returns
    -------
    points : GeoDataFrame of Point
        The provided set of points with two new attributes:

        - 'selected_labelgrid': Set to True if the point has been selected.
        - 'column_total': The sum of the column values if provided.
        - 'column_mean': The mean of the column values if provided.
    
    grid : GeoDataFrame of Polygon, optional

    See Also
    --------
    labelgrid_selection :
        Reduce a set of points by selection using the Label Grid method.
    labelgrid_simplification :
        Reduce a set of points by simplification using the Label Grid method.

    References
    ----------
    .. footbibliography::
    """
    return __reduce_labelgrid(points, width, height, shape, mode='aggregation', column=column, grid=grid)

def __reduce_labelgrid(
        points, width, height,
        shape='square', mode='simplification',
        column=None, grid=False
    ):
    """
    Reduce a set of points using the Label Grid method.

    This algorithm was proposed by Gröbe & Burghardt. :footcite:p:`grobe:2021`
    It assigns each point to a grid cell of a given height, width and shape
    before reducing the points per cells either by selection or aggregation.
    
    Parameters
    ----------
    points : GeoDataFrame of Point
        The points to reduce.
    width : float
        Width of the grid cells.
    height : float
        Height of the grid cells.
    shape : str, optional
        Shape of the grid cells, can be 'square', 'diamond', 'hexagonal'.
    mode : str, optional
        There are three available modes:

        - *'selection'*: for one cell, the algorithm retains
          the point with the largest value in the chosen column. This option requires
          the column parameter to be provided.
        - *'simplification'*: the point retained in the cell
          is the closest to the center of the cell.
        - *'aggregation'*: the points are all aggregated to
          the centroid of the cell. The count of point is added as a new attribute.
          If a column name is provided, also adds the sum of the attribute.

    column : str, optional
        Name of the column to use.
    grid : bool, optional
        If set to True, returns a tuple with the points and the grid.
        In aggregation mode, the returned grid cells contain
        the count of point aggregated along with sum and mean.

    Returns
    -------
    GeoDataFrame of Point or tuple (GeoDataFrame of Point, GeoDataFrame of Polygon)
    
    References
    ----------
    .. footbibliography::
    """
    if shape not in ['square', 'diamond', 'hexagonal']:
        raise Exception('{0} shape is not recognized.'.format(shape))

    if mode not in ['selection', 'simplification', 'aggregation']:
        raise Exception('{0} mode is not recognized.'.format(mode))

    if mode == 'selection' and column is None:
        raise Exception('Selection mode requires an attribute.')

    # if column not in list(points):  #-> always raise exception when column not provided
    #     raise Exception('Column {0} not in the provided GeoDataFrame.'.format(column))

    if column is not None and mode == 'simplification':
        warnings.warn("Warning: There is no need to indicate a column name in simplification mode.")

    lg = LabelGrid(points, width, height, shape, mode, column)
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
    
    def __init__(self, points, width, height, shape='square', mode='simplification', column=None):
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
        shape : str
            Form of a cell. The default is 'square'.
        """
        self.__column = column
        self.__width = width
        self.__height = height
        self.__shape = shape
        self.__mode = mode
        self.__points = points
        self.__points_res = None
        self.__grid = None
        self.__crs = points.crs

    def __createGrid(self):
        """
        Create the grid of the label grid.
        """
        polygons = tessellate(self.__points.total_bounds, self.__width, self.__height, self.__shape)
        return gpd.GeoDataFrame({'geometry': polygons}, crs=self.__crs)

    def set_point_label_grid(self):
        """Set the attributes points_res and grid of the LabelGrid object."""
        self.__grid = self.__createGrid()
        
        lst_in_value = [self.__points.loc[cell.contains(self.__points['geometry'])] 
                        for cell in self.__grid['geometry']]

        lst_inter_value = [self.__points.loc[cell.touches(self.__points['geometry'])] 
                           for cell in self.__grid['geometry']]
        
        lst_all_value = [lst_inter_value[i] 
                         if lst_in_value[i].empty 
                         else pd.concat([lst_in_value[i],lst_inter_value[i]]) 
                         for i in range(len(lst_in_value))]
               
        if self.__mode == 'simplification':
            result = []
            for c, centroid in enumerate(self.__grid.centroid):
                points = lst_all_value[c]
                if not points.empty:
                    mindist = float("inf")
                    nearest = None
                    for index, point in points.iterrows():
                        dist = point['geometry'].distance(centroid)
                        if dist < mindist:
                            mindist = dist
                            nearest = index

                    for index, p in points.iterrows():
                        select = False
                        if index == nearest:
                            select = True
                        p['cell'] = c
                        p['selected_labelgrid'] = select
                        result.append(p)

            self.__points_res = gpd.GeoDataFrame(result, crs=self.__crs)

        elif self.__mode == 'selection':
            result = []
            for c, cell in self.__grid.iterrows():
                points = lst_all_value[c]
                if not points.empty:
                    # retain the largest value in each cell for the chosen attribute of the point
                    selected = None
                    largest = 0
                    for index, point in points.iterrows():
                        if point[self.__column] > largest:
                            selected = index
                            largest = point[self.__column]
                    
                    for index, p in points.iterrows():
                        select = False
                        if index == selected:
                            select = True
                        p['cell'] = c
                        p['cell_count'] = len(points)
                        p['selected_labelgrid'] = select
                        result.append(p)

            self.__points_res = gpd.GeoDataFrame(result, crs=self.__crs)
            
        elif self.__mode == 'aggregation':
            aggregation = []
            for c, centroid in enumerate(self.__grid.centroid):
                p = lst_all_value[c]
                entry = { "cell_count": len(p) }
                if self.__column is not None:
                    entry["{0}_total".format(self.__column)] = p[self.__column].sum()
                    entry["{0}_mean".format(self.__column)] = p[self.__column].mean()

                entry["geometry"] = centroid
                entry["cell"] = c
                if len(p) > 0: #only keep centroid of cells that contains points
                    aggregation.append(entry)
            self.__points_res = gpd.GeoDataFrame(aggregation, crs=self.__crs)

    def getPointResults(self):
        """Get points results."""
        return self.__points_res
    
    def getGrid(self):
        """Get the grid."""
        if self.__mode == 'aggregation':
            grid = self.__grid.to_dict('records')
            for p, point in self.__points_res.iterrows():
                grid[point['cell']]["count"] = point["cell_count"]
                if self.__column is not None:
                    grid[point['cell']]["{0}_total".format(self.__column)] = point["{0}_total".format(self.__column)]
                    grid[point['cell']]["{0}_mean".format(self.__column)] = point["{0}_mean".format(self.__column)]
            self.__grid = gpd.GeoDataFrame(grid, crs=self.__crs)
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
            self.__points_res.plot(ax=ax2, color='chocolate', markersize=self.__width*self.__points_res["count"])