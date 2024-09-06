import math
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.strtree import STRtree

def partition_grid(objects, width, height=None, shape='square'):
    """
    Partition objects using a grid of a given shape.

    This algorithm divides the extent of the provided objects
    using a regular grid of the specified shape and size and assign
    the objects to the cell they intersects.

    Parameters
    ----------
    objects : GeoDataFrame
        The objects to partition.
    width : float
        Width of the grid.
    height : float, optional
        Height of the grid.
        If set to None, the height equals the width.
    shape : str, optional
        Shape of the grid cells, can be 'square', 'diamond', 'hexagonal'.

    Returns
    -------
    partition : tuple
        A tuple containing two elements :

        #. A list of lists of index ordered by the grid cells
        #. A list of the geometry of the grid cells

    See Also
    --------
    tesselate :
        Create a tessellation of cells of a given shape.
    partition_networks :
        Partition objects using one or multiple networks.

    Examples
    --------
    >>> points = gpd.GeoDataFrame(geometry=[ Point(2, 1), Point(4, 2) ])
    >>> partition_grid(points, 1)
    ([[0], [1]], [<POLYGON ((2 1, 3 1, 3 2, 2 2, 2 1))>, <POLYGON ((3 1, 4 1, 4 2, 3 2, 3 1))>])
    """
    partition = ([], [])
    grid = tessellate(objects.total_bounds, width, height, shape)

    centroids = []
    for i, o in objects.iterrows():
        if o['geometry'].geom_type == 'Polygon':
            centroids.append(o['geometry'].point_on_surface())
        else:
            centroids.append(o['geometry'].centroid)

    done = []
    for cell in grid:
        contained = []
        for i, c in enumerate(centroids):
            # Using intersects to ensure every vertex is assigned
            if cell.intersects(c):
                # Check if it is not already added
                if i not in done:
                    contained.append(i)
                    done.append(i)
        partition[0].append(contained)
        partition[1].append(cell)
    
    return partition

def tessellate(extent, width, height=None, shape='square'):
    """
    Create a tessellation of cells of a given shape.

    This function creates a tessellation on the provided extent
    with cells of the provided shape.

    Parameters
    ----------
    extent : list of float
        The extent of the area to tessellate.
        Must be of the form [ xmin, ymin, xmax, ymax ]
    width : float
        Width of the grid.
    height : float, optional
        Height of the grid.
        If set to None, the height equals the width.
    shape : str, optional
        Shape of the grid cells, can be 'square', 'diamond', 'hexagonal'.

    Returns
    -------
    list of Polygon

    See Also
    --------
    partition_grid :
        Partition objects using a grid of a given shape.

    Examples
    --------
    >>> tessellate([0, 0, 2, 1], 1)
    [<POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))>, <POLYGON ((1 0, 2 0, 2 1, 1 1, 1 0))>]
    """
    def __create_hexagons(coords_min, coords_max):
        xstep = width*3
        ystep = width*2
        cols = list(np.arange(coords_min[0], coords_max[0]+xstep, xstep))
        rows = list(np.arange(coords_min[1], coords_max[1]+ystep, ystep))
        
        hexagons = [
            Polygon([ 
                (x, y),
                (x + width, y),
                (x + width*3/2, y + width),
                (x + width, y + 2*width),
                (x, y + 2*width),
                (x - width/2, y + width)
            ]) 
            for x in cols[:-1] for y in rows[:-1]
        ]

        return hexagons, len(rows) - 1, len(cols) - 1

    xmin, ymin, xmax, ymax = extent

    if height is None:
        height = width

    if shape == 'square':
        cols = list(np.arange(xmin, xmax + width, width))
        rows = list(np.arange(ymin, ymax + height, height))
    
        return [
            Polygon([
                (x, y),
                (x + width, y),
                (x + width, y + height),
                (x, y + height)
            ]) 
            for x in cols[:-1] for y in rows[:-1]
        ]

    elif shape == 'hexagonal':
        odd_coords_min = (xmin - width, ymin - width)
        odd_coords_max = (xmax + width, ymax)
        
        odd_poly, odd_rows, odd_cols = __create_hexagons(odd_coords_min, odd_coords_max)
        
        even_coords_min = (odd_coords_min[0] + 3/2*width, odd_coords_min[1] - width)
        even_coords_max = (odd_coords_max[0], odd_coords_max[1] + width)
        
        even_poly, even_rows, even_cols = __create_hexagons(even_coords_min, even_coords_max)

        odd_chunks = [ odd_poly[i:i + odd_rows] for i in range(0, len(odd_poly), odd_rows) ]
        even_chunks = [ even_poly[i:i + even_rows] for i in range(0, len(even_poly), even_rows) ]

        result = [None]*(len(odd_chunks) + len(even_chunks))
        result[::2] = odd_chunks
        result[1::2] = even_chunks

        ordered = []
        for r in result:
            for e in r:
                ordered.append(e)
        
        return ordered
    
    elif shape == 'diamond':
        cols = list(np.arange(xmin, xmax + width*2, width))
        rows = list(np.arange(ymin - width, ymax + width, width))
        
        return [
            Polygon([
                (x, y), 
                (x + width / 2, y + width / 2), 
                (x, y + width), 
                (x - width / 2, y + width / 2)
            ]) 
            for x in cols[:-1] for y in rows[:-1]
        ]

    else:
        raise Exception('{0} shape not recognized.'.format(shape))

class HexagonalTessellation:
    corner = None # the top left corner of the tessellation
    width = 0 # the width of the cells in the tessellation
    envelope = None # the envelope of the tessellation
    nb_columns = 0
    nb_rows = 0
    cells = []

    def __init__(self, envelope, width):
        self.envelope = envelope
        self.width = width
        env_width = self.envelope.bounds[2] - self.envelope.bounds[0]
        env_height = self.envelope.bounds[3] - self.envelope.bounds[1]
        self.corner = Point(envelope.centroid.x
        - (env_width / 2), envelope.centroid.y
        + (env_height / 2))
        self.__compute_row_col_nb()
        self.__create_cells()
    
    def __compute_row_col_nb(self):
        col_size = self.width * 3 / 4
        env_width = self.envelope.bounds[2] - self.envelope.bounds[0]
        env_height = self.envelope.bounds[3] - self.envelope.bounds[1]
        self.nb_columns = round(env_width / col_size) + 2
        self.nb_rows = round(env_height / col_size) * 2 + 3
        return
    
    def __create_cells(self):
        for i in range(0,self.nb_rows):
            for j in range(0,self.nb_columns):
                if (i % 2 == 0) and (j % 2 != 0):
                    continue
                if (i % 2 != 0) and (j % 2 == 0):
                    continue
                # now compute the center for hexagon (i,j)
                xCenter = self.corner.x + (0.75 * self.width * (j - 1))
                yCenter = self.corner.y + (math.sqrt(3) * self.width / 4) - (math.sqrt(3) * self.width * (i - 1) / 4)
                self.cells.append(HexagonalCell(self,i,j,Point(xCenter,yCenter),self.width))

        return
    
    def get_containing_cells(self, point):
        cells = []
        hexagons = []
        for cell in self.cells:
            hexagons.append(cell.get_geometry())
            
        tree = self.__get_STRtree(hexagons)
        hids = tree.query(Point(point))

        for hid in hids:
            cells.append(self.cells[hid])
            
        return cells
    
    def __get_STRtree(self, hexagons):
        tree = STRtree(hexagons)
        return tree

class HexagonalCell:
    row = 0
    column = 0
    tessellation = None
    center = None
    radius = 0
    width = 0
    segment_size = 0

    def __init__(self, tessellation, row, column, center, width):
        self.tessellation = tessellation
        self.row = row
        self.column = column
        self.center = center
        self.width = width
        self.radius = math.sqrt(3) / 2 * width
        self.segment_size = width / 2
    
    def get_geometry(self):
        coords = []
        coords.append((self.center.x - (self.width / 2), self.center.y))
        coords.append((self.center.x - (self.width / 4), (math.sqrt(3)* self.width / 4) + self.center.y))
        coords.append((self.center.x + (self.width / 4), (math.sqrt(3)* self.width / 4) + self.center.y))
        coords.append((self.center.x + (self.width / 2), self.center.y))
        coords.append((self.center.x + (self.width / 4), self.center.y - (math.sqrt(3) * self.width / 4)))
        coords.append((self.center.x - (self.width / 4), self.center.y - (math.sqrt(3) * self.width / 4)))
        coords.append((self.center.x - (self.width / 2), self.center.y))
        return Polygon(coords)
    
    def __eq__(self, other): 
        if not isinstance(other, HexagonalCell):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.row == other.row and self.column == other.column
