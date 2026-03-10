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
    tessellate :
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
    def __init__(self, envelope, width):
        self.envelope = envelope
        self.width = width
        env_width = self.envelope.bounds[2] - self.envelope.bounds[0]
        env_height = self.envelope.bounds[3] - self.envelope.bounds[1]
        self.corner = Point(
            envelope.centroid.x - (env_width / 2), 
            envelope.centroid.y + (env_height / 2)
        )
        self.__compute_row_col_nb()
        self.__create_cells()
        # Create the strtree
        self._strtree = None
        self._geometries = None
    
    def __compute_row_col_nb(self):
        col_size = self.width * 3 / 4
        env_width = self.envelope.bounds[2] - self.envelope.bounds[0]
        env_height = self.envelope.bounds[3] - self.envelope.bounds[1]
        self.nb_columns = round(env_width / col_size) + 2
        self.nb_rows = round(env_height / col_size) * 2 + 3
    
    def __create_cells(self):
        self.cells = []
        for i in range(self.nb_rows):
            for j in range(self.nb_columns):
                # Skip pattern for hex grid
                if (i % 2 == 0 and j % 2 != 0) or (i % 2 != 0 and j % 2 == 0):
                    continue
                
                # Calculate hexagon center
                x_center = self.corner.x + (0.75 * self.width * (j - 1))
                y_center = (self.corner.y + (math.sqrt(3) * self.width / 4) - 
                           (math.sqrt(3) * self.width * (i - 1) / 4))
                
                self.cells.append(
                    HexagonalCell(self, i, j, Point(x_center, y_center), self.width)
                )
    
    def _build_strtree(self):
        """Build the STRtree, only call it once to avoid long calculation time"""
        if self._strtree is None:
            self._geometries = [cell.get_geometry() for cell in self.cells]
            self._strtree = STRtree(self._geometries)
    
    def get_containing_cells(self, point):
        """Find cells containing a given point"""
        self._build_strtree()
        
        point_geom = Point(point) if not isinstance(point, Point) else point
        indices = self._strtree.query(point_geom)
        
        return [self.cells[idx] for idx in indices]


class HexagonalCell:
    def __init__(self, tessellation, row, column, center, width):
        self.tessellation = tessellation
        self.row = row
        self.column = column
        self.center = center
        self.width = width
        self.radius = math.sqrt(3) / 2 * width
        self.segment_size = width / 2
        self._geometry = None  # Cache to avoid memory consumption
    
    def get_geometry(self):
        """Get hexagon geometry (use the cache when needed)"""
        if self._geometry is None:
            w = self.width
            cx, cy = self.center.x, self.center.y
            sqrt3_4 = math.sqrt(3) * w / 4
            
            coords = [
                (cx - w/2, cy),
                (cx - w/4, cy + sqrt3_4),
                (cx + w/4, cy + sqrt3_4),
                (cx + w/2, cy),
                (cx + w/4, cy - sqrt3_4),
                (cx - w/4, cy - sqrt3_4),
                (cx - w/2, cy)
            ]
            from shapely.geometry import Polygon
            self._geometry = Polygon(coords)
        
        return self._geometry
    
    def __eq__(self, other): 
        if not isinstance(other, HexagonalCell):
            return NotImplemented
        return self.row == other.row and self.column == other.column
    
    def __hash__(self):
        """Allows the use of sets instead of lists"""
        return hash((self.row, self.column))