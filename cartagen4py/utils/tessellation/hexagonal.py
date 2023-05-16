# classes to compute on hexagonal tessellation on top of a map

from shapely.geometry import Point, Polygon
import math
from shapely.strtree import STRtree

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
