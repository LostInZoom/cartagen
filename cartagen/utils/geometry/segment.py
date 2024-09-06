# this file define the Segment class, which representation a mathematical segment

from shapely.geometry import Point, LineString
from cartagen.utils.geometry.angle import angle_3_pts
import math
import numpy as np

# get the list of segments of a linestring or a polygon
def get_segment_list(geometry):
    coords = []
    if(geometry.geom_type == 'LineString'):
        coords = list(geometry.coords)
    elif(geometry.geom_type == 'Polygon'):
        coords = list(geometry.exterior.coords)
    segment_list = []
    for i in range(1, len(coords)):
        segment_list.append(Segment(coords[i-1],coords[i]))
    return segment_list

# get the list of segments in a polygon, including the inner rings. Returns a list of pairs (segment, n) n being -1 for segments 
# in the outer ring, and i for the segments in inner rings (i is the index of the ring in the list of inner rings).
def get_segment_list_polygon(geometry):
    segment_list = []
    coords_outer = geometry.exterior.coords
    for i in range(1,len(coords_outer)):
        segment_list.append((Segment(coords_outer[i-1],coords_outer[i]), -1))
    if(len(geometry.interiors)>0):
        index = 0
        for interior in geometry.interiors:
            for i in range(1,len(interior.coords)):
                segment_list.append((Segment(interior.coords[i-1],interior.coords[i]),index))
            index += 1
    return segment_list

class Segment:
    point1 = []
    point2 = []
    __coefA = 0
    __coefB = 0
    __coefC = 0

    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
        self.__coefA = point2[1] - point1[1]
        self.__coefB = point1[0] - point2[0]
        self.__coefC = point1[1] * (point2[0]-point1[0]) + point1[0] * (point1[1] - point2[1])
    
    def get_coefA(self):
        return self.__coefA
    
    def get_coefB(self):
        return self.__coefB
    
    def get_coefC(self):
        return self.__coefC
    
    def get_geom(self):
        return LineString([self.point1, self.point2])

    def orientation(self):
        xAxisPt = Point(self.point1[0] + 10.0, self.point1[1])
        return angle_3_pts(xAxisPt,Point(self.point1),Point(self.point2))
    
    def length(self):
        return math.sqrt((self.point2[0]-self.point1[0])*(self.point2[0]-self.point1[0]) + (self.point2[1]-self.point1[1])*(self.point2[1]-self.point1[1]))

    # compute the intersection point of two segments extended as straight lines
    def straight_line_intersection(self, segment2):
        matrix = np.zeros((2,2))
        matrix[(0, 0)] = self.__coefA
        matrix[(1, 0)] = segment2.__coefA
        matrix[(0, 1)] = self.__coefB
        matrix[(1, 1)] = segment2.__coefB

        # check if straight lines are parallel, i.e. matrix determinant is zero
        det = np.linalg.det(matrix)
        if (det==0.0):
            return None

        # inverse the matrix
        inverse = np.linalg.inv(matrix)

        x_intersection = - inverse.item(0,0)*self.__coefC - inverse.item(0,1)*segment2.__coefC
        y_intersection = - inverse.item(1,0)*self.__coefC - inverse.item(1,1)*segment2.__coefC

        return Point(x_intersection, y_intersection)
    
