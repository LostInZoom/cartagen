# this file define the Segment class, which representation a mathematical segment

from shapely import Point
from util import angle_operations
import math
import numpy as np

# get the list of segments of a linestring or a polygon
def get_segment_list(geometry):
    segment_list = []
    for i in range(1,len(geometry.coords)-1):
        segment_list.append(Segment(geometry.coords[i-1],geometry.coords[i]))
    return segment_list

class Segment:
    point1
    point2
    __coefA = 0
    __coefB = 0
    __coefC = 0

    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
        self.__coefA = point2[1] - point1[1]
        self.__coefB = point1[0] - point2[0]
        self.__coefC = point1[1] * (point2[0]-point1[0]) + point1[0] * (point1[1] - point2[1])
    
    def orientation(self):
        xAxisPt = Point(point1[0] + 10.0, point1[1])
        return angle_operations.angle_3_pts(xAxisPt,self.point1,self.point2)
    
    def length(self):
        return math.sqrt((self.point2[0]-self.point1[0])*(self.point2[0]-self.point1[0]) + (self.point2[1]-self.point1[1])*(self.point2[1]-self.point1[1]))

    # compute the intersection point of two segments extended as straight lines
    def straight_line_intersection(self, segment2):
        matrix = np.zeros(2,2)
        matrix.set(0,0,self.__coefA)
        matrix.set(1,0,segment2.__coefA)
        matrix.set(0,1,self.__coefB)
        matrix.set(1,1,segment2.__coefB)

        # check if straight lines are parallel, i.e. matrix determinant is zero
        # TODO

        # inverse the matrix
        # TODO
        inverse = None

        x_intersection = - inverse.get(0,0)*self.__coefC - inverse.get(0,1)*segment2.__coefC
        y_intersection = - inverse.get(1,0)*self.__coefC - inverse.get(1,1)*segment2.__coefC

        return Point(x_intersection, y_intersection)
