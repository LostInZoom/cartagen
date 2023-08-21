# This file contains functions that provide angle operations

from shapely.geometry import Point
import numpy as np
from math import pi

def angle_3_pts(point1, point2, point3):
    angle1 = angle_2_pts(point2, point1)
    angle2 = angle_2_pts(point2, point3)

    return angle2 - angle1

def angle_2_pts(point1, point2):
    x = point2.x - point1.x
    y = point2.y - point1.y
    return np.arctan2(y, x)

# converts an angle between 0 and 2Pi or in [-Pi, Pi] to the [0, Pi] interval
def angle_to_zero_pi(angle):
    if angle > pi:
        return angle - pi
    if angle < 0:
        return -angle
    return angle

# Computes the angle between two touching polylines.
def angle_two_connected_lines(line1, line2):
    # check that lines are touching
    if(line1.touches(line2) == False):
        return None
    # first get the intersection point
    intersection_pt = line1.intersection(line2)

    # get the following vertex on line1
    vertex1 = Point(line1.coords[len(line1.coords)-2])
    first1 = Point(line1.coords[0])
    if first1.equals(intersection_pt):
        vertex1 = Point(line1.coords[1])

    # get the following vertex on line2
    vertex2 = Point(line2.coords[len(line2.coords)-2])
    first2 = Point(line2.coords[0])
    if first2.equals(intersection_pt):
        vertex2 = Point(line2.coords[1])

    return angle_3_pts(vertex1, intersection_pt, vertex2)