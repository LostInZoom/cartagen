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

def angle_from_3points_coordinates(p1, p2, p3):
    """
    Returns the angle formed by three points. The returned value is in degrees always positive,
    it doesn't provide any information about orientation.
    p1, p2, p3 must be python iterables of two elements representing the coordinates of the points.
    """
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]
    x3, y3 = p3[0], p3[1]

    deg1 = (360 + np.rad2deg(np.arctan2(x1 - x2, y1 - y2))) % 360
    deg2 = (360 + np.rad2deg(np.arctan2(x3 - x2, y3 - y2))) % 360

    return deg2 - deg1 if deg1 <= deg2 else 360 - (deg1 - deg2)