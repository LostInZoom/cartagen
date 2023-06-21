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