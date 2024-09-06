from shapely.geometry import Point
import numpy as np
from math import pi

def angle_2_pts(point1, point2):
    """Returns the angle formed by two points"""
    x = point2.x - point1.x
    y = point2.y - point1.y
    return np.arctan2(y, x)

def angle_3_pts(point1, point2, point3):
    """Returns the angle formed by three points"""
    angle1 = angle_2_pts(point2, point1)
    angle2 = angle_2_pts(point2, point3)
    return angle2 - angle1

# converts an angle between 0 and 2Pi or in [-Pi, Pi] to the [0, Pi] interval
def angle_to_zero_pi(angle):
    if angle > pi:
        return angle - pi
    if angle < 0:
        return - angle
    return angle

def angle_between_2lines(line1, line2):
    """
    Return the angle between two lines that crosses at their end or start point. Takes shapely LineStrings as input.
    """
    # retrieve coordinates of lines nodes, number of node in lines and start and end node
    coords1 = line1.coords
    l1 = len(coords1)
    start1, end1 = coords1[0], coords1[-1]
    coords2 = line2.coords
    l2 = len(coords2)
    start2, end2 = coords2[0], coords2[-1]

    sgeom1 = sgeom2 = None
    current = None

    if start1 == start2:
        current = start1
        sgeom1 = True
        sgeom2 = True
    if start1 == end2:
        current = start1
        sgeom1 = True
        sgeom2 = False
    if end1 == end2:
        current = end1
        sgeom1 = False
        sgeom2 = False
    if end1 == start2:
        current = end1
        sgeom1 = False
        sgeom2 = True

    # In this case, lines are not crossing
    if current is None:
        return None

    previous = following = None

    if l1 > 2:
        if sgeom1:
            previous = coords1[2]
        else:
            previous = coords1[-3]
    else:
        if sgeom1:
            previous = coords1[1]
        else:
            previous = coords1[-2]

    if l2 > 2:
        if sgeom2:
            following = coords2[2]
        else:
            following = coords2[-3]
    else:
        if sgeom2:
            following = coords2[1]
        else:
            following = coords2[-2]

    angle = angle_3_pts(Point(previous), Point(current), Point(following))

    return angle

def get_curvature(p1, p2, p3):
    """
    Return the curvature, i.e, the inverse of the radius of the circumcircle
    passing by the three given points.
    """
    return 1 / abs(p2.distance(p3) / (2 * np.sin(angle_3_pts(p3, p2, p1))))