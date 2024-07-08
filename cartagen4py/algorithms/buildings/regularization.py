import shapely
import numpy as np
from shapely import Point
from shapely.geometry.polygon import orient

import matplotlib.pyplot as plt

from cartagen4py.utils.geometry.angle import angle_2_pts

from cartagen4py.utils.debug import plot_debug

def rectangle_transformation(polygon, factor):
    """
    Transform a polygon to a rectangle.
    """
    mbr = polygon.minimum_rotated_rectangle
    return shapely.affinity.scale(mbr, xfact=factor, origin=mbr.centroid)

def recursive_regression(polygon, sigma):
    """
    Regularize a polygon using recursive linear regression.
    """

    if shapely.is_ccw(polygon):
        polygon = polygon.reverse()

    # Calculate the minimum rotated rectangle
    mbr = shapely.minimum_rotated_rectangle(polygon)

    # Retrieve first two vertex
    mbr_coords = list(mbr.exterior.coords)
    p1, p2 = Point(mbr_coords[0]), Point(mbr_coords[1])
    # Calulcate the angle of the first segment
    angle = angle_2_pts(p1, p2)

    # Rotate the polygon to reach horizontality
    rotated_mbr = shapely.affinity.rotate(mbr, -angle, origin=p1, use_radians=True)
    rotated = shapely.affinity.rotate(polygon, -angle, origin=p1, use_radians=True)

    if not shapely.is_ccw(rotated_mbr):
        rotated_mbr = rotated_mbr.reverse()

    # Get the coordinates of the four corner points
    mbr_coords = list(rotated_mbr.exterior.coords)
    extent = mbr_coords[:-1]

    # Get rotated polygon vertexes
    coords = list(rotated.exterior.coords)[:-1]

    corners = []
    for e in extent:
        distance = rotated_mbr.exterior.length 
        closest = None
        pe = Point(e)
        for i, vertex in enumerate(coords):
            d = shapely.distance(pe, Point(vertex))
            if d < distance:
                distance = d
                closest = i
        corners.append(closest)

    sides = []
    side = []
    for i, vertex in enumerate(coords):
        side.append(vertex)
        if i in corners:
            sides.append(side)
            side = [vertex]
    sides.append(side)

    last = sides.pop()
    first = sides.pop(0)
    sides.insert(0, last + first)
    
    for i, side in enumerate(sides):
        j = i + 1
        if j == len(extent):
            j = 0

        cstart, cend = extent[i], extent[j]

        side_angle = angle_2_pts(Point(cstart), Point(cend))

        if side_angle > 0:
            side_angle = side_angle - np.pi
        
        # Handle vertical sides
        if side_angle < 0:
            print(side_angle)

        start, end = np.array(side[0]), np.array(side[-1])
        x, y = [], []
        for vertex in side:
            x.append(vertex[0])
            y.append(vertex[1])

            v = np.array(vertex)
            distance = abs(np.cross(end - start, v - start) / np.linalg.norm(end - start))
        
        # Calculate the linear regression coefficients
        X = np.column_stack((np.ones_like(x), x))
        b, a = np.linalg.inv(X.T.dot(X)).dot(X.T.dot(y))

        b2, a2 = cstart[1] - (-1/a * cstart[0]), -1/a
        b3, a3 = cend[1] - (-1/a * cend[0]), -1/a

        x1 = (b2 - b) / (a - a2)
        y1 = a*x1+b
        x2 = (b3 - b) / (a - a3)
        y2 = a*x2+b

        distance = ((cstart[0] - x1)**2 + (cstart[1] - y1)**2)**0.5

        # print(distance)

        plot_debug(rotated, Point(cstart), Point(cend), shapely.LineString([(x1, y1), (x2, y2)]), [ Point(v) for v in side ])

    return polygon

def feature_edge_reconstruction(polygon):
    """
    Regularize a polygon using feature edge reconstruction.
    """