import numpy as np
import math
import shapely

import geopandas as gpd

from cartagen4py.utils.geometry.line import *

def offset_curve(line, offset, crs):
    """
    Offset a line using dilation on its left (positive offset) or right side (negative offset).
    """

    if offset == 0:
        return line

    oline = list(line.coords)
    pline = []

    previous = None
    
    for i, n in enumerate(oline):
        if i >= (len(oline) - 1):
            break

        # Get points
        a, b = n, oline[i + 1]
        xa, ya, xb, yb = a[0], a[1], b[0], b[1]

        # Get points as numpy arrays
        va, vb = np.array(a), np.array(b)

        # Calculate direction vector
        direction = va - vb
        # Normalize it
        direction /= np.linalg.norm(direction)

        # Calculate the orthogonal vector to get the projected vector
        proj = np.array([-direction[1], direction[0]])

        # Calculate the normalized buffer vector
        nbv = offset * proj / np.linalg.norm(proj)

        # Create the two new points
        a1 = tuple(a + nbv)
        b1 = tuple(b + nbv)

        # If it's the first node
        if i == 0:
            # Calculate the rounded extremity of the offset curve
            ab = np.sqrt((yb - ya)**2 + (xb -xa)**2)
            c = (xb - (((xb - xa) / ab) * (ab + abs(offset))), yb - (((yb - ya) / ab) * (ab + abs(offset))))
            # Interpolate points between the projected a1 point and the rounded extremity c
            ip = circle_interpolation(a, c, a1)
            # Add the interpolated points as the first entry...
            pline.append(ip)
            # ...and the end of the segment projected as the second
            pline.append([b1])
        # If it's a point between the first and the last
        else:
            # Retrieve the last projected point a and b
            a0, b0 = pline[-2][-1], pline[-1][-1]

            seg0 = shapely.LineString([a0, b0])
            seg1 = shapely.LineString([a1, b1])

            # Check is the current and previous segments are crossing
            if seg0.crosses(seg1):
                # If they are, it means the angle is acute and the last point needs to be recalculated
                # Calculate the intersection of both segments
                crosspoint = list(shapely.intersection(seg0, seg1).coords)[0]
                # Replace last point with this one
                pline[-1][-1] = (crosspoint[0], crosspoint[1])               
            else:
                pline[-1].pop()
                # If they are not crossing, making a circle interpolation between those two points.
                ip = circle_interpolation(a, b0, a1)
                pline.append(ip)

            pline.append([b1])


    final = []
    for nodes in pline:
        for node in nodes:
            final.append(node)

    # buffer = line.buffer(offset, quad_segs=8, cap_style="round", join_style='mitre')
    # boundary = buffer.boundary

    # print(boundary)

    result = extend_line_by_length(line, offset)

    points = []
    for i, nodes in enumerate(pline):
        for node in nodes:
            points.append({
                'id': i,
                'geometry': shapely.Point(node)
                })

    pgdf = gpd.GeoDataFrame(points, crs=crs)
    pgdf.to_file("cartagen4py/data/points.geojson", driver='GeoJSON')

    return shapely.LineString(final)

def circle_interpolation(a, b, c, quad_segs=8):
    """
    Given b and c two points at equal distance from a third point a, interpolates n points between
    b and c along the circle of center a and of radius ab. The number of provided point depends on the
    number of segments allowed per circle quadrant (default to 8).
    """
    # Create vectors
    ab = np.array(b) - np.array(a)
    ac = np.array(c) - np.array(a)

    # Calculate circle radius
    radius = np.linalg.norm(ab)

    if not np.isclose(np.linalg.norm(ab), np.linalg.norm(ac)):
        raise ValueError("Points b and c are not equidistant from a.")

    # Calculate the angle 
    angle = np.arccos(np.dot(ab / radius, ac / radius))
    n_points = int(angle / (np.pi / 2) * quad_segs) + 1

    print(n_points)

    interpolated_points = [b]

    if n_points > 2:
        angle_between_points = angle / (n_points - 1)
        for i in range(1, n_points - 1):
            theta = angle_between_points * i
            rotated_vector = np.cos(theta) * (ab / radius) + np.sin(theta) * (ac / radius)
            interpolated_point = a + radius * rotated_vector
            interpolated_points.append(tuple(interpolated_point))
    
    interpolated_points.append(c)

    return interpolated_points