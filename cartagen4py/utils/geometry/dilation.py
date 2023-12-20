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
        # Get the current point
        a = n

        if i >= (len(oline) - 1):
            # Get the previous point
            b = oline[i - 1]

            # Calculate the rounded extremity of the offset curve
            ba = np.sqrt((ya - yb)**2 + (xa - xb)**2)
            c = (xa - (((xa - xb) / ba) * (ba + abs(offset))), ya - (((ya - yb) / ba) * (ba + abs(offset))))

            # Interpolate points between the last projected point and the rounded extremity c
            rotation = 'cw' if offset < 0 else 'ccw'
            ip = circle_interpolation(a, pline[-1][-1], c, rotation=rotation)

            # Replace the last point with the interpolation
            pline[-1] = ip

            # Break the loop as it's the last point
            break

        # Get the next point
        b = oline[i + 1]
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
            ab = np.sqrt((yb - ya)**2 + (xb - xa)**2)
            c = (xb - (((xb - xa) / ab) * (ab + abs(offset))), yb - (((yb - ya) / ab) * (ab + abs(offset))))
            # Interpolate points between the projected a1 point and the rounded extremity c
            rotation = 'cw' if offset < 0 else 'ccw'
            ip = circle_interpolation(a, c, a1, rotation=rotation)
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
                pline.pop()
                # If they are not crossing, making a circle interpolation between those two points.
                rotation = 'cw' if offset < 0 else 'ccw'
                ip = circle_interpolation(a, b0, a1, rotation=rotation)
                pline.append(ip)

            pline.append([b1])

    segments = __create_segments(pline)

    # Create the full line object from groups of points
    fline = []
    # Loop through groups of points
    for nodes in pline:
        # Loop through each point
        for node in nodes:
            fline.append(node)
    
    # Loop through all vertex of the full line
    breakpoints = []
    breakstart = False
    parts = []
    part = []
    for n in range(0, len(fline) - 1):
        n1, n2 = fline[n], fline[n + 1]
        segment = shapely.LineString([n1, n2])

        cp = []
        for s in segments:
            if segment.crosses(s):
                cross = shapely.intersection(segment, s).coords[0]
                cp.append(cross)
                if cross not in breakpoints:
                    breakpoints.append(cross)

        part.append(n1)

        if len(cp) > 0:
            dpart = part[1:] if breakstart else part
            distance = round(__segment_distance(dpart, line), 8)

            if distance >= abs(offset):
                part.append(cp[0])
                parts.append(part)

            if len(cp) > 1:
                parts.append(cp)

            part = [cp[-1]]
            breakstart = True

        if (n + 1) == len(fline) - 1:
            part.append(n2)

    parts.append(part)

    groups = __merge_connected(parts)                  

    return parts

def circle_interpolation(a, b, c, rotation='cw', quad_segs=8):
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

    tangle = np.arccos(np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac)))
    
    n_points = int(tangle / (np.pi / 2) * quad_segs) + 1

    interpolated_points = [b]

    if n_points > 2:
        if rotation == 'cw':
            angle = - (tangle / (n_points - 1))
        elif rotation == 'ccw':
            angle = tangle / (n_points - 1)
        else:
            raise ValueError("Rotation value muste be clockwise ('cw') or counterclockwise ('ccw').")

        rotation_matrix = np.array(
            [[np.cos(angle), -np.sin(angle)],
            [np.sin(angle), np.cos(angle)]]
        )
        rotated_vector = ab / np.linalg.norm(ab)

        for i in range(1, n_points - 1):
            rotated_vector = np.dot(rotation_matrix, rotated_vector)
            interpolated_point = a + radius * rotated_vector
            interpolated_points.append(tuple(interpolated_point))
    
    interpolated_points.append(c)

    return interpolated_points

def __create_segments(points):
    """
    Create individual segments from the full nested list of points.
    """
    segments = []
    pp = None
    # Loop through groups of points
    for group in points:
        # Loop through each point
        for node in group:
            if pp is not None:
                segments.append(shapely.LineString([pp, node]))
            pp = node
    return segments

def __segment_distance(coords, line):
    """
    Return the minimum distance between the line and the coordinates provided.
    """
    distance = None

    for c in coords:
        point = shapely.Point(c)
        d = shapely.distance(point, line)

        if distance is None:
            distance = d
            continue

        if d < distance:
            distance = d
    
    return distance

def __merge_connected(groups):
    """
    Merge lines that have one or several points in common.
    """
    def __get_connected(group, groups):
        connected = []
        for point in group:
            for ig, g in enumerate(groups):
                if point in g and ig not in connected:
                    connected.append(ig)
        return connected

    for g in groups:
        c = __get_connected(g, groups)
        print(c)
