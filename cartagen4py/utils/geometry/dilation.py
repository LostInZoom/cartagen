import numpy as np
import math
import shapely

import geopandas as gpd

from cartagen4py.utils.geometry.line import *

def offset_curve(line, offset, cap_style='round', quad_segs=8):
    """
    Offset a line using dilation on its left (positive offset) or right side (negative offset).
    Parameters
    ----------
    line : shapely LineString
        The line to offset.
    offset : float
        The length of the offset to apply. Negative value for left-side dilation, positive for right-side.
    cap_style : str optional
        The type of caps at the start and end of the provided line. Possible values are 'round' or 'flat'.
    quad_segs : int optional
        The number of point allowed per circle quadrant when interpolating points using round method.
    """

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
                # Add the segment if the previous point exists
                if pp is not None:
                    segments.append(shapely.LineString([pp, node]))
                # Set the previous point as the current
                pp = node

        return segments

    def __segment_distance(coords, line):
        """
        Return the minimum distance between the line and the coordinates provided.
        """
        # Set the distance to None
        distance = None

        # Loop through coordinates
        for c in coords:
            # Create the point
            point = shapely.Point(c)
            # Calculate the distance from the line
            d = shapely.distance(point, line)

            # If the distance is none, it's the first iteration
            if distance is None:
                # Set the distance to the current and continue
                distance = d
                continue

            # Update the distance if this one is lower
            if d < distance:
                distance = d
        
        return distance

    def __merge_connected(groups):
        """
        From a list of groups, recreate groups by merging connected groups and
        removing duplicates.
        """
        result = []

        # Loop through each groups
        for igroup, group in enumerate(groups):
            # Storage for already merged groups
            merged = []
            # Loop through result groups
            for sub in result:
                # If common coordinates are found, add the group to the list to merge
                if any(coord in sub for coord in group):
                    merged.append(sub)
            
            # If no groups were found with common coordinates, add the group to the result
            if not merged:
                result.append(group)
            # If groups are found with common coordinates
            else:
                # Storage for the new group
                new = []
                # Loop through groups with common coordinates
                for sub in merged:
                    # Loop through coordinates inside the group
                    for coord in sub:
                        # Add them to the new group
                        new.append(coord)
                    # Loop through coordinates inside the current group
                    for coord in group:
                        # Add them to the new group
                        new.append(coord)

                # Rebuild the result by removing groups merged inside the new one
                result = [g for g in result if g not in merged]

                # Add the new group to the list
                result.append(new)

        for i, group in enumerate(result):
            if len(group) > 2:
                # Destructuring the list to get start, middle and end
                start, *middle, end = group

                # Get the list without duplicates
                unique = list(set(group))
                # If it's only three coordinates
                if len(unique) == 3:
                    # Replace the result group without the duplicates and add the start at the end
                    # This should handle double crossing segments in holes
                    result[i] = unique + [unique[0]]
                else:
                    # Create a set to store coordinates
                    seen = set()
                    # Recreate a list keeping the start and the end untouched as they can't be problematic
                    # Remove all duplicates in between start and end
                    unique_middle = [start] + [coord for coord in middle if coord not in seen and not seen.add(coord)] + [end]
                    # Replace the resulting group
                    result[i] = unique_middle

        return result

    # OFFSET CURVE ALGORITHM

    # Check if provided cap style is allowed
    acaps = ['round', 'flat']
    if cap_style not in acaps:
        raise Exception('Offset curve cap style unknown: {0}.'.format(cap_style))

    # If offset is 0, return the provided line
    if offset == 0:
        return line

    # Get the original line as a list of vertex coordinates
    oline = list(line.coords)
    # Storage for the projected line
    pline = []

    # Loop through vertex of the original line
    for i, n in enumerate(oline):
        # Get the current point
        a = n

        # If it's the last point of the line
        if i >= (len(oline) - 1):
            # If round cap style is chosen
            if cap_style == 'round':
                # Get the previous point
                b = oline[i - 1]

                # Calculate the rounded extremity of the offset curve
                ba = np.sqrt((ya - yb)**2 + (xa - xb)**2)
                c = (xa - (((xa - xb) / ba) * (ba + abs(offset))), ya - (((ya - yb) / ba) * (ba + abs(offset))))

                # Interpolate points between the last projected point and the rounded extremity c
                rotation = 'cw' if offset < 0 else 'ccw'
                ip = circle_interpolation(a, pline[-1][-1], c, rotation=rotation, quad_segs=quad_segs)

                # Replace the last point with the interpolation
                pline[-1] = ip

            # Do nothing if the cap style is flat

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
            
            # If round cap style is chosen
            if cap_style == 'round':
                # Calculate the rounded extremity of the offset curve
                ab = np.sqrt((yb - ya)**2 + (xb - xa)**2)
                c = (xb - (((xb - xa) / ab) * (ab + abs(offset))), yb - (((yb - ya) / ab) * (ab + abs(offset))))

                # Interpolate points between the projected a1 point and the rounded extremity c
                rotation = 'cw' if offset < 0 else 'ccw'
                ip = circle_interpolation(a, c, a1, rotation=rotation, quad_segs=quad_segs)

                # Add the interpolated points as the first entry
                pline.append(ip)

            # If a flat cap style is chosen
            elif cap_style == 'flat':
                # Add the projection of a
                pline.append([a1])
            
            # Add the projection of b
            pline.append([b1])

        # If it's a point between the first and the last (middle point)
        else:
            # Retrieve the last projected point a and b
            a0, b0 = pline[-2][-1], pline[-1][-1]

            # Create the previous and the current segment
            seg0 = shapely.LineString([a0, b0])
            seg1 = shapely.LineString([a1, b1])

            # Check if the current and previous segments are crossing
            if seg0.crosses(seg1):
                # If they are, it means the angle is acute and the last point needs to be recalculated
                # Calculate the intersection of both segments
                crosspoint = list(shapely.intersection(seg0, seg1).coords)[0]
                # Replace last point with this one
                pline[-1][-1] = (crosspoint[0], crosspoint[1])    

            # If they are not           
            else:
                # Remove the last projection of b as the circle interpolation will add it again
                pline.pop()
                # Making a circle interpolation between those two points.
                rotation = 'cw' if offset < 0 else 'ccw'
                ip = circle_interpolation(a, b0, a1, rotation=rotation, quad_segs=quad_segs)
                # Add the interpolated points to the line
                pline.append(ip)

            # Add the projection of b to the list
            pline.append([b1])

    # Get the full line as single segments
    segments = __create_segments(pline)

    # Create the full line object from groups of points
    fline = []
    # Loop through groups of points
    for nodes in pline:
        # Loop through each point
        for node in nodes:
            fline.append(node)
    
    # Flag to see if the current line starts with a point where segment crosses
    breakstart = False
    # Storage for the final groups of individual parts (below)
    parts = []
    # Storage for individual groups of nodes forming a continuous line
    part = []

    # Loop through all vertex of the full line
    for n in range(0, len(fline) - 1):
        # Get the current node and the following
        n1, n2 = fline[n], fline[n + 1]
        # Create the segment
        segment = shapely.LineString([n1, n2])

        # Storage for points where segments crosses
        cp = []
        # Loop through segments
        for s in segments:
            # If the current segment crosses an other
            if segment.crosses(s):
                # Calculated the intersection point
                cross = shapely.intersection(segment, s).coords[0]
                # Add the point to the list
                cp.append(cross)

        # Add the first point as the start of the line part
        part.append(n1)

        # If the segment crosses an other
        if len(cp) > 0:
            # If the segment starts with a cross point, remove the first node
            # This avoids having a distance smaller than the offset value when it shouldn't
            dpart = part[1:] if breakstart else part
            # Calculate the min distance between the nodes and the original line
            distance = round(__segment_distance(dpart, line), 8)

            # If the distance is above or equal to the offset value (should be equal)
            if distance >= abs(offset):
                # Add the cross point to the group...
                part.append(cp[0])
                # ...and add the group to the list of parts
                parts.append(part)

            # Particular case, the segment crosses multiple segments
            if len(cp) > 1:
                # Recalculate the distance between the line and the segment between the two crossing points.
                # TODO: This might cause issues if the new segment is too close to the line.
                # Maybe, the rounding value can be tweaked for better results.
                cpdist = round(__segment_distance(cp, line), 8)
                if cpdist >= abs(offset):
                    # Add the full list of cross point (should be 2 max) as a full group
                    parts.append(cp)

            # Restart a new group with the start being a cross point
            part = [cp[-1]]
            # Set the cross point start flag to true
            breakstart = True

        # If the last point of the full line is reached, add the last point to the last individual part
        if (n + 1) == len(fline) - 1:
            part.append(n2)

    # Add the last created part to the full list of parts
    parts.append(part)

    # Merge parts that have a common set of coordinates
    groups = __merge_connected(parts)

    return groups

def circle_interpolation(a, b, c, rotation='cw', quad_segs=8):
    """
    Given b and c two points at equal distance from a third point a, interpolates n points between
    b and c along the circle of center a and of radius ab. The number of provided point depends on the
    number of segments allowed per circle quadrant (default to 8).
    Parameters
    ----------
    a : tuple
        Point A coordinates
    b : tuple
        Point B coordinates
    c : tuple
        Point C coordinates
    rotation : str optional
        Define the rotation direction, clockwise ('cw') or counterclockwise ('ccw').
        Default is clockwise.
    quad_segs : int optional
        The number of point allowed on a quarter circle. This defines the number of interpolated points.
        Default to 8.
    """

    # Create vectors
    ab = np.array(b) - np.array(a)
    ac = np.array(c) - np.array(a)

    # Calculate circle radius
    radius = np.linalg.norm(ab)

    # Check that b and c are equidistant from a
    if not np.isclose(np.linalg.norm(ab), np.linalg.norm(ac)):
        raise ValueError("Points b and c are not equidistant from a.")

    # Check if vector ab and ac are equal, if so return b and c
    # This can happen if b and c are really close to each other
    if np.allclose(ab, ac):
        return [b, c]

    # Calculate the full angle formed by ab and ac
    tangle = np.arccos(np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac)))
    
    # Calculate the number of needed point to interpolate between b and c
    # This is based on the quad_segs argument
    n_points = int(tangle / (np.pi / 2) * quad_segs) + 1

    # Start the interpolated point by b
    result = [b]

    # If more than 2 points are needed, it means the interpolation is required
    if n_points > 2:
        # Calculate the individual angle value depending on the rotation direction
        if rotation == 'cw':
            angle = - (tangle / (n_points - 1))
        elif rotation == 'ccw':
            angle = tangle / (n_points - 1)
        else:
            raise ValueError("Rotation value muste be clockwise ('cw') or counterclockwise ('ccw').")

        # Calculate the rotation matrix
        rmatrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

        # Calculated the rotated vector
        rotated = ab / np.linalg.norm(ab)

        # Loop through the number of point needed
        for i in range(1, n_points - 1):
            # Apply the rotation matrix to the vector
            rotated = np.dot(rmatrix, rotated)
            # Calculate the new point coordinates
            interpolated = a + radius * rotated
            # Add the point to the list
            result.append(tuple(interpolated))
    
    # Add c to the result
    result.append(c)

    return result