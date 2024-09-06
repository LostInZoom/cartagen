import numpy as np
import shapely

def dilate_line(line, offset, cap_style='round', quad_segs=8):
    """
    Dilate a line on one side.

    This algorithm proposed by Mustière :footcite:p:`mustiere:2001` dilates
    a line on one side by a given distance in meters. It is the basis of many
    mountain roads generalisation algorithm.
    
    Parameters
    ----------
    line : LineString
        The line to offset.
    offset : float
        The length of the offset to apply in meters. Negative value for left-side dilation, positive for right-side.
    cap_style : str, optional
        The type of caps at the start and end of the provided line. Possible values are 'round' or 'flat'.
    quad_segs : int, optional
        The number of point allowed per circle quadrant when interpolating points using round method.

    Returns
    -------
    list of LineString

    See Also
    --------
    offset_line : 
        This function preserves the relationship between the provided list of points and the result.
    circle_interpolation :
        The function used to interpolate point along the quadrant of a circle.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(1, 1), (2, 3), (5, 2)])
    >>> dilate_line(line, 1.0)
    [<LINESTRING (0.553 0.106, 0.404 0.197, 0.274 0.312, 0.165 0.449, 0.082 0.602...>]
    """

    # Check if provided cap style is allowed
    acaps = ['round', 'flat']
    if cap_style not in acaps:
        raise Exception('Offset curve cap style unknown: {0}.'.format(cap_style))

    # If offset is 0, return the provided line
    if offset == 0:
        return line

    # Calculate the offset points along the line
    groups = offset_line(line, offset, cap_style, quad_segs)

    # Reconstruct the line into parts
    parts, breaks = reconstruct_line(groups, line, offset)

    # Merge parts that have a common set of coordinates
    groups = merge_connected_parts(parts)

    return groups

def offset_line(line, offset, cap_style='round', quad_segs=8):
    """
    Offset the vertices of a line on one side.

    Offset the vertex of the line by a given distance and keeps the relationship
    between the line vertices and the result.

    Parameters
    ----------
    line : LineString
        The line to offset.
    offset : float
        The length of the offset to apply in meters. Negative value for left-side dilation, positive for right-side.
    cap_style : str, optional
        The type of caps at the start and end of the line. Possible values are 'round' or 'flat'.
    quad_segs : int, optional
        The number of point allowed per circle quadrant when interpolating points using round method.

    Returns
    -------
    list of dict
        The index of each dict correpond to the index of the provided points, with the keys:
        
        - *'type'*: The type of point, this value can be:

          - *'start'*: The start of the list. The number of coordinates depends on ``cap_style`` and ``quad_segs``.
          - *'concave'*: The offset vertex is inside a bend.
          - *'convex'*: The offset vertex is outside a bend.
          - *'end'*: The start of the list. Same as **start**.
        
        - *'original'*: The index of the vertex inside the input line
        - *'projected'*: The list of projected coordinates

    See Also
    --------
    dilate_line : 
        This function loses the relationship between the provided line vertex and the offset points.
    circle_interpolation :
        The function used to interpolate point along the quadrant of a circle.
    """
    def calculate_extremity(x1, y1, x2, y2):
        # Calculate the rounded extremity of the offset curve
        ab = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        c = (x2 - (((x2 - x1) / ab) * (ab + abs(offset))), y2 - (((y2 - y1) / ab) * (ab + abs(offset))))
        return c

    points = list(line.coords)

    # Storage for the projected line
    pline = []

    # Loop through vertex of the original line
    for i, n in enumerate(points):
        # Get the current point
        a = n

        # If it's the last point of the line
        if i >= (len(points) - 1):
            # If round cap style is chosen
            if cap_style == 'round':
                # Calculate the rounded extremity of the offset curve
                c = calculate_extremity(xb, yb, xa, ya)

                # Interpolate points between the last projected point and the rounded extremity c
                rotation = 'ccw' if offset < 0 else 'cw'

                # Replace the last point with the interpolation
                pline[-1]['type'] = 'end'
                pline[-1]['projected'] = circle_interpolation(a, pline[-1]['projected'][-1], c, rotation=rotation, quad_segs=quad_segs)
            
            # Break the loop as it's the last point
            # Do nothing for flat cap style
            break

        # Get the next point
        b = points[i + 1]
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
        nbv = - offset * proj / np.linalg.norm(proj)

        # Create the two new points
        a1 = tuple(a + nbv)
        b1 = tuple(b + nbv)

        # If it's the first node
        if i == 0:
            # If round cap style is chosen
            if cap_style == 'round':
                # Calculate the rounded extremity of the offset curve
                c = calculate_extremity(xa, ya, xb, yb)

                # Interpolate points between the projected a1 point and the rounded extremity c
                rotation = 'ccw' if offset < 0 else 'cw'
                ip = circle_interpolation(a, c, a1, rotation=rotation, quad_segs=quad_segs)

                # Add the interpolated points as the first entry
                p = ip

            # If a flat cap style is chosen
            else:
                p = [a1]

            pline.append({
                'type': 'start',
                'original': i,
                'projected': p
            })

        # If it's a point between the first and the last (middle point)
        else:
            # TODO: Better handle concavities on the line, for now, when the projected segment
            # crosses the original line, some weird things happen, but these are removed after the line reconstruction
            # Retrieve the last projected point a and b
            a0, b0 = pline[-2]['projected'][-1], pline[-1]['projected'][-1]

            # Create the previous and the current segment
            seg0 = shapely.LineString([a0, b0])
            seg1 = shapely.LineString([a1, b1])

            # Check if the current and previous segments are crossing
            if seg0.crosses(seg1):
                # If they are, it means the angle is concave, add both points
                t = 'concave'
                p = [b0, a1]

            # If they are not           
            else:
                # Making a circle interpolation between those two points.
                rotation = 'ccw' if offset < 0 else 'cw'
                t = 'convex'
                p = circle_interpolation(a, b0, a1, rotation=rotation, quad_segs=quad_segs)

            # Replace last point with this one
            pline[-1]['type'] = t
            pline[-1]['original'] = i
            pline[-1]['projected'] = p

        # Add the projection of b to the list
        pline.append({
            'type': 'end',
            'original': i + 1,
            'projected': [b1]
        })

    return pline

def circle_interpolation(a, b, c, rotation='cw', quad_segs=8):
    """
    Interpolate points along a circle quadrant.

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
    rotation : str, optional
        Define the rotation direction, clockwise ('cw') or counterclockwise ('ccw').
        Counter-clockwise interpolation is the plane symetry of the clockwise interpolation
        using the line ``ab``.
    quad_segs : int, optional
        The number of point allowed on a quarter circle.
        This defines the number of interpolated points.

    Returns
    -------
    list of tuple :
        The list of coordinates of the interpolated points.

    Examples
    --------
    >>> a = (1, 1)
    >>> b = (2, 1)
    >>> c = (1, 2)
    >>> circle_interpolation(a, b, c, quad_segs=2)
    [(2, 1), (1.985, 0.826), (1.94, 0.658), (1.866, 0.5), (1.766, 0.357), (1.643, 0.234), (1.5, 0.134), (1.342, 0.06), (1.174, 0.015), (1, 2)]
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
    n_points = int(tangle / (np.pi / 2) * quad_segs) + 2

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

def merge_connected_parts(groups):
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

    return [shapely.LineString(c) for c in result]

def reconstruct_line(groups, line, offset):
    """
    Reconstruct the projected line from the offset of points.
    Return the line parts and the list of breakpoints.
    """

    # Return the minimum distance between the line and the coordinates provided.
    def __segment_distance(coords, line):
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

    # Create individual segments from the full nested list of points.
    def __create_segments(points):
        segments = []

        nn = None
        # Loop through groups of points
        for igroup, group in enumerate(points):
            if nn is not None:
                segments.append({
                    'original': [igroup - 1, igroup],
                    'geometry': shapely.LineString([nn, group['projected'][0]])
                })

            gn = None
            # Loop through each point of the group
            for inode, node in enumerate(group['projected']):
                # Add the segment if the previous point exists
                if gn is not None:
                    segments.append({
                        'original': [igroup],
                        'geometry': shapely.LineString([gn, node])
                    })

                # Set the previous point as the current
                gn = node

            nn = node

        return segments

    # Handle concavities inside the line by calculating the intersection formed by the previous and following lines
    def __handle_concavity(groups):
        # Storage to update geometries
        geom = []

        # Loop through groups
        for i, group in enumerate(groups):
            # Check that group is a concavity
            if group['type'] == 'concave':
                # Calculate previous and following segment
                seg1 = shapely.LineString([groups[i - 1]['projected'][-1], group['projected'][0]])
                seg2 = shapely.LineString([group['projected'][-1], groups[i + 1]['projected'][0]])
                
                # Calculate intersection between both segments
                intersection = shapely.intersection(seg1, seg2)

                # If the intersection is a Point, replace the projected points of the group
                if intersection.geom_type == 'Point':
                    geom.append((i, intersection.coords[0]))

        # Loop through the new geoms
        for g in geom:
            # Update the geometry of the corresponding group
            groups[g[0]]['projected'] = [g[1]]

        return groups

    # Handle concavities
    groups = __handle_concavity(groups)

    # Get the full line as single segments
    segments = __create_segments(groups)

    # Flag to see if the current line starts with a point where segment crosses
    breakstart = False
    # Storage for the break points
    breakpoints = []
    breaks = []
    # Storage for the final groups of individual parts (below)
    parts = []
    # Storage for individual groups of nodes forming a continuous line
    part = []

    # Loop through all vertex of the full line
    for i, group in enumerate(groups):
        nodes = group['projected']
        for j, node in enumerate(nodes):
            # Get the current node
            n1 = node
            # If it's the last node of the group
            if j == (len(nodes) - 1):
                # If it's not the last group of the line
                if i < (len(groups) - 1):
                    # Get the first point of the next group of nodes
                    n2 = groups[i + 1]['projected'][0]
                    o1 = [i, i + 1]
            else:
                # Get the next point of the group
                n2 = nodes[j + 1]
                o1 = [i]

            # Here, n1 and n2 are assigned, creating the segment
            segment = shapely.LineString([n1, n2])

            # Storage for points where segments crosses
            cp = []
            # Loop through segments
            for s in segments:
                # If the current segment crosses an other
                if segment.crosses(s['geometry']):
                    # Calculate the intersection point
                    cross = shapely.intersection(segment, s['geometry']).coords[0]

                    # Get associated original segment
                    o2 = s['original']

                    # Create the object to add
                    obj = {'o1': o1, 'o2': o2, 'geometry': cross}

                    # If a crossing point already exists
                    if len(cp) > 0:
                        # Create previous and current crossing point
                        pp = shapely.Point(cp[0]['geometry'])
                        pc = shapely.Point(cross)
                        # Create starting point of the projected segment
                        p1 = shapely.Point(n1)

                        # If the current crossing point if closer from the segment start
                        if shapely.distance(pp, p1) > shapely.distance(pc, p1):
                            # Insert the new crossing point at the start of the list to preserve the order
                            cp.insert(0, obj)
                        else:
                            # Add the point at the end
                            cp.append(obj)
                    else:
                        # Append the point to the list
                        cp.append(obj)

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
                    part.append(cp[0]['geometry'])
                    # ...and add the group to the list of parts
                    parts.append(part)

                    # Add the break point to the list if it doesn't already exists
                    if cp[0]['geometry'] not in breakpoints:
                        breakpoints.append(cp[0]['geometry'])
                        breaks.append(cp[0])

                # Particular case, the segment crosses multiple segments
                if len(cp) > 1:
                    # Recalculate the distance between the line and the segment between the two crossing points.
                    # TODO: This might cause issues if the new segment is too close to the line.
                    # Maybe, the rounding value can be tweaked for better results.
                    geomlist = [x['geometry'] for x in cp]
                    cpdist = round(__segment_distance(geomlist, line), 8)
                    if cpdist >= abs(offset):
                        # Add the full list of cross point as a full group
                        parts.append(geomlist)

                        # Add all breaks if they do not already exists
                        for b in cp:
                            if b['geometry'] not in breakpoints:
                                breakpoints.append(b['geometry'])
                                breaks.append(b)

                # Restart a new group with the start being a cross point
                part = [cp[-1]['geometry']]
                # Set the cross point start flag to true
                breakstart = True

            # If the last point of the full line is reached, add the last point to the last individual part
            if (i == (len(groups) - 1)) and (j == (len(nodes))):
                part.append(n2)

    # Add the last created part to the full list of parts only if its last point is not too close to the line
    if round(shapely.distance(shapely.Point(part[-1]), line), 8) >= abs(offset):
        parts.append(part)

    return parts, breaks