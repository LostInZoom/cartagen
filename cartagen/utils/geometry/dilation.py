import numpy as np
import shapely
from shapely import Point, LineString
from shapely.ops import unary_union

from cartagen.utils.partitioning.network import network_faces

def dilate_line(line, offset, cap_style='round', quad_segs=8):
    """
    Dilate a line on both sides.

    This algorithm proposed by Mustière :footcite:p:`mustiere:2001-a` dilates
    a line on both sides by a given distance in meters. It is the basis of many
    mountain roads generalisation algorithm.
    
    Parameters
    ----------
    line : LineString
        The line to offset.
    offset : float
        The length of the offset to apply in meters.
    cap_style : str, optional
        The type of caps at the start and end of the provided line. Possible values are 'round' or 'flat'.
    quad_segs : int, optional
        The number of point allowed per circle quadrant when interpolating points using round method.

    Returns
    -------
    left, right: tuple of list of LineString
        A tuple of two list of LineString, the left and the right side of the dilation.

    See Also
    --------
    offset_line : 
        This function preserves the relationship between the provided list of points and the result.
    circle_interpolation :
        The function used to interpolate point along the quadrant of a circle.

    Notes
    -----
    This algorithm returns dilation as lists of LineString. It often is composed of only one line but
    sometimes, it can creates multiple lines when the dilation creates holes.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(1, 1), (2, 3), (5, 2)])
    >>> dilate_line(line, 1.0)
    ([<LINESTRING (0.553 0.106, 0.715 0.042, 0.886 0.007, 1.06 0.002, 1.232 0.027,...>], [<LINESTRING (0.553 0.106, 0.404 0.197, 0.274 0.312, 0.165 0.449, 0.082 0.602...>])
    """

    # Check if provided cap style is allowed
    acaps = ['round', 'flat']
    if cap_style not in acaps:
        raise Exception('Offset curve cap style unknown: {0}.'.format(cap_style))

    # If offset is 0, return the provided line
    if offset == 0:
        return line

    # Calculate the offset points along the line
    groups1 = offset_line(line, offset, cap_style, quad_segs)
    groups2 = offset_line(line, -offset, cap_style, quad_segs)

    # Reconstruct the line into parts
    left, right = reconstruct_line(groups1, groups2, line, offset)

    return left[0], right[0]

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
        ab = float(np.sqrt((y2 - y1)**2 + (x2 - x1)**2))
        c = (float(x2 - (((x2 - x1) / ab) * (ab + abs(offset)))), float(y2 - (((y2 - y1) / ab) * (ab + abs(offset)))))
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
                interpolation = circle_interpolation(a, pline[-1]['projected'][-1], c, rotation=rotation, quad_segs=quad_segs)
                pline[-1]['projected'] = [ shapely.set_precision(Point(x), 0.000000001).coords[0] for x in interpolation ]
            
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

        a1 = (float(a1[0]), float(a1[1]))
        b1 = (float(b1[0]), float(b1[1]))

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
                p = [ shapely.set_precision(Point(x), 0.000000001).coords[0] for x in ip ]

            # If a flat cap style is chosen
            else:
                p = [ shapely.set_precision(Point(a1), 0.000000001).coords[0] ]

            pline.append({
                'type': 'start',
                'original': i,
                'projected': p
            })

        # If it's a point between the first and the last (middle point)
        else:
            # Retrieve the last projected point a and b
            a0, b0 = pline[-2]['projected'][-1], pline[-1]['projected'][-1]
            a0 = (float(a0[0]), float(a0[1]))
            b0 = (float(b0[0]), float(b0[1]))

            # Create the previous and the current segment
            seg0 = shapely.LineString([a0, b0])
            seg1 = shapely.LineString([a1, b1])

            # Check if the current and previous segments are crossing
            if seg0.crosses(seg1):
                # If they are, it means the angle is concave, add both points
                t = 'concave'
                p = [ shapely.set_precision(Point(b0), 0.000000001).coords[0], shapely.set_precision(Point(a1), 0.000000001).coords[0] ]

            # If they are not crossing        
            else:
                # Making a circle interpolation between those two points.
                rotation = 'ccw' if offset < 0 else 'cw'
                t = 'convex'
                interpolation = circle_interpolation(a, b0, a1, rotation=rotation, quad_segs=quad_segs)
                p = [ shapely.set_precision(Point(x), 0.000000001).coords[0] for x in interpolation ]

            # Replace last point with this one
            pline[-1]['type'] = t
            pline[-1]['original'] = i
            pline[-1]['projected'] = p

        # Add the projection of b to the list
        pline.append({
            'type': 'end',
            'original': i + 1,
            'projected': [ shapely.set_precision(Point(b1), 0.000000001).coords[0] ]
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
            result.append((float(interpolated[0]), float(interpolated[1])))
    
    # Add c to the result
    result.append(c)

    return result

def reconstruct_line(groups1, groups2, line, offset):
    """
    Reconstruct lines using provided groups of dilated points.
    """
    def __polygonize(line, groups, offset):
        # Get the list of vertex of the original line 
        coords = list(line.coords)

        # Reverse the line coordinates
        reverse = [ Point(c) for c in coords][::-1]

        # Get the vertexes of the original line along with the projected line
        vertexes = reverse.copy()

        for group in groups:
            vertexes.extend([ Point(p) for p in group['projected'] ])
        vertexes.append(reverse[0])
        
        # Create the individual segments from the vertexes
        segments = []
        for p in range(len(vertexes) - 1):
            segments.append(shapely.LineString([vertexes[p], vertexes[p + 1]]))

        # Polygonize the "network"
        polygons = network_faces(segments, convex_hull=False)

        # Retrieve the holes formed by the projection if they exists
        # by retrieving polygon that have their interior further from the line than the offset value
        holes = []
        for polygon in polygons:
            if shapely.distance(polygon.point_on_surface(), line) >= abs(offset):
                if offset > 0:
                    holes.append(list(polygon.exterior.coords)[::-1])
                else:
                    holes.append(list(polygon.exterior.coords))

        # Merge all the polygons to form one
        polygon = unary_union(polygons)

        # Sometimes, a MultiPolygon can be created and create issues
        # Keeping only the largest sub polygon
        if polygon.geom_type == 'MultiPolygon':
            largest = None
            area = 0
            for sub in polygon.geoms:
                if sub.area > area:
                    area = sub.area
                    largest = sub
            polygon = largest
        
        # Get the boundary of the merged polygon as a list of coords and remove the last vertex (same as first)
        return list(polygon.boundary.coords)[:-1], holes

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

                if shapely.crosses(seg1, seg2):
                    geom.append((i, shapely.intersection(seg1, seg2).coords[0]))

        # Loop through the new geoms
        for g in geom:
            # Update the geometry of the corresponding group
            groups[g[0]]['projected'] = [g[1]]

        return groups

    def __handle_wrapping(groups, polygon, iprojected):
        # Get the start and end vertex of the projected line
        start, end = groups[0]['projected'][0], groups[-1]['projected'][-1]

        first, last = None, None
        # Check if the start of the projected line is in the polygon boundary
        if start in polygon:
            first = polygon.index(start)
        else:
            for igroup, inode in iprojected:
                g = groups[igroup]['projected'][inode]
                if g in polygon:
                    first = polygon.index(g)
                    break 

        # Same check for the end of the projected line
        if end in polygon:
            last = polygon.index(end)
        else:
            for igroup, inode in iprojected[::-1]:
                g = groups[igroup]['projected'][inode]
                if g in polygon:
                    last = polygon.index(g)
                    break 

        return first, last

    def __treat_part(coords, breaks, groups, iprojected, hole):
        if len(breaks) == 0:
            return coords, []

        breakpoints = []
        part = []

        previous_node = None
        crosspoint = None

        minindex = len(iprojected)
        mincoords = 0

        if not hole:
            if breaks[0] == 0:
                breaks.pop(0)

            if breaks[-1] == (len(coords) - 1):
                breaks.pop()

        for index in breaks:
            b = Point(coords[index])
            
            s1, s2 = None, None
            o1, o2 = None, None
            ip1 = None
            for ip in range(len(iprojected) - 1):
                igroup = iprojected[ip][0]
                inode = iprojected[ip][1]
                current = groups[igroup]['projected'][inode]

                if inode == len(groups[igroup]['projected']) - 1:
                    igroup2, inode2 = igroup + 1, 0
                    g = [ igroup, igroup + 1 ]
                else:
                    igroup2, inode2 = igroup, inode + 1
                    g = [ igroup ]
                
                following = groups[igroup2]['projected'][inode2]

                segment = {
                    'n1': {
                        'igroup': igroup,
                        'inode': inode,
                        'geometry': current
                    },
                    'n2': {
                        'igroup': igroup2,
                        'inode': inode2,
                        'geometry': following
                    },
                    'geometry': shapely.LineString([current, following])
                }

                if shapely.intersects(segment['geometry'], b.buffer(0.00000001)):
                    if s1 is None:
                        s1 = segment
                        o1 = g
                        ip1 = ip
                    else:
                        if s2 is None:
                            s2 = segment
                            o2 = g
            
            if ip1 < minindex:
                minindex = ip1
                mincoords = index

            crosspoint = shapely.intersection(s1['geometry'], s2['geometry']).coords[0]
            breakpoints.append({ 's1': s1, 's2':s2, 'o1': o1, 'o2': o2, 'geometry': crosspoint })

            if hole:
                coords[index] = crosspoint
            else:
                if s1['n1']['geometry'] in coords:
                    start = 0 if previous_node is None else previous_node
                    end = coords.index(s1['n1']['geometry'])
                    part.extend(coords[start:end + 1])
                part.append(crosspoint)
                
                if s2['n2']['geometry'] in coords:
                    previous_node = coords.index(s2['n2']['geometry'])

        if hole:
            part = coords[mincoords:] + coords[0:mincoords] + [coords[mincoords]]
            break_index = breaks.index(mincoords)
            breakpoints = breakpoints[break_index:] + breakpoints[0:break_index]
        else:
            part.extend(coords[previous_node:])

        return part, breakpoints
    
    # Remove small concavities from the projected line
    groups1 = __handle_concavity(groups1)
    groups2 = __handle_concavity(groups2)

    # Polygonize both sides of the dilation
    polygon1, holes1 = __polygonize(line, groups1, offset)
    polygon2, holes2 = __polygonize(line, groups2, offset)

    # Retrieve the full list of projected points
    iprojected1 = [ (igroup, inode) for igroup, group in enumerate(groups1) for inode, node in enumerate(group['projected']) ]
    iprojected2 = [ (igroup, inode) for igroup, group in enumerate(groups2) for inode, node in enumerate(group['projected']) ]

    first1, last1 = __handle_wrapping(groups1, polygon1, iprojected1)
    last2, first2 = __handle_wrapping(groups2, polygon2, iprojected2)
    
    # Build the dilated line by getting the merged polygon boundary between start and end vertex
    if first1 > last1:
        dilated1 = polygon1[first1:] + polygon1[:last1 + 1]
    else:
        dilated1 = polygon1[first1:last1 + 1]

    # Build the dilated line by getting the merged polygon boundary between start and end vertex
    if first2 > last2:
        dilated2 = polygon2[first2:] + polygon2[:last2 + 1]
    else:
        dilated2 = polygon2[first2:last2 + 1]
    dilated2.reverse()

    line1 = LineString(dilated1)
    line2 = LineString(dilated2)

    # Check if the two sides crosses
    if shapely.crosses(line1, line2):
        index1, index2, intersection = None, None, None
        for i1 in range(len(dilated1) - 1):
            s1 = LineString([dilated1[i1], dilated1[i1 + 1]])
            if shapely.crosses(s1, line2):
                for i2 in range(len(dilated2) - 1):
                    s2 = LineString([dilated2[i2], dilated2[i2 + 1]])
                    if shapely.crosses(s1, s2):
                        index1, index2, intersection = i1 + 1, i2 + 1, shapely.intersection(s1, s2)
        
        dilated1 = [intersection] + dilated1[index1:]
        dilated2 = [intersection] + dilated2[index2:]

    line1 = LineString(dilated1)
    line2 = LineString(dilated2)

    # Retrieve the list of break points,
    # i.e. points that are present on the dilated line
    # but absent from the projected points (intersection points)
    projected1 = [ p for group in groups1 for p in group['projected'] ]
    projected2 = [ p for group in groups2 for p in group['projected'] ]
    breakpoints1 = [ i for i, p in enumerate(dilated1) if p not in projected1 ]
    breakpoints2 = [ i for i, p in enumerate(dilated2) if p not in projected2 ]

    # from cartagen.utils.debug import plot_debug
    # plot_debug(line, LineString(projected1), LineString(projected2), [Point(dilated1[x]) for x in breakpoints1])

    part1, breaks1 = __treat_part(dilated1, breakpoints1, groups1, iprojected1, False)
    part2, breaks2 = __treat_part(dilated2, breakpoints2, groups2, iprojected2, False)

    parts1 = [ part1 ]
    parts2 = [ part2 ]

    for hole1 in holes1:
        hole1.pop()
        hbreakpoints = [ i for i, h in enumerate(hole1) if h not in projected1 ]
        hpart, hbreaks = __treat_part(hole1, hbreakpoints, groups1, iprojected1, True)
        parts1.append(hpart)
        breaks1.extend(hbreaks)

    for hole2 in holes2:
        hole2.pop()
        hbreakpoints = [ i for i, h in enumerate(hole2) if h not in projected2 ]
        hpart, hbreaks = __treat_part(hole2, hbreakpoints, groups2, iprojected2, True)
        parts2.append(hpart)
        breaks2.extend(hbreaks)
    
    return ([ shapely.LineString(x) for x in parts2 ], breaks2), ([ shapely.LineString(x) for x in parts1 ], breaks1)