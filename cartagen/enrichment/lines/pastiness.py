import shapely
from shapely.ops import nearest_points

from cartagen.utils.geometry.dilation import offset_line, merge_connected_parts, reconstruct_line

def detect_pastiness(line, tolerance, cap_style='round', quad_segs=8):
    """
    Detect the pastiness of a series of bends.

    This algorithm proposed by Mustière :footcite:p:`mustiere:2001`
    subdivide the provided line into multiple chunks, thus modifying the geometry,
    it is not a data enrichment function stricto sensu.

    Parameters
    ----------
    line : LineString
        The line to detect the pastiness from.
    tolerance : float
        The width of the offset used to detect the pastiness.
    cap_style : str, optional
        The type of caps at the start and end of the provided line.
        Possible values are 'round' or 'flat'.
    quad_segs : int, optional
        The number of point allowed per circle quadrant
        when interpolating points using round method.

    Returns
    -------
    list of dict
        The line subdivided into chunks depending of the paste.
        
        Dict keys are as following:

        - *'paste'* represents the number of conflicts, it can be:

          * *0* when no conflicts are detected.
          * *1* when a conflict exists on one side only.
          * *2* when conflicts are on both side of the line).

        - *'geometry'* is the geometry of the line section.

    References
    ----------
    .. footbibliography::
    """
    # Handle individual sides
    def __treat_side(line, tolerance, cap_style, quad_segs):
        # Prepare the break points
        def __prepare_breaks(breaks, coords, lines, oline):
            # Project the point on a segment or a vertex
            def __project_point(o, coords, oline):
                # If the length is 1
                if len(o) == 1:
                    # The point is the node
                    p = coords[o[0]]
                    point = shapely.Point(p)
                # If it's not (must be 2)
                else:
                    # Create the segment formed by both coordinates
                    segment = shapely.LineString([coords[o[0]], coords[o[1]]])
                    # Project the point on the segment
                    point = nearest_points(shapely.Point(b['geometry']), segment)[1]
                    p = point.coords[0]

                if oline is not None:
                    return p, oline.project(point)
                else:
                    return p

            # Loop through break points
            for b in breaks:
                # Invert o1 and o2 which can be in the wrong order when the offset line
                # was crossing itself with a section of itself closer to its start
                if b['o1'][-1] > b['o2'][0]:
                    b['o1'], b['o2'] = b['o2'], b['o1']

                # Retrieve breakpoint geometry
                geom = shapely.Point(b['geometry'])
                # Loop through group lines
                for iline, gline in enumerate(lines):
                    # If the breakpoint intersects the line
                    if shapely.intersects(geom, gline['geometry']):
                        # Add the line index as the group of the breakpoint
                        b['group'] = iline

                # Project the point on both associated segments or vertex
                b['p1'], b['length'] = __project_point(b['o1'], coords, oline)
                b['p2'] = __project_point(b['o2'], coords, None)

            # Order the breaks using the length from the start of the original line to the projection of the break on the line
            # This is useful when dealing with multiple breaks projected on the same segment.
            return sorted(breaks, key=lambda d: d['length'])

        # Create a conflict
        def __create_conflict(breaks, i1, i2, coords, ctype, otype):
            # Retrieve the index of the first concerned node of the main line
            n = breaks[i1]['o{0}'.format(i2)][0]

            # Create the list entry
            e = { 'type': ctype, 'node': n, 'distance': 0, 'order': otype }

            # Get coordinates of the projection of the node on the line and the associated vertex
            c, cl = breaks[i1]['p{0}'.format(i2)], coords[n]

            # If the projected point is different from the vertex
            if c != cl:
                # Add the distance as an attribute of the dict
                e['distance'] = shapely.distance(shapely.Point(cl), shapely.Point(c))
                e['coords'] = c
            
            return e

        coords = list(line.coords)

        # Calculate the offset points along the line
        groups = offset_line(line, tolerance, cap_style, quad_segs)

        # Reconstruct the line into parts
        parts, breaks = reconstruct_line(groups, line, tolerance)

        # Merge parts that have a common set of coordinates
        groups = merge_connected_parts(parts)

        # Prepare lines
        lines = __prepare_lines(groups)

        # Prepare the break points
        breaks = __prepare_breaks(breaks, coords, lines, line)

        # This flag handles projected points right after a hole conflict
        # and skips the point processing
        skip = False
        # This storage will contain pairs of indexes of projected line parts
        # between which a hole conflict exists. This allows to avoid treating break points
        # that are not hole conflicts.
        done = []

        # Storage for holes and pastes
        conflicts = []

        # Loop through breaking points
        # This version handles nested holes formed by a dilatation
        for i, b in enumerate(breaks):
            if skip:
                skip = False
                continue

            # Retrieve the location of the break (on main line or inside a hole)
            location1 = lines[breaks[i]['group']]['type']

            # If it's not the last break
            if i < (len(breaks) - 1):
                if breaks[i + 1]['group'] != b['group']:
                    # If the index of the group of the following breakpoint and of the current breakpoint
                    # is in the done list, it means that the breaks should not be treated
                    # Note: This should handle fine nested holes formed by the dilatation of a complex weaving mountain roads
                    if [breaks[i + 1]['group'], b['group']] in done:
                        pass
                    else:
                        # Here, there is a conflict of hole inside the symbol
                        # Adding the group indicator to the done list
                        done.append([b['group'], breaks[i + 1]['group']])

                        # Create the first conflict and add it to the list
                        h1 = __create_conflict(breaks, i, 1, coords, 'hole', 'start')
                        h2 = __create_conflict(breaks, i + 1, 1, coords, 'hole', 'end')
                        conflicts.extend([h1, h2])

                        # Create the second conflict and add it to the list
                        h3 = __create_conflict(breaks, i + 1, 2, coords, 'hole', 'start')
                        h4 = __create_conflict(breaks, i, 2, coords, 'hole', 'end')
                        conflicts.extend([h3, h4])

                        skip = True
                        continue

            # Create the break point
            bgeom = shapely.Point(b['geometry'])

            n1, n2 = b['o1'], b['o2']

            # Define a max distance to zero
            maxdist = 0
            # Loop through the nodes of the concerned section of the line
            for index in range(n1[-1], n2[0] + 1):
                # Create the line node point
                ogeom = shapely.Point(coords[index])
                # Calculate the distance between the breakpoint and the line node
                d = shapely.distance(bgeom, ogeom)
                # If the distance is above the max distance, update it
                if d > maxdist:
                    maxdist = d

            # If the distance is above 1.7 times the width of the symbol
            if maxdist > (1.7 * abs(tolerance)):
                # Here, it is a strict pastiness conflict
                p1 = __create_conflict(breaks, i, 1, coords, 'paste', 'start')
                p2 = __create_conflict(breaks, i, 2, coords, 'paste', 'end')

                conflicts.extend([p1, p2])

        return conflicts

    # Prepare the lines generated from the dilatation
    def __prepare_lines(groups):
        lines = []
        # Loops through groups of nodes
        for group in groups:
            coords = list(group.coords)

            # Create individual lines
            line = []
            for node in coords:
                line.append(node)
            
            # If the line has at least one node
            if len(line) > 0:
                # Set the type of the line, main being the main dilatation line
                ltype = 'main'
                # If the first and last node of the line is the same, it's a hole
                if line[0] == line[-1]:
                    ltype = 'hole'
                
                # Add the line to the list
                lines.append({
                    'type': ltype,
                    'geometry': shapely.LineString(line)
                })
        
        return lines

    # Calculate conflicts on both sides of the line
    lconflicts = __treat_side(line, tolerance, cap_style, quad_segs)
    rconflicts = __treat_side(line, -tolerance, cap_style, quad_segs)

    # Sort the conflicts by node number and distance from the node
    conflicts = sorted(rconflicts + lconflicts, key=lambda s: (s['node'], s['distance']))

    # Retrieve coordinates of the original line
    coords = list(line.coords)

    # Storage for the final lines
    result = []
    # The state flag identify the current pastiness level
    state = 0

    # This will be the first line
    current = []

    # Loop through the coordinates of the original line
    for i, c in enumerate(coords):
        # Add the current node of the original line
        current.append(c)

        for o in conflicts:
            if o['node'] == i:
                otype = o['order']

                if 'coords' in o:
                    current.append(o['coords'])
                    result.append({ 'paste': state, 'geometry': shapely.LineString(current) })
                    n = o['coords']
                else:
                    if len(current) > 0:
                        result.append({ 'paste': state, 'geometry': shapely.LineString(current) })
                    n = c
                
                current = [n]

                if otype == 'start':
                    state += 1
                else:
                    state -= 1

    result.append({ 'paste': state, 'geometry': shapely.LineString(current) })

    return result