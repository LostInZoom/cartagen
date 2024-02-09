import shapely
import geopandas as gpd
import pprint

from cartagen4py.utils.geometry.dilation import *
from cartagen4py.utils.geometry.line import *

from test_functions import *

def detect_pastiness(line, width, cap_style='flat', quad_segs=8):
    """
    Detect pastiness of a line object. The result is a list of sections of the original line.
    """

    # Handle individual side
    def __treat_side(line, width, cap_style, quad_segs):
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
            return sorted(breaks, key=lambda d: d['length'])

        # Create a conflict
        def __create_conflict(breaks, i1, i2, coords, ctype):
            # Retrieve the index of the first concerned node of the main line
            n = breaks[i1]['o{0}'.format(i2)][0]

            # Create the list entry
            e = { 'type': ctype, 'node': n }

            # Get coordinates of the projection of the node on the line and the associated vertex
            c, cl = breaks[i1]['p{0}'.format(i2)], coords[n]

            # If the projected point is different from the vertex
            if c != cl:
                # Add the distance as an attribute of the dict
                e['distance'] = shapely.distance(shapely.Point(cl), shapely.Point(c))
                e['coords'] = c
            
            return e

        # Get the coordinates of the line vertex
        coords = list(line.coords)

        # Calculate the offset points along the line
        groups = offset_points(coords, width, cap_style, quad_segs)

        # Reconstruct the line into parts
        parts, breaks = reconstruct_line(groups, line, width)

        # Merge parts that have a common set of coordinates
        groups = merge_connected_parts(parts)

        # Prepare lines
        lines = __prepare_lines(groups)

        # Prepare the break points
        breaks = __prepare_breaks(breaks, coords, lines, line)

        # make_gdf(breaks, 'breaks')

        # This flag handles projected points right after a hole conflict
        # and skips the point processing
        skip = False
        # This storage will contain pairs of indexes of projected line parts
        # between which a hole conflict exists. This allows the avoid treating break points
        # that are not hole conflict.
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

                        # Create the conflict for the first part of the hole conflict
                        h1 = __create_conflict(breaks, i, 1, coords, 'hole')
                        h4 = __create_conflict(breaks, i, 2, coords, 'hole')

                        # Create the conflict for the second part of the hole conflict
                        h2 = __create_conflict(breaks, i + 1, 1, coords, 'hole')
                        h3 = __create_conflict(breaks, i + 1, 2, coords, 'hole')

                        conflicts.extend([(h1, h2), (h3, h4)])

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
            if maxdist > (1.7 * abs(width)):
                # Here, it is a strict pastiness conflict
                p1 = __create_conflict(breaks, i, 1, coords, 'paste')
                p2 = __create_conflict(breaks, i, 2, coords, 'paste')

                conflicts.append((p1, p2))

        return sorted(conflicts, key=lambda d: d[0]['node'])

    # Prepare the lines generated from the dilatation
    def __prepare_lines(groups):
        lines = []
        # Loops through groups of nodes
        for group in groups:
            # Create individual lines
            line = []
            for node in group:
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
    lconflicts = __treat_side(line, width, cap_style, quad_segs)
    rconflicts = __treat_side(line, -width, cap_style, quad_segs)

    # Retrieve coordinates of the original line
    coords = list(line.coords)

    # Copy the list of coordinates of the original line
    newcoords = coords.copy()

    result = []

    # Loop through conflicts
    for conflict in rconflicts + lconflicts:       
        # Retrieve start and end object
        start, end = conflict

        # Retrieve the type of conflict
        ctype = start['type']

        # Retrieve the start and end node index
        n1, n2 = start['node'], end['node']

        nodes = []
        # Loop through the original line coordinates
        for i, c in enumerate(coords):
            # Append the index if it is between the start and end node index
            if n1 <= i <= n2:
                nodes.append(c)
        
        if 'distance' in start:
            nodes[0] = start['coords']

        if 'distance' in end:
            nodes.append(end['coords'])

        result.append({'geometry': shapely.LineString(nodes)})

    make_gdf(result, 'sections')

    return None