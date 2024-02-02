import shapely
import geopandas as gpd

from cartagen4py.utils.geometry.dilation import *
from cartagen4py.utils.geometry.line import *

def detect_pastiness(line, width, cap_style='flat', quad_segs=8):
    """
    Detect pastiness of a line object. The result is a list of sections of the original line.
    """

    # Prepare the break points
    def __prepare_breaks(breaks, coords, lines): 
         # Project the point on a segment or a vertex
        def __project_point(o, coords):
            # If the length is 1
            if len(o) == 1:
                # The point if the node
                p = coords[o[0]]
            # If it's not (must be 2)
            else:
                # Create the segment formed by both coordinates
                segment = shapely.LineString([coords[o[0]], coords[o[1]]])
                # Project the point on the segment
                p = nearest_points(shapely.Point(b['geometry']), segment)[1].coords[0]
            return p
        
        # Loop through break points
        for b in breaks:
            # Retrieve breakpoint geometry
            geom = shapely.Point(b['geometry'])
            # Loop through group lines
            for iline, gline in enumerate(lines):
                # If the breakpoint intersects the line
                if shapely.intersects(geom, gline):
                    # Add the line index as the group of the breakpoint
                    b['group'] = iline

            # Project the point on both associated segments or vertex
            b['p1'] = __project_point(b['o1'], coords)
            b['p2'] = __project_point(b['o2'], coords)

        return breaks

    def __get_max_distance(point, line):
        pass

    coords = list(line.coords)

    # Calculate the offset points along the line
    groups = offset_points(coords, width, cap_style, quad_segs)

    # Reconstruct the line into parts
    parts, breaks = reconstruct_line(groups, line, width)

    # Merge parts that have a common set of coordinates
    groups = merge_connected_parts(parts)

    glines = []
    for nodes in groups:
        gline = []
        for node in nodes:
            gline.append(node)
        if len(gline) > 0:
            glines.append(shapely.LineString(gline))

    breaks = __prepare_breaks(breaks, coords, glines)

    # This will store the new sections
    sections = []

    inhole = False
    skip = False
    # Loop through breaking points
    # TODO: Allow having nested holes within the symbol
    for i, b in enumerate(breaks):
        # If it's not the last break
        if i < (len(breaks) - 1):
            if skip:
                skip = False
                continue

            # If this break is not on the same group as the following
            if breaks[i + 1]['group'] != b['group'] and inhole == False:
                # Here, there is a conflict of hole inside the symbol
                p1 = shapely.Point(b['p1'])
                p2 = shapely.Point(breaks[i + 1]['p1'])

                middle = split_line_at_point(split_line_at_point(line, p1)[1], p2)[0]
                sections.append({
                    'type': 'hole',
                    'geometry': middle    
                })

                p3 = shapely.Point(b['p2'])
                p4 = shapely.Point(breaks[i + 1]['p2'])

                middle = split_line_at_point(split_line_at_point(line, p4)[1], p3)[0]
                sections.append({
                    'type': 'hole',
                    'geometry': middle    
                })

                inhole = True
                skip = True
                continue

            # Set inside hole flag to False when leaving the hole
            if breaks[i + 1]['group'] != b['group'] and inhole:
                inhole = False

        # Retrieve both associated projected points on the line
        p1 = shapely.Point(b['p1'])
        p2 = shapely.Point(b['p2'])

        # Split the original line to get the section between both projected point
        start, end = split_line_at_point(line, p1)
        middle, end = split_line_at_point(end, p2)

        maxdist = 0
        # Loop through the section's vertices
        for s in list(middle.coords):
            # Calculate the distance between the breakpoint and the vertex
            d = shapely.distance(shapely.Point(s), shapely.Point(b['geometry']))
            # If the distance is above the max distance, update it
            if d > maxdist:
                maxdist = d

        # If the distance is above 1.7 times the width of the symbol
        if maxdist > (1.7 * abs(width)):
            # Here, it is a strict pastiness conflict
            # Divide the line into three pieces start, middle, end
            sections.append({
                'type': 'pasty',
                'geometry': middle
            })

    if len(sections) > 0:
        sgdf = gpd.GeoDataFrame(sections, crs='EPSG:3857')
        sgdf.to_file("cartagen4py/data/sections.geojson", driver='GeoJSON')

    return groups