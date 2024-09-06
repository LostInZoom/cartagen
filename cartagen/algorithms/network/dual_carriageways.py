import shapely, numpy
import geopandas as gpd

from cartagen.utils.attributes import attributes_from_longest
from cartagen.utils.network.roads import Crossroad
from cartagen.utils.network.faces import NetworkFace
from cartagen.utils.geometry.skeletonization import SkeletonTIN
from cartagen.utils.geometry.line import extend_line_with_point, merge_linestrings

def collapse_dual_carriageways(roads, carriageways, sigma=None, propagate_attributes=None):
    """
    Collapse dual carriageways using a TIN skeleton.

    This algorithm proposed by Thom :footcite:p:`thom:2005`
    collapses the network faces considered as dual carriageways
    using a skeleton calculated from a Delaunay triangulation.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        The road network.
    carriageways : GeoDataFrame of Polygon
        The polygons representing the faces of the network detected as dual carriageways.
    sigma : float, optional
        If not None, apply a gaussian smoothing to the collapsed dual carriageways to
        avoid jagged lines that can be created during the TIN skeleton creation.
    propagate_attributes : list of str, optional
        Propagate the provided list of column name to the resulting network.
        The propagated attribute is the one from the longest line.

    See Also
    --------
    detect_dual_carriageways : 
        Detect dual carriageways inside a road network.
    skeletonize_network :
        Blends a TIN skeleton inside a road network.

    References
    ----------
    .. footbibliography::
    """
    # Retrieve crs for output
    crs = roads.crs

    carriageways = carriageways.to_dict('records')

    if len(carriageways) == 0:
        return roads

    # Convert geodataframe to list of dicts
    roads = roads.to_dict('records')
    
    # Create a list of all the roads geometry of the network
    network = []
    for n in roads:
        network.append(n['geometry'])

    # Calculate the spatial index on roads
    tree = shapely.STRtree(network)

    # This list will store the indexes of roads to throw away
    originals = []
    # This will store the new collapsed geometries
    collapsed = []

    # Stores carriageways that touches either on their long side or short side
    longside = []
    shortside, shortside_list = [], []
    # Stores carriageways that touches on a single point
    pointside, pointside_list = [], []

    # Here, find touching dual carriageways
    # --------------------------------------------------------
    # Loop through all carriageways polygons
    for pindex1, p1 in enumerate(carriageways):
        geom1 = p1['geometry']
        # Retrieve first boundary
        b1 = geom1.boundary
        # Loop though all carriageways polygons
        for pindex2, p2 in enumerate(carriageways):
            geom2 = p2['geometry']
            # Check if it's not the same polygon
            if geom1 != geom2:
                # If both boundaries overlap
                b2 = geom2.boundary
                if b1.overlaps(b2):
                    # Calculate intersection between the two boundaries
                    i = shapely.intersection(b1, b2)
                    # Keep only lines forming the intersection
                    if i.geom_type == 'MultiLineString':
                        line = shapely.ops.linemerge(i.geoms)
                    elif i.geom_type == 'LineString':
                        line = i
                    else:
                        # If it's not a line, continue the loop
                        continue

                    # Here, the two polygons share an edge
                    # Calculate the polygon properties
                    face = NetworkFace(geom1)
                    # Find whether the length or the width is closer to the length of the shared edge
                    if abs(line.length - face.length) < abs(line.length - face.width):
                        # Here the shared edge is the long side
                        if pindex1 not in longside:
                            longside.append(pindex1)
                        if pindex2 not in longside:
                            longside.append(pindex2)
                    else:
                        # Here, the shared edge is the short side
                        add = True
                        for shorts in shortside:
                            if shapely.equals(shorts[2], line):
                                add = False
                        if add:
                            shortside.append([pindex1, pindex2, line])
                            if pindex1 not in shortside_list:
                                shortside_list.append(pindex1)
                            if pindex2 not in shortside_list:
                                shortside_list.append(pindex2)

                # Here, stores carriageways that touches at a single point
                elif b1.crosses(b2):
                    pointside.append([pindex1, pindex2])
                    if pindex1 not in pointside_list:
                        pointside_list.append(pindex1)
                    if pindex2 not in pointside_list:
                        pointside_list.append(pindex2)

    # This will store future skeletons
    skeletons = []

    # Here, calculate the tin skeleton if applicable
    # --------------------------------------------------------
    for cid, carriageway in enumerate(carriageways):
        # Get the geometry of the face
        polygon = carriageway['geometry']

        # Calculate the crossroad object
        crossroad = Crossroad(network, tree, polygon)

        # If there are more than one external road to the crossroad
        if len(crossroad.externals) > 1 and cid not in longside:
            # Retrieve the id of the external roads that have not been changed by the conversion to a crossroad object
            unchanged = crossroad.get_unchanged_roads('externals')

            # Check if attribute propagation is wanted
            attributes = None
            if propagate_attributes is not None:
                # Retrieve roads entry to calculate longest attributes
                proads = []
                boundary = polygon.boundary
                # Loop through the original roads id of the crossroad
                for rid in crossroad.original:
                    # Get the geometry of the road
                    rgeom = roads[rid]["geometry"]
                    # Check that the road is contained by the boundary of the carriageway or overlaps it
                    if boundary.contains(rgeom) or boundary.overlaps(rgeom):
                        # If so, append it
                        proads.append(roads[rid])
                
                # Retrieve the longest attributes value
                attributes = attributes_from_longest(proads, propagate_attributes)

            # Retrieve incoming roads, i.e. the external network of the crossroad
            incoming = []
            # Looping through external roads
            for ext in crossroad.externals:
                original = None
                # Retrieve the geometry of the line
                egeom = crossroad.network[ext]
                # Looping through unchanged network
                for u in unchanged:
                    # Retrieve geometry
                    ugeom = network[u]
                    # If the line equals an unchanged line
                    if shapely.equals(egeom, ugeom):
                        # Set original to be the road object with its attributes
                        original = roads[u]
                # If an unchanged line match has been found, add the object to the list
                if original is not None:
                    incoming.append(original)
                # Else, create a new object without attributes
                else:
                    # Check if attribute propagation is wanted
                    if attributes is not None:
                        incoming.append({ **attributes, **{"geometry": egeom} })
                    else:
                        incoming.append({ "geometry": egeom })

            # Calculate the skeleton
            skeleton = SkeletonTIN(polygon)
            skeleton.add_incoming_lines(incoming)
            skeleton.create_network()
            skeleton.blend(attributes, sigma=sigma)
            skeletons.append(skeleton)

            # Storing the original geometries of the crossroad
            originals.extend(crossroad.original)

        else:
            # Add None to keep the order of the indexes
            skeletons.append(None)

    # Here, carriageways connected by their short sides are treated
    # --------------------------------------------------------------

    # Get the bone junction between the two provided entry points
    def get_junction(points, skeleton):
        # Get the first junction starting from an entry point
        def get_next_joint(point, network):
            for ni, n in enumerate(network):
                # Retrieve start and end point of bone
                start, end = n.coords[0], n.coords[-1]
                if start == point:
                    return ni, end
                elif end == point:
                    return ni, start
            return None, None

        ni1, j1 = get_next_joint(points[0], skeleton.network)
        ni2, j2 = get_next_joint(points[1], skeleton.network)

        # If that junction is the same point, return the point
        if j1 == j2:
            return [ni1, ni2], j1
        # Else, return None
        else:
            return None, None

    # Remove 'fake' incoming lines from both skeletons object
    # i.e. incoming lines from one skeleton objects that is contained or overlap the other polygon geometry
    def remove_fake_incoming(skeleton1, skeleton2):
        def keep_loop(incoming, boundary):
            keep = []
            for i in incoming:
                if boundary.contains(i['geometry']) or boundary.overlaps(i['geometry']):
                    continue
                keep.append(i)
            return keep

        b1 = skeleton1.polygon.boundary
        b2 = skeleton2.polygon.boundary

        return keep_loop(skeleton1.incoming, b2), keep_loop(skeleton2.incoming, b1)

    # Retrieve incoming lines
    def retrieve_incoming_lines(incoming, entries, boundary):
        iresult = []
        for i, inc in enumerate(incoming):
            add = False
            # Add only if the line intersects an entry point
            for e in entries:
                if shapely.intersects(inc['geometry'], shapely.Point(e)):
                    add = True
            if add:
                iresult.append([i, inc])
        return iresult
    
    # Get the 'middle' line of a skeleton network
    def get_middle_line(junction, shortentries, skeleton):
        for ni, n in enumerate(skeleton.network):
            start, end = n.coords[0], n.coords[-1]
            if junction == start:
                if end not in shortentries:
                    return ni, n, 'start'
            elif junction == end:
                if start not in shortentries:
                    return ni, n, 'end'
        return None, None, None

    # Treating short side connections between carriageways
    for shorts in shortside:
        # Retrieve carriageways index and the shortside geometry
        cid1, cid2, shortline = shorts[0], shorts[1], shorts[2]

        if cid1 in longside or cid2 in longside:
            continue
        
        # Create a list containing entry points intersecting the short side
        shortentries = list(filter(lambda x: shapely.intersects(x, shortline), skeletons[cid1].entries))

        # Remove those entries from skeletons list of entries
        skeletons[cid1].entries = [x for x in skeletons[cid1].entries if x not in shortentries]
        skeletons[cid2].entries = [x for x in skeletons[cid2].entries if x not in shortentries]

        # Converts this list to tuples of coordinates
        shortentries = [x.coords[0] for x in shortentries]

        # Retrieve the junction between both entries inside both skeletons
        diamond1, j1 = get_junction(shortentries, skeletons[cid1])
        diamond2, j2 = get_junction(shortentries, skeletons[cid2])

        # Continue only if both junctions were found
        if j1 is not None and j2 is not None:
            # Remove 'fake' incoming lines
            skeletons[cid1].incoming, skeletons[cid2].incoming = remove_fake_incoming(skeletons[cid1], skeletons[cid2])

            # Retrieve incoming lines
            incoming = retrieve_incoming_lines(skeletons[cid1].incoming, shortentries, carriageways[cid2]['geometry'].boundary)

            # Get the middle line of both skeletons along with the direction of the line
            index1, line1, pos1 = get_middle_line(j1, shortentries, skeletons[cid1])
            index2, line2, pos2 = get_middle_line(j2, shortentries, skeletons[cid2])
            # Create the line between both junctions
            linejunction = shapely.LineString([j1, j2])

            # If there are incoming lines intersecting the short side entries
            if len(incoming) > 0:
                # Stores for distances and positions
                distances = []
                positions = []

                # Loop through incoming lines
                for i in incoming:
                    # Stores geometry
                    index = i[0]
                    geom = i[1]['geometry']

                    # Get start and end point
                    start, end = shapely.Point(geom.coords[0]), shapely.Point(geom.coords[-1])
                    # Calculate distance between start point and the line between junctions
                    startdist = shapely.distance(start, linejunction)
                    # Same for end point
                    enddist = shapely.distance(end, linejunction)

                    # Here, project the start or end point of the incmoming line on the junction line
                    # The stored distance represents the distance between the start of the junction line and the projected point on this same line
                    if startdist < enddist:
                        distances.append(linejunction.project(start))
                        positions.append('start')
                    else:
                        distances.append(linejunction.project(end))
                        positions.append('end')

                # Calculate the mean of the list of distances
                meandist = numpy.mean(distances)
                # This value is used to create the new intersection point on the junction line
                point = shapely.Point(linejunction.interpolate(meandist).coords)

                # Extend both interior lines with the new intersection point
                line1 = extend_line_with_point(line1, point, position=pos1)
                line2 = extend_line_with_point(line2, point, position=pos2)

                # Update each incoming line geometry with a new extended line with the intersection point
                for i, inc in enumerate(incoming):
                    inc[1]['geometry'] = extend_line_with_point(inc[1]['geometry'], point, position=positions[i])

                # Update the geometry of the 'middle' lines in their skeleton
                skeletons[cid1].network[index1] = line1
                skeletons[cid2].network[index2] = line2

            # Here, there is no incoming line, which means both interior lines can be merged along with the junction line
            else:
                # TODO: handle this situation by updating geometries somewhere...
                # This line is one long line which is both skeletons middle lines merged together.
                # This situation should not appear.
                merged = merge_linestrings(merge_linestrings(line1, linejunction), line2)

            # Remove the diamond shaped skeletons part
            skeletons[cid1].network = [x for i, x in enumerate(skeletons[cid1].network) if i not in diamond1]
            skeletons[cid2].network = [x for i, x in enumerate(skeletons[cid2].network) if i not in diamond2]
    
    shortdone, pointdone = [], []

    for skeleton in skeletons:
        if skeleton is not None:
            # Storing the blended skeleton
            collapsed.extend(skeleton.blended)

    result = []
    remove = []
    for i, c in enumerate(collapsed):
        cgeom = c['geometry']
        add = True
        for o in originals:
            if shapely.equals(cgeom, network[o]):
                remove.append(o)
                add = False
        if add:
            result.append(c)

    removeroad = []
    for o in originals:
        if o not in remove:
            removeroad.append(o)

    for rid, road in enumerate(roads):
        if rid not in removeroad:
            result.append(road)

    return gpd.GeoDataFrame(result, crs=crs)