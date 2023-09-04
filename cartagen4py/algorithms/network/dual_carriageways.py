import shapely
import geopandas as gpd

from cartagen4py.utils.geometry import *
from cartagen4py.utils.network import *

def collapse_dual_carriageways(roads, carriageways):
    """
    Collapse dual carriageways using the polygon skeleton made from a Delaunay Triangulation
    """

    # Retrieve crs for output
    crs = roads.crs

    # Convert geodataframe to list of dicts
    roads = roads.to_dict('records')
    carriageways = carriageways.to_dict('records')

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

    # This will store future skeletons, polygons, and minimum rectangles
    skeletons = []
    polygons = []

    for carriageway in carriageways:
        # Get the geometry of the face
        polygon = carriageway['geometry']

        # Calculate the crossroad object
        crossroad = Crossroad(network, tree, polygon)

        # If there are more than one external road to the crossroad
        if len(crossroad.externals) > 1:
            # Retrieve the id of the external roads that have not been changed by the conversion to a crossroad object
            unchanged = crossroad.get_unchanged_roads('externals')

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
                    incoming.append({ "geometry": egeom })

            # Calculate the skeleton
            skeleton = SkeletonTIN(polygon, incoming=incoming, distance_douglas_peucker=3)

            skeletons.append(skeleton.network)
            polygons.append(polygon)



            # # Storing the original geometries of the crossroad
            # originals.extend(crossroad.original)
            # # Storing the blended skeleton
            # collapsed.extend(skeleton.blended)

    # result = []
    # remove = []
    # for c in collapsed:
    #     cgeom = c['geometry']
    #     add = True
    #     for o in originals:
    #         if shapely.equals(cgeom, network[o]):
    #             remove.append(o)
    #             add = False
    #     if add:
    #         result.append(c)

    # removeroad = []
    # for o in originals:
    #     if o not in remove:
    #         removeroad.append(o)

    # for rid, road in enumerate(roads):
    #     if rid not in removeroad:
    #         result.append(road)

    result = []
    for skeleton in skeletons:
        for bone in skeleton:
            result.append( {"geometry": bone} )

    # Stores carriageways that touches either on their long side or short side
    touching = [[], []]
    # Here, find touching dual carriageways
    # Loop through all carriageways polygons
    for pindex, p1 in enumerate(polygons):
        # Retrieve first boundary
        b1 = p1.boundary
        # Loop though all carriageways polygons
        for p2 in polygons:
            # Check if it's not the same polygon
            if p1 != p2:
                # If both boundaries overlap
                b2 = p2.boundary
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

                    # Here, the two polygon share an edge
                    # Calculate the polygon properties
                    face = NetworkFace(p1)
                    # Find whether the length or the width is closer to the length of the shared edge
                    if abs(line.length - face.length) < abs(line.length - face.width):
                        # Here the shared edge is the long side
                        add = True
                        for longside in touching[0]:
                            if shapely.equals(longside, line):
                                add = False
                        if add:
                            touching[0].append(line)
                    else:
                        # Here, the shared edge is the short side
                        add = True
                        for shortside in touching[1]:
                            if shapely.equals(shortside, line):
                                add = False
                        if add:
                            touching[1].append(line)

    tl = []
    for longside in touching[0]:
        tl.append({"geometry": longside})

    tl_gdf = gpd.GeoDataFrame(tl, crs=crs)
    tl_gdf.to_file("cartagen4py/data/touching_long.geojson", driver="GeoJSON")

    ts = []
    for shortside in touching[1]:
        ts.append({"geometry": shortside})

    ts_gdf = gpd.GeoDataFrame(ts, crs=crs)
    ts_gdf.to_file("cartagen4py/data/touching_short.geojson", driver="GeoJSON")

    return gpd.GeoDataFrame(result, crs=crs)

    # nodes = []
    # for i, n in enumerate(skeleton.nodes):
    #     nodes.append({
    #         'nid': i,
    #         'geometry': shapely.Point(n)
    #         })
    # n = gpd.GeoDataFrame(nodes, crs=crs)
    # n.to_file("cartagen4py/data/nodes.geojson", driver="GeoJSON")

    # edges = []
    # for j, s in enumerate(skeleton.edges):
    #     edges.append({
    #         'eid': j,
    #         'start': s[0],
    #         'end': s[1],
    #         'geometry': shapely.LineString([skeleton.nodes[s[0]], skeleton.nodes[s[1]]])
    #         })
    # e = gpd.GeoDataFrame(edges, crs=crs)
    # e.to_file("cartagen4py/data/edges.geojson", driver="GeoJSON")

    # joints = []
    # for j in skeleton.joints:
    #     joints.append({
    #         'geometry': j
    #     })
    # j = gpd.GeoDataFrame(joints, crs=crs)
    # j.to_file("cartagen4py/data/joints.geojson", driver="GeoJSON")

    # bones = []
    # for b in skeleton.bones:
    #     bones.append({
    #         'geometry': b
    #     })
    # b = gpd.GeoDataFrame(bones, crs=crs)
    # b.to_file("cartagen4py/data/bones.geojson", driver="GeoJSON")

    # blend = []
    # for i, b in enumerate(blended):
    #     blend.append({
    #         'nid': i,
    #         'geometry': b
    #         })
    # bl = gpd.GeoDataFrame(blend, crs=crs)
    # bl.to_file("cartagen4py/data/blended.geojson", driver="GeoJSON")