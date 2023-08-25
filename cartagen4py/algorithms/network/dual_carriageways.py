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
    # This list will store the geometries of internal and external roads.
    # It is important when crossroads are too close to each other and external and internal roads might overlap.
    internals = []
    externals = []
    # This will store the new collapsed geometries
    collapsed = []

    for carriageway in carriageways:
        print(carriageway["cid"])
        polygon = carriageway['geometry']
        boundary = polygon.boundary
        lines = tree.query(polygon)

        # cnetwork = []
        # incoming = []
        # # Retrieve the entry points of the dual carriageway
        # for coords in boundary.coords:
        #     point = shapely.Point(coords)
        #     add = False
        #     for l in lines:
        #         line = network[l]
        #         if shapely.contains(boundary, line) == False:
        #             if shapely.intersects(line, point):
        #                 if l not in incoming:
        #                     cnetwork.append(line)
        #                     incoming.append(l)

        # Calculate the crossroad object
        crossroad = Crossroad(network, tree, polygon)

        # If there are more than one external road to the crossroad
        if len(crossroad.externals) > 1:
            incoming = []
            for ext in crossroad.externals:
                incoming.append(crossroad.network[ext])

            # Calculate the skeleton
            skeleton = SkeletonTIN(polygon, incoming=incoming, distance_douglas_peucker=3)
            i = []
            for internal in crossroad.internals:
                i.append(crossroad.network[internal])
            e = []
            for external in crossroad.externals:
                e.append(crossroad.network[external])
            internals.extend(i)
            externals.extend(e)
            originals.extend(crossroad.original)
            collapsed.extend(skeleton.blended)

    result = []
    for rid, road in enumerate(roads):
        if rid not in originals:
            result.append(road)

    for col in collapsed:
        result.append({
            'geometry': col
        })

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
            # break