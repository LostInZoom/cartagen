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
            for ext in crossroad.externals:
                original = None
                egeom = crossroad.network[ext]
                for u in unchanged:
                    ugeom = network[u]
                    if shapely.equals(egeom, ugeom):
                        original = roads[u]
                if original is not None:
                    incoming.append(original)
                else:
                    incoming.append({ "geometry": egeom })

            # Calculate the skeleton
            skeleton = SkeletonTIN(polygon, incoming=incoming, distance_douglas_peucker=3)
            # Storing the original geometries of the crossroad
            originals.extend(crossroad.original)
            # Storing the blended skeleton
            collapsed.extend(skeleton.blended)

    result = []
    remove = []
    for c in collapsed:
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