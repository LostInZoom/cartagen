import shapely
import geopandas as gpd

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

    for carriageway in carriageways:
        polygon = carriageway['geometry']

        skeleton = SkeletonTIN(polygon)
    
        if carriageway['cid'] == 4:
            nodes = []
            for i, n in enumerate(skeleton.nodes):
                nodes.append({
                    'nid': i,
                    'geometry': shapely.Point(n)
                    })
            n = gpd.GeoDataFrame(nodes, crs=crs)
            n.to_file("cartagen4py/data/nodes.geojson", driver="GeoJSON")

            edges = []
            for j, s in enumerate(skeleton.edges):
                edges.append({
                    'eid': j,
                    'start': s[0],
                    'end': s[1],
                    'geometry': shapely.LineString([skeleton.nodes[s[0]], skeleton.nodes[s[1]]])
                    })
            e = gpd.GeoDataFrame(edges, crs=crs)
            e.to_file("cartagen4py/data/edges.geojson", driver="GeoJSON")

            # triangles = []
            # for s in skeleton.triangles:
            #     print(s)
            #     triangles.append({
            #         'e1': s[0],
            #         'e2': s[1],
            #         'e3': s[2],
            #         'geometry': shapely.Polygon([
            #             skeleton.nodes[skeleton.edges[s[0]][0]], skeleton.nodes[skeleton.edges[s[0]][1]],
            #             skeleton.nodes[skeleton.edges[s[1]][0]], skeleton.nodes[skeleton.edges[s[1]][1]],
            #             skeleton.nodes[skeleton.edges[s[2]][0]], skeleton.nodes[skeleton.edges[s[2]][1]]
            #         ])
            #     })
            # t = gpd.GeoDataFrame(triangles, crs=crs)
            # t.to_file("cartagen4py/data/triangles.geojson", driver="GeoJSON")
            break