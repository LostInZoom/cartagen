import networkx as nx
import numpy as np
import shapely, geopandas
from shapely.ops import linemerge, unary_union
from cartagen4py.utils.geometry import *

def find_line(start, end, linestrings):
    """
    Find a linestring inside a sequence of linestrings given a start and end point as coordinates
    """
    line = None
    for linestring in linestrings:
        l = linestring.coords
        lstart, lend = l[0], l[-1]
        if (start == lstart) and (end == lend):
            line = linestring
        elif (start == lend) and (end == lstart):
            line = linestring
    return line

def create_graph_from_face(face, network):
    """
    Return a list of nodes with their degree and a list a link between the index of the nodes.
    """
    coordinates = []
    nodes = []
    index = 0

    def update_nodes(node, index, coordinates):
        if node in coordinates:
            i = coordinates.index(node)
            nodes[i][0] += 1
        else:
            if shapely.intersects(face, shapely.Point(node)):
                nodes.append([1, node])
                coordinates.append(node)
                index += 1

    for link in network:
        linenodes = link.coords
        firstnode = linenodes[0]
        lastnode = linenodes[-1]
        update_nodes(firstnode, index, coordinates)
        update_nodes(lastnode, index, coordinates)

    links = []
    previous = None
    linear_ring = face.exterior.coords
    for i, cc in enumerate(linear_ring):
        if cc in coordinates:
            cp = linear_ring[-1] if i == 0 else linear_ring[i - 1]
            cn = linear_ring[0] if i == len(linear_ring) - 1 else linear_ring[i + 1]
            xcp, ycp = cp[0], cp[1]
            xcc, ycc = cc[0], cc[1]
            xcn, ycn = cn[0], cn[1]
            alpha = np.arctan2(ycn - ycc, xcn - xcc) - np.arctan2(xcp - xcc, ycp - ycc)
            
            node = coordinates.index(cc)
            nodes[node].append(alpha)
            if previous is not None:
                links.append([previous, node])
            previous = node

    if len(nodes) > 5:
        nodeslist = []
        for i, node in enumerate(coordinates):
            nodeslist.append({
                'id': i,
                'degree': nodes[i][0],
                'angle': nodes[i][2],
                'geometry': shapely.Point(node)
            })
        ngdf = geopandas.GeoDataFrame(nodeslist, crs='EPSG:2154')
        ngdf.to_file("cartagen4py/data/network_nodes.geojson", driver="GeoJSON")

        linkslist = []
        for i, link in enumerate(links):
            linkslist.append({
                'id': i,
                'origin': link[0],
                'destination': link[1],
                'geometry': shapely.LineString([coordinates[link[0]], coordinates[link[1]]])
            })
        lgdf = geopandas.GeoDataFrame(linkslist, crs='EPSG:2154')
        lgdf.to_file("cartagen4py/data/network_links.geojson", driver="GeoJSON")
    
    return nodes, links