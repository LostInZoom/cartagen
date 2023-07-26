import networkx as nx
import numpy as np
import shapely, geopandas
from shapely.geometry.polygon import orient
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
    linear_ring = orient(face).exterior.coords
    for i, cc in enumerate(linear_ring):
        if cc in coordinates:
            node = coordinates.index(cc)
            if len(nodes[node]) < 3:
                cp = linear_ring[-2] if i == 0 else linear_ring[i - 1]
                cn = linear_ring[1] if i == len(linear_ring) - 1 else linear_ring[i + 1]
                alpha = int(angle_from_3points_coordinates(cp, cc, cn))
                nodes[node].append(alpha)

            if previous is not None:
                links.append([previous, node])
            previous = node

    return nodes, links