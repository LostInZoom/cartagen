import networkx as nx
import geopandas as gpd
import shapely
from shapely.ops import linemerge, unary_union
from shapely import STRtree, Point

from cartagen.utils.geometry.conversion import multi_to_simple
from cartagen.utils.debug import plot_debug, geojson_debug

def make_planar(network):
    """
    Make the network planar and keep the relationship between sections.
    
    This function fully merges the input network and reconstruct a planar network, i.e,
    a network having its edges intersect only at their endpoints.
    As an example, a road network will lose its tunnels and bridges as topological features
    and merging lanes will split the main road it merge into.

    The resulting network will always have more network edges as provided as this function
    keeps the division of the original network even if it creates two-degree nodes. That way,
    the resulting sections will keep the attributes of the provided network.

    Parameters
    ----------
    network : GeoDataFrame of LineString
        The network make planar.

    Returns
    -------
    GeoDataFrame of LineString
    """
    # Retrieve the crs
    crs = network.crs

    # Make sure the geomtry is simple
    simple = multi_to_simple(network)

    # Convert the network to a list of records
    records = simple.to_dict('records')

    # Retrieve a list of geometries
    geometries = [ e['geometry'] for e in records ]

    # Fully dissolve and node the network
    union = unary_union(geometries)
    # Merge all contiguous lines
    lines = linemerge(union)

    # Retrieve individual lines if the result is a multilinestring
    if lines.geom_type == 'MultiLineString':
        lines = [ l for l in lines.geoms ]
    
    # Creating the spatial index on the original lines
    tree = STRtree(geometries)

    planar = []

    # Loop through resulting simple lines
    for resulting in lines:
        # Get the original lines intersecting the bbox of the current line
        close = tree.query(resulting)
        # Loop through intersecting lines
        for i in close:
            original = geometries[i]
            # If the resulting line is equal to the original one
            if shapely.equals(resulting, original):
                # Add the original entry to the result
                planar.append(records[i])
            # If the resulting line contains the original
            elif shapely.contains(resulting, original):
                planar.append(records[i])
            # If the resulting line is contained inside the original
            elif shapely.contains(original, resulting):
                e = records[i].copy()
                e['geometry'] = resulting
                planar.append(e)
            # If the resulting line and the original one overlaps
            elif shapely.overlaps(resulting, original):
                e = records[i].copy()

                # This whole section is just a way of returning the overlapping geometry
                # because Shapely returns weird chunks of lines when using shapely.intersection() with two lines

                # Union and merge both lines into one
                union = unary_union([resulting, original])
                merge = linemerge(union)
                # Make sure it outputs a line
                if merge.geom_type == 'LineString':
                    # Storage for start and end of overlapping part
                    s1, s2 = None, None
                    
                    # The first and last vertex are the one that doesn't touch the merged line
                    # because touching nodes are the ones at extremities
                    start, end = Point(resulting.coords[0]), Point(resulting.coords[-1])
                    s1 = end if shapely.touches(merge, start) else start
                    start, end = Point(original.coords[0]), Point(original.coords[-1])
                    s2 = end if shapely.touches(merge, start) else start

                    # Split the line using the first vertex
                    split = list(shapely.ops.split(merge, s1).geoms)

                    # Set the final line as the one intersecting the last vertex and split it at the vertex
                    if split[0].intersects(s2):
                        final = list(shapely.ops.split(split[0], s2).geoms)
                    else:
                        final = list(shapely.ops.split(split[1], s2).geoms)

                    # Loop through splitted sections of line
                    for f in final:
                        # The right one is the one intersecting both vertexes
                        if f.intersects(s1) and f.intersects(s2):
                            e['geometry'] = f
                            planar.append(e)
    
    return gpd.GeoDataFrame(planar, crs=crs)


def longest_path(graph):
    """
    Get the longest path of the graph.
    """

def create_graph(roads, cost=None):
    """
    Create a graph object from a road network.

    Create a networkx graph object from the provided road network.
    The nodes of the returned graph has coordinates as attributes.
    The edges of the returned graph has cost and the index of the roads as attributes.
    
    Parameters
    ----------
    roads : list of dict from geopandas entries or geopandas GeoDataFrame of LineStrings.
        The road network to create the graph from.
    cost : str optional
        The name of the attribute giving the cost of the road section. Make sure the attribute is a number.
        Default to None, which means the length of the road is used as the cost.
    """
    roads = roads.to_dict('records')

    # Create storages for nodes and links of the graph
    nodes = []
    edges = []

    # Create the graph from the network
    for rid, road in enumerate(roads):
        # Retrieve geometry
        geom = road['geometry']

        # Retrieve start and end point of the road
        startcoords, endcoords = geom.coords[0], geom.coords[-1]
        start, end = None, None

        # If start coordinates already in nodes
        if startcoords in nodes:
            # Set the start as the index of the node
            start = nodes.index(startcoords)
        # Else...
        else:
            # Set the start as the length of the node list and add the new node
            start = len(nodes)
            nodes.append(startcoords)

        # Same as starting node for the end node
        if endcoords in nodes:
            end = nodes.index(endcoords)
        else:
            end = len(nodes)
            nodes.append(endcoords)

        # Set the weight value depending if cost is set or not
        weight = road[cost] if cost is not None else geom.length

        # Create the link with the start and end point along with attributes
        edges.append((start, end, { 'weight': weight, 'rid': rid }))

    # Create the graph and add its edges
    graph = nx.Graph()
    graph.add_edges_from(edges)

    # Add coordinates to the graph nodes
    for node in graph.nodes:
        graph.nodes[node]['coords'] = nodes[node]

    return graph