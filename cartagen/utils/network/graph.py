import networkx as nx

def make_planar():
    """TODO: Make a functino to create a planar road network"""


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