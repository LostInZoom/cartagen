import networkx as nx
import geopandas as gpd
import shapely
from shapely.ops import linemerge, unary_union
from shapely import STRtree, LineString

from cartagen.utils.geometry.conversion import multi_to_simple

def make_planar(network):
    """
    Make the network planar and keep the relationship between sections.
    
    This function fully merges the input network and reconstruct a planar network, i.e,
    a network having its edges intersect only at their endpoints.
    As an example, a road network will lose its tunnels and bridges as topological features
    and merging lanes will split the main road it merges into.

    This function is also suitable for river networks as it preserves the
    direction of the LineStrings.

    Parameters
    ----------
    network : GeoDataFrame of LineString
        The network to make planar.

    Returns
    -------
    GeoDataFrame of LineString

    Notes
    -----
    The resulting network will always have more network edges than provided as this function
    keeps the division of the original network even if it creates two-degree nodes. That way,
    the resulting edges will keep the attributes of the provided network.

    Examples
    --------
    >>> line1 = LineString([(0, 0), (0, 2)])
    >>> line2 = LineString([(1, 1), (0, 1)])
    >>> network = GeoDataFrame(geometry=[line1, line2])
    >>> make_planar(network)
        geometry                planar_id
    0   LINESTRING (0 0, 0 1)   0
    1   LINESTRING (1 1, 0 1)   1
    2   LINESTRING (0 1, 0 2)   2
    """
    def has_same_direction(line1, line2):
        """
        Returns True is the line1 has the same direction as the line2, else False.
        """
        direction = True
        # Get list of coords for each line
        coords1, coords2 = list(line1.coords), list(line2.coords)
        # Get start and end vertex of line2
        start, end = coords2[0], coords2[-1]

        # Boolean to check if the end has been found first
        endfirst = False
        # Loop through vertex of line1
        for vertex in coords1:
            # If the current vertex if the start of line2
            if vertex == start:
                # Break the loop and keep direction=True
                break
            # If the current vertex is the end of line2
            if vertex == end:
                # Set direction to False and break the loop
                direction = False
                break

        return direction

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
    pid = 0

    # Loop through resulting simple lines
    for resulting in lines:
        # Get the original lines intersecting the bbox of the current line
        close = tree.query(resulting)
        # Loop through intersecting lines
        for i in close:
            original = geometries[i]
            # If the resulting line is equal to the original one
            if shapely.equals(resulting, original):
                e = dict(records[i])
                e['planar_id'] = pid
                # Add the original entry to the result
                planar.append(e)
                pid += 1
            # If the resulting line contains the original
            elif shapely.contains(resulting, original):
                e = dict(records[i])
                e['planar_id'] = pid
                planar.append(e)
                pid += 1
            # If the resulting line is contained inside the original
            elif shapely.contains(original, resulting):
                e = dict(records[i])
                e['planar_id'] = pid
                geom = resulting if has_same_direction(original, resulting) else resulting.reverse()
                e['geometry'] = geom
                planar.append(e)
                pid += 1
            # If the resulting line and the original one overlaps
            elif shapely.overlaps(resulting, original):
                intersection = shapely.intersection(resulting, original)

                if intersection.geom_type == 'MultiLineString':
                    parts = [ i for i in intersection.geoms ]
                    union = unary_union(parts)
                    intersection = linemerge(union)

                # The intersection being contained in the original line
                # we can check its direction and reverse it if needed
                geom = intersection if has_same_direction(original, intersection) else intersection.reverse()
                e = dict(records[i])
                e['planar_id'] = pid
                e['geometry'] = geom
                planar.append(e)
                pid += 1

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