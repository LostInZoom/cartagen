import networkx as nx
import shapely, geopandas
from shapely.ops import linemerge, unary_union

def create_graph_from_linestrings(linestrings):
    """
    Create an undirected networkx graph from a set of linestrings or multilinestrings.
    Each node represents an intersection.
    """

    # Fully dissolve and node the network
    unioned = unary_union(linestrings)
    # Merge all contiguous lines
    merged = linemerge(unioned)

    network = []
    if merged.geom_type == 'LineString':
        network.append(merged)
    elif merged.geom_type == 'MultiLineString':
        for line in merged.geoms:
            network.append(line)

    print(network)



    lines = []
    for i, line in enumerate(network):
        lines.append({
            'id': i,
            'geometry': line
        })
    cgdf = geopandas.GeoDataFrame(lines, crs='EPSG:2154')
    cgdf.to_file("cartagen4py/data/test.geojson", driver="GeoJSON")