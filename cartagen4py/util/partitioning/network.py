import geopandas, shapely
from cartagen4py.util.extent import *
from shapely.ops import linemerge, unary_union, polygonize

# Calculates the faces of one or multiple networks and return a list of polygons
def network_partition(*networks):
    if len(networks) < 1:
        raise Exception('No networks provided, network partition cannot be created.')

    network = []
    for i, n in enumerate(networks):
        for object in n:
            if object.geom_type == 'LineString':
                network.append(object)
            elif object.geom_type == 'MultiLineString':
                for line in object.geoms:
                    network.append(line)

    # Calculating the networks convex hull
    multilines = unary_union(network)
    multilines = linemerge(network)
    hull = shapely.convex_hull(multilines)
    # Adding the convex hull boundary as a linestring to the network
    network.append(hull.boundary)
    
    # Fully dissolve and node the network
    network = unary_union(network)
    # Merge all contiguous lines
    network = linemerge(network)
    # Return a list of polygons representing the faces of the network
    return polygonize(network)