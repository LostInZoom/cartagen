import shapely
from shapely.ops import linemerge, unary_union, polygonize
from cartagen4py.util.extent import *

# Calculates the faces of one or multiple networks and return a list of polygons
def calculate_network_faces(*networks):
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

# Returns a tuple with objects and faces polygons
# TODO: Manage centroids that happens to be on the edges of the networks faces
def network_partition(objects, *networks):
    # Create an empty tuple to store future partitions
    partition = ([], [])

    # Calculate the network faces from the networks provided
    faces = calculate_network_faces(*networks)

    # Calculate the centroids of each objects and store them in a list
    centroids = []
    for obj in objects:
        centroid = shapely.centroid(obj)
        centroids.append(centroid)

    # Create the spatial index for the objects centroids
    tree = shapely.STRtree(centroids)

    # Loop through all network faces
    for face in faces:
        # Retrieve objects that intersects the considered network face
        intersects = tree.query(face)
        # If objects are intersecting the considered network face
        if len(intersects) > 0:
            group = []
            # Loop through those objects
            for i in intersects:
                group.append(objects[i])
            partition[0].append(group)
            partition[1].append(face)
    return partition