import shapely
from shapely.ops import linemerge, unary_union, polygonize

def network_faces(*networks, convex_hull=True):
    """
    Calculates the faces of one or multiple networks.

    This function is often called polygonize.
    It creates polygons in place of a network of lines.

    Parameters
    ----------
    *networks : list of LineString
        The networks to generate the faces from.
        All the provided lists are merged into one list of geometries.
    convex_hull : bool, optional
        If True, add the convex hull of the provided list of lines to the
        network, thus including the borders.

    Returns
    -------
    list of Polygon

    See Also
    --------
    partition_networks :
        Partition objects using one or multiple networks.

    Notes
    -----
    When the network is not completely planar (*e.g.* tunnels
    and bridges inside a road network with crossing
    sections not intersecting), the network
    is fully unioned and merged.

    Examples
    --------
    >>> lines = [LineString([(0, 0), (1, 1), (2, 0), (3, 1)])]
    >>> list(network_faces(lines))
    [<POLYGON ((1 1, 3 1, 2 0, 1 1))>, <POLYGON ((2 0, 0 0, 1 1, 2 0))>]
    """
    if len(networks) == 0:
        raise Exception('No networks provided, network partition cannot be created.')

    network = []
    for i, n in enumerate(networks):
        for object in n:
            if object.geom_type == 'LineString':
                network.append(object)
            elif object.geom_type == 'MultiLineString':
                for line in object.geoms:
                    network.append(line)

    if convex_hull:
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
def partition_networks(objects, *networks):
    """
    Partition objects using one or multiple networks.

    Create the faces from one or multiple networks and
    partition the provided objects into groups. Objects are
    grouped when their centroid intersects the same
    network face. This algorithm rely on a
    :class:`STRtree <shapely.STRtree>` indexing solution. 

    Parameters
    ----------
    objects : GeoDataFrame of Geometry
        The objects to partition. Their centroids are used
        to assign them to a network face.
    *networks : GeoDataFrame of LineString
        The networks used to partition the objects.

    Returns
    -------
    partition : tuple
        A tuple containing two elements :

        #. A list of lists of index ordered by the network faces
        #. A list of the geometry of the network faces

    See Also
    --------
    network_faces :
        Calculates the faces of one or multiple networks.
    partition_grid :
        Partition objects using a grid of a given shape.

    Examples
    --------
    >>> points = GeoDataFrame(geometry=[ Point(2, 1), Point(4, 1) ])
    >>> network = GeoDataFrame(geometry=[ LineString([(0, 0), (3, 0), (3, 2), (6, 2)]) ])
    >>> partition_networks(points, network)
    ([[1], [0]], [<POLYGON ((3 2, 6 2, 3 0, 3 2))>, <POLYGON ((3 2, 3 0, 0 0, 3 2))>])
    """
    obj = objects.to_dict('records')

    # Create an empty tuple to store future partitions
    partition = ([], [])

    shapes = []
    for network in networks:
        shapes.append(network.geometry)

    # Calculate the network faces from the networks provided
    faces = network_faces(*shapes)

    # Calculate the centroids of each objects and store them in a list
    centroids = []
    for o in obj:
        if o['geometry'].geom_type == 'Polygon':
            centroids.append(o['geometry'].point_on_surface())
        else:
            centroids.append(o['geometry'].centroid)

    # Create the spatial index for the objects centroids
    tree = shapely.STRtree(centroids)

    # Loop through all network faces
    for face in faces:
        # Retrieve objects that intersects the considered network face
        intersects = list(tree.query(face, predicate='intersects'))

        indexes = []
        # Make sure it truly intersects before adding
        for i in intersects:
            if shapely.intersects(centroids[i], face):
                indexes.append(int(i))

        # If objects are intersecting the considered network face
        if len(indexes) > 0:
            partition[0].append(indexes)
            partition[1].append(face)

    return partition