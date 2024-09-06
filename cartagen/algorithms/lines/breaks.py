import shapely, networkx

from cartagen.utils.geometry.dilation import dilate_line, offset_line, merge_connected_parts, reconstruct_line
from cartagen.utils.geometry.line import gaussian_smoothing, get_bend_side, merge_linestrings
from cartagen.utils.geometry.skeletonization import SkeletonTIN

def max_break(line, offset, exaggeration=1.0):
    """
    Spread a road bend and keep its shape.

    This algorithm was proposed by  Mustière. :footcite:p:`mustiere:2001`
    It first finds the side of the road bend and
    dilates the line on its exterior side.

    Parameters
    ----------
    line : LineString
        The line to dilate.
    offset : float
        The distance to dilate.
    exaggeration : float, optional
        A multiplicator for the offset parameter.

    Returns
    -------
    LineString

    See Also
    --------
    detect_pastiness :
        Detect the pastiness of a series of bends.
    min_break :
        Spread a road bend and minimize its extent.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (3, 3), (6, 0)])
    >>> max_break(line, 1.0)
    <LINESTRING (-0.707 0.707, 2.293 3.707, 2.426 3.819...)>
    """
    # Get the side of the bend
    side = get_bend_side(line)

    # Change the offset in case of left sided bend
    if side == 'left':
        offset = -offset

    # Dilate the bend
    dilated = dilate_line(line, offset*exaggeration, cap_style='flat', quad_segs=8)
    
    return dilated[0]

def min_break(line, offset, sigma=30, sample=None):
    """
    Spread a road bend and minimize its extent.

    This algorithm was proposed by  Mustière. :footcite:p:`mustiere:2001`
    It mimics the cartographic technique of overlapping the internal border of the bend.

    It creates a polygon from the linestring by closing it at it extremities,
    then it calculates the TIN skeleton, dilates it apply a gaussian smoothing.

    Parameters
    ----------
    line : LineString
        The line to dilate.
    offset : float
        The distance to dilate the skeleton.
    sigma : float, optional
        Gaussian smoothing strength.
    sample : int, optional
        Gaussian smoothing sample size.

    Returns
    -------
    LineString

    See Also
    --------
    detect_pastiness :
        Detect the pastiness of a series of bends.
    max_break :
        Spread a road bend and keep its shape.
    gaussian_smoothing :
        Smooth a line and attenuate its inflexion points.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (3, 3), (6, 0)])
    >>> min_break(line, 1.0)
    <LINESTRING (5 0, 5 1, 5 2, 4.97 2.347, 4.879 2.684...)>
    """

    # Create the offset of the skeleton
    def __offset(line, offset):
        # Calculate the offset points along the line
        groups = offset_line(line, offset)

        # Remove the circle projection from the first node
        groups[0]['projected'] = groups[0]['projected'][-1:]

        # Reconstruct the line into parts
        parts, breaks = reconstruct_line(groups, line, offset)

        # Merge parts that have a common set of coordinates
        groups = merge_connected_parts(parts)

        return groups

    # Get the nodes of the linestring
    coordinates = list(line.coords)
    start, end = coordinates[0], coordinates[-1]
    base = shapely.LineString([start, end])

    # Add the center of the line formed by the last and first node of the list
    entry = base.interpolate(base.length / 2)
    middle = entry.coords[0]

    coordinates.append(middle)

    # Create the polygon from the linestring nodes
    polygon = shapely.Polygon(coordinates)

    # Create a new skeleton object
    sk = SkeletonTIN(polygon)
    
    # Add a new entry which is the entroid of the line connecting the start and point of the bend
    sk.add_incoming_points([entry], connection='interior')

    # Remove start and end node of the bend from the skeleton
    sk.remove_entries([shapely.Point(start), shapely.Point(end)])

    # Create the network from the skeleton
    network = sk.create_network()

    # Check that a network exists
    if len(network) > 0:
        if len(network) == 1:
            skeleton = network[0]
        else:
            # Here the network is composed of more than one section
            # We nee to find the longest path
            # As there is more than one sections, it means the skeleton has interiors nodes
            # Interior nodes and entries represent the nodes of the graph
            nodes = [ i.coords[0] for i in sk.interiors ] + [ e.coords[0] for e in sk.entries ]
            edges = []

            for n in network:
                # Retrieve start and end point of the road
                startcoords, endcoords = n.coords[0], n.coords[-1]
                # Get the index of the corresponding node
                start, end = nodes.index(startcoords), nodes.index(endcoords)
                # Create the edge object
                edges.append((start, end, { 'length': n.length }))

            # Create the graph and add its edges
            graph = networkx.Graph()
            graph.add_edges_from(edges)

            # Get index of the root node, should be the last
            root = nodes.index(middle)

            # Storage for the longest path (length and actual node succession)
            longest = 0
            path = None
            # Loop through all the other nodes
            for target, tcoords in enumerate(nodes):
                if target != root:
                    # Calculate the shortest path between the root and the target node
                    dijkstra = networkx.shortest_path(graph, source=root, target=target, weight='length')
                    # Calculate the length of the path
                    length = networkx.path_weight(graph, dijkstra, 'length')
                    # If it's longer than the longest, overwrites it
                    if length > longest:
                        longest = length
                        path = dijkstra

            # Recreate links as pairs of nodes
            links = []
            if path is not None:
                for i in range(0, len(path) - 1):
                    # Append the tuple of coordinates of start and end nodes
                    links.append((nodes[path[i]], nodes[path[i + 1]]))

            skeleton = None
            # Loop through lines of the network
            for lid, l in enumerate(network):
                # Get start and end coords
                start, end = l.coords[0], l.coords[-1]
                # Check if the link is kept before merging it to the new skeleton
                if (start, end) in links or (end, start) in links:
                    if skeleton is None:
                        skeleton = l
                    else:
                        skeleton = merge_linestrings(skeleton, l)

        # Get the nodes of the skeleton without the last one
        coords = list(skeleton.coords)

        # Reverse the skeleton direction if the entry point is a the end of the skeleton
        if middle == coords[-1]:
            coords.reverse()
        else:
            # Raise an error if the entry point doesn't match the start of the skeleton
            if middle != coords[0]:
                raise Exception('There is problem with the generated skeleton.')

        # Remove the last node of the skeleton as it represents the other entry point
        coords.pop()

        # Reconstruc the line
        newline = shapely.LineString(coords)

        # Apply a gaussian smoothing to the line before offset
        newline = gaussian_smoothing(newline, sigma, sample)

        # Calculate the offset points along the skeleton
        left = list(__offset(newline, -offset)[0].coords)
        right = list(__offset(newline, offset)[0].coords)

        # Reverse the order on the right side
        right.reverse()

        # Return the line formed by the merging of left and right dilation
        return shapely.LineString(left + right)

    else:
        return line