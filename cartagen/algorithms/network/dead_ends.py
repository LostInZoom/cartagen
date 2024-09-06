import itertools

import networkx
import geopandas as gpd

def eliminate_dead_ends(roads, length=0, keep_longest=True):
    """
    Eliminate dead-ends groups.

    Eliminates roads section detected as dead-ends inside a road network
    if the length of their main component is below a given threshold.
        
    Parameters
    ----------
    roads : GeoDataFrame of LineString
        The road network with dead-ends properly characterized.
    length : float, optional
        Length below which dead-ends are eliminated.
        If left to 0, shorter dead-ends are not eliminated.
    keep_longest : bool, optional
        If set to True, in case of complex dead-end, keep only the main component (*c.f.* description).
        If False, eliminate the whole dead-end group.
        This concerns also dead-ends with a length above the provided threshold.

    Returns
    -------
    GeoDataFrame of LineString

    Warning
    -------
    This function requires the detection of dead-ends beforehand.

    See Also
    --------
    detect_dead_ends :
        Characterize dead-ends groups inside a road network.

    Notes
    -----
    The **main component** of a dead-end is defined as:
    
    - If the dead-end is simple (*i.e.* just one road), the main component is the road.
    - If the dead-end contains multiple ramification of roads,
      the main component represents the path between the entry and the longest ramification.
    - If the dead-end contains inner network faces (*i.e.* enclosed roads),
      the main component represents the longest of the shortest paths
      between the entry and all the vertex of the dead ends.
    """

    # Retrieve crs for output
    crs = roads.crs

    # Convert geodataframe to list of dicts
    roads = roads.to_dict('records')

    deadends = []
    # Reconstruct a list with only dead ends
    for r in roads:
        if r['deadend']:
            deadends.append(r)

    if len(deadends) == 0:
        return gpd.GeoDataFrame(roads, crs=crs)

    # Storage for indexes to remove from the results
    remove = []

    # Create filters
    facefilter = lambda d: d['face']
    deadendfilter = lambda d: d['deid']
    # Loop through faces
    deadendgroups = []
    deadengeom = []
    for fid, facegroup in itertools.groupby(sorted(deadends, key=facefilter), key=facefilter):
        group = []
        geomgroup = []
        # Loop through dead end groups
        for deid, deadend in itertools.groupby(sorted(facegroup, key=deadendfilter), key=deadendfilter):
            # if True: # deid == 12:
            # Get the list of entries
            deadendlist = list(deadend)
            group.append(deadendlist)
            # Retrieve individual geometries
            geoms = []
            for d in deadendlist:
                geoms.append(d['geometry'])
            geomgroup.append(geoms)
        # Add entries to both lists
        deadendgroups.append(group)
        deadengeom.append(geomgroup)

    # Now, deadendgroups is a list of a list of roads with all their initial attributes
    # and deadendgeom is the same but with only geometries
    # Loop through each faces
    for fid, face in enumerate(deadengeom):
        # Loop through each group inside this face
        for gid, group in enumerate(face):
            # Here, treating unconnected dead ends
            if deadendgroups[fid][gid][0]['connected'] is False:
                # Add them to the remove list
                remove = __add_to_remove(deadendgroups[fid][gid], remove)
            else:
                # Check that the group have more than one road 
                if len(group) > 1:
                    # Create storages for nodes and links of the graph
                    nodes = []
                    edges = []

                    # Loop through roads inside the group
                    for rid, road in enumerate(group):
                        # Retrieve start and end point of the road
                        startcoords, endcoords = road.coords[0], road.coords[-1]
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

                        # Create the link with the start and end point along with the length of the road
                        edges.append((start, end, { 'length': road.length, 'root': deadendgroups[fid][gid][rid]['root'] }))

                    # Create the graph and add its edges
                    graph = networkx.Graph()
                    graph.add_edges_from(edges)

                    # Here, we want to find the root node
                    # Loop through the graph nodes
                    root = None
                    current = 0
                    for node in graph.nodes:
                        # Get the degree of the node
                        degree = graph.degree(node)
                        # Get the number of connected edges that are root roads
                        redges = len([x for x in graph.edges(node, 'root') if x[2]])

                        # If the degree is equal to the number of connected root roads
                        # It is the root node
                        if degree == redges:
                            # Update only if the degree is higher than the previous one
                            if degree > current:
                                root = node
                                # Update the current one
                                current = degree

                    # Check if root is none
                    if root is None:
                        # If so, add it to the remove list
                        remove = __add_to_remove(deadendgroups[fid][gid], remove)
                    else:
                        # Check if a hole is in the dead end group
                        if __has_hole(deadendgroups[fid][gid]):
                            targets = graph.nodes
                        else:
                            targets = [ n for n in graph.nodes if graph.degree(n) == 1 ]

                        longest = 0
                        path = None
                        # Loop through target nodes
                        for target in targets:
                            if networkx.has_path(graph, source=root, target=target):
                                # Calculate the shortest path between the root and the target node
                                dijkstra = networkx.shortest_path(graph, source=root, target=target, weight='length')
                                # Calculate the length of the path
                                plength = networkx.path_weight(graph, dijkstra, 'length')
                                # If it's longer than the longest, overwrites it
                                if plength > longest:
                                    longest = plength
                                    path = dijkstra

                        # If the longest is below the given length threshold, remove all roads
                        if longest < length:
                            # Add it to the remove list
                            remove = __add_to_remove(deadendgroups[fid][gid], remove)
                        else:
                            # Check if the longest path only should be kept
                            if keep_longest:
                                # Storage for links
                                links = []
                                if path is not None:
                                    # Loop through the path nodes
                                    for inode, node in enumerate(path):
                                        if inode < len(path) - 1:
                                            # Append the tuple of coordinates of start and end nodes
                                            links.append((nodes[node], nodes[path[inode + 1]]))
                                    # Loop through roads of the group
                                    for rid, road in enumerate(group):
                                        # Get start and end coords
                                        start, end = road.coords[0], road.coords[-1]
                                        # Check if the link is to be kept
                                        if (start, end) in links:
                                            continue
                                        elif (end, start) in links:
                                            continue
                                        else:
                                            # Add to remove if not to be kept
                                            remove.append(deadendgroups[fid][gid][rid]['rid'])

                # If the group has only one road
                else:
                    # If the its length is below the threshold
                    if group[0].length < length:
                        # Add it to the remove list
                        remove = __add_to_remove(deadendgroups[fid][gid], remove)

    result = []
    # Loop through roads
    for r in roads:
        # Check if the road is to be kept
        if r['rid'] not in remove:
            result.append(r)

    return gpd.GeoDataFrame(result, crs=crs)

# Return true if the group has a hole
def __has_hole(group):
    hashole = False
    # Loop through the group
    for road in group:
        # If it's a hole road
        if road['hole']:
            # Set to true and break
            hashole = True
            break
    # Return if a hole has been found
    return hashole

# Add the group for removal
def __add_to_remove(group, remove):
    # Loop through the group roads
    for road in group:
        # Append the rid to remove
        remove.append(road['rid'])
    # Return the list to remove
    return remove