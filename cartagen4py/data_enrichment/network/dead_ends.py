import shapely
import geopandas as gpd
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

def detect_dead_ends(roads):
    """
    This function detects dead ends inside a road network and returns their geometries.
    """

    crs = roads.crs
    roads = roads.to_dict('records')

    network = []
    for road in roads:
        network.append(road['geometry'])

    faces = calculate_network_faces(network)

    # Create a storage for the dead ends road index
    deadends = []
    # Storage for future enclosed faces
    enclosed = []

    # Create a tree for the network roads and for the network faces
    netree = shapely.STRtree(network)
    facetree = shapely.STRtree(faces)
    
    # Loop through network faces
    for fid, face in enumerate(faces):
        # Retrieve roads inside the considered face
        # All those roads are part of dead ends groups
        contained = netree.query(face, predicate='contains').tolist()

        # Find holes inside the face
        holes = __find_holes(fid, faces, facetree)

        # Loop through each hole inside the considered face
        for hole in holes:
            # Add the roads intersecting the boundary of the hole to the list of dead ends roads
            contained.extend(netree.query(faces[hole].boundary, predicate='contains').tolist())

        # Retrieve the groups of roads
        connected, unconnected = __topological_grouping(contained, network, face)

        # Add the dead ends to the result
        for gid, group in enumerate(connected + unconnected):
            connect = True
            if gid >= len(connected):
                connect = False
            for cid in group:
                deadend = roads[cid]
                deadend['face'] = fid
                deadend['deid'] = gid
                deadend['connected'] = connect
                deadends.append(deadend)

    if len(deadends) > 0:
        return gpd.GeoDataFrame(deadends, crs=crs)
    else:
        return None

def __find_holes(index, faces, tree):
    """
    Returns a list of faces indexes if those faces are completely contained inside the given face index.
    Returns an empty list if none are found.
    """
    holes = []

    # Loop through faces intersecting the considered one
    for pid in tree.query(faces[index], predicate='intersects'):
        # Make sure the index is not the considered face
        if pid != index:
            count = 0
            # Loop through faces intersecting this one
            for pid1 in tree.query(faces[pid], predicate='intersects'):
                # Make sure it's not the same
                if pid1 != pid:
                    # Make sure the boundaries are not only crossing
                    if shapely.crosses(faces[pid].boundary, faces[pid1].boundary) == False:
                        count += 1
            
            # If count is 1, it only intersects one face, the considered
            # thus, it is a hole inside the face, append it to the result
            if count == 1:
                holes.append(pid)

    return holes

def __topological_grouping(indexes, geometries, face):
    """
    Create groups of roads depending on their topology, starting on the edge of the face.
    """

    # Recursively create groups
    def __recursive_grouping(current, group, indexes, leave, groups):
        # Storage for the following roads
        following = []
        # Loop through current roads
        for c in current:
            # Loop through all indexes
            for i in indexes:
                # If the indexes is not to be left
                if i not in leave:
                    # Check if the roads is following the current one
                    if shapely.intersects(geometries[c], geometries[i]):
                        # If so, add it as a following road, add it to be left alone, and add it to the current group
                        following.append(i)
                        leave.append(i)
                        group.append(i)

        # If roads are following the current ones
        if len(following) > 0:
            # Launch the function once again
            groups = __recursive_grouping(following, group, indexes, leave, groups)
        else:
            # If no roads are following, add the group to the list of groups
            # Here, the recursion is over
            groups.append(group)

    # Storage for the final groups of roads
    connected = []

    # Find the starting road(s) for each group
    entries = []
    for i in indexes:
        # Check if the geometry touches the exterior ring of the face
        if shapely.touches(face.exterior, geometries[i]):
            # Append it to the starting roads list
            entries.append(i)

    # A storage for already treated roads
    leave = []

    # Loop through each first road
    for e in entries:
        if e not in leave:
            group = [e]
            leave.append(e)

            # Retrieve an other touching starting road if exists
            # It can happen if the dead end starts with a junction or a hole
            for e1 in entries:
                if shapely.touches(geometries[e], geometries[e1]):
                    # Append the index to the list of entry roads to leave
                    group.append(e1)
                    leave.append(e1)
            
            # Launch the recursion to create groups
            __recursive_grouping(group, group, indexes, leave, connected)
    
    # Handle roads not connected to the network face boundary.
    # i.e. roads not connected to the rest of the network
    unconnected = []
    for i in indexes:
        # Here, avoid already treated roads
        if i not in leave:
            # Starting a new group
            group = [i]
            leave.append(i)

            # Launch the recursion to fill up the unconnected groups
            __recursive_grouping(group, group, indexes, leave, unconnected)

    return connected, unconnected
