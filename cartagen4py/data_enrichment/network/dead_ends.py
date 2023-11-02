import shapely
import geopandas as gpd
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

def detect_dead_ends(roads):
    """
    This function detects dead ends inside a road network and returns their geometries.
    """

    crs = roads.crs

    network = []
    for road in roads.geometry:
        network.append(road)

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
        if fid == 84:
            # Retrieve roads inside the considered face
            # All those roads are part of dead ends groups
            contained = netree.query(face, predicate='contains').tolist()

            # Find holes inside the face
            holes = __find_holes(fid, faces, facetree)

            # Loop through each hole inside the considered face
            for hole in holes:
                # Add the roads intersecting the boundary of the hole to the list of dead ends roads
                contained.extend(netree.query(faces[hole].boundary, predicate='contains').tolist())

            groups = __topological_grouping(contained, network, holes, face)

            for cid in groups:
                deadends.append({ "geometry": network[cid] })

            # Creating a group for the dead ends
            group = []

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

def __topological_grouping(indexes, geometries, holes, face):
    """
    Create groups of roads depending on their topology, starting on the edge of the face.
    """

    def __recursive_grouping():
        pass

    firsts = []
    # Find the starting road(s) for each group
    for i in indexes:
        # Check if the geometry touches the exterior ring of the face
        if shapely.touches(face.exterior, geometries[i]):
            # Append it to the starting roads list
            firsts.append(i)

    # A storage for already treated starting roads
    leave = []
    # A storage for already treated road
    done = []

    # Loop through each first road
    for f in firsts:
        # Check that the road has not already been treated
        if f not in leave:
            # Retrieve an other touching starting road if exists
            # It can happen if the dead end starts with a junction.
            for f1 in firsts:
                if shapely.touches(geometries[f], geometries[f1]):
                    # Append the index to the list of road to leave
                    leave.append(f1)

            current = f

            while pursue:
                group = []

                remove = None

                for i in indexes:
                    if shapely.intersects(geometries[current], geometries[i]):
                        group.append(i)
                        remove = i

                if remove is not None:
                    i.pop(remove)
    
    return firsts
