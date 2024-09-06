import shapely
import geopandas as gpd
from cartagen.utils.partitioning.network import network_faces
from cartagen.utils.network import *

def detect_dead_ends(roads, outside_faces=False):
    """
    Detect dead-ends groups.

    This function detects dead-ends groups inside a road network.
    A dead-end group is a connected group of road sections that
    overlaps the border of only one network face.
    A connected group of road sections not connected to the rest
    of the network are considered dead-ends.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        The road network.
    outside_faces : bool, optional
        Whether dead-ends should be calculated on the outside faces
        of the road network. This can induce wrong characterization
        on the border of the provided dataset.

    Returns
    -------
    roads : GeoDataFrame of LineString
        The input road network with new attributes:

        - *'deadend'*: boolean indicating whether the road is part of a dead-end group. 
        - *'face'*: Index of the network face it belongs to.
        - *'deid'*: Index of the dead end group inside a given face.
        - *'connected'*: Set to true if the dead end group is connected to the network.
        - *'root'*: Set to true if the road section is the root of the dead end group,
          *i.e.* the section connecting the dead end group to the road network.
        - *'hole'*: Set to true if the road section touches a hole inside the dead end group.
    
    See Also
    --------
    eliminate_dead_ends
    """

    crs = roads.crs
    roads = roads.to_dict('records')

    network = []
    for road in roads:
        network.append(road['geometry'])

    faces = network_faces(network, convex_hull=outside_faces)

    hull = None
    if outside_faces:
        hull = shapely.convex_hull(shapely.MultiPolygon(faces)).boundary

    # Create a storage for the dead ends road index
    deadends = []

    # Create a tree for the network roads and for the network faces
    netree = shapely.STRtree(network)
    facetree = shapely.STRtree(faces)

    # This list will store indexes of faces that are holes
    # Those faces won't be treated individually, but as part of the face they are inside of
    leave = []
    # This list list will store the list of faces contained inside each face at the right index
    # It will contain an empty list where no faces are contained
    iholes = []

    # Loop through network faces
    for fid, face in enumerate(faces):
        # Retrieve the list of holes inside the face
        holes = list(faces[fid].interiors)

        # If there are holes inside the face
        if len(holes) > 0:
            ihole = []
            # Loop through individual holes
            for hole in holes:
                # Retrieve the face(s) contained inside the hole
                faceholes = facetree.query(shapely.Polygon(hole), predicate='contains').tolist()
                # Add the faces to the storage of faces not to treat
                leave.extend(faceholes)
                # Extend the list of faces contained inside this face
                ihole.extend(faceholes)
            iholes.append(ihole)
        else:
            iholes.append([])
                
    deadindex = []

    # Loop through network faces
    for fid, face in enumerate(faces):
        # Make sure the face is not inside the hole of an other face
        if fid not in leave:
            # Retrieve roads inside the considered face
            # All those roads are part of dead ends groups
            contained = netree.query(face, predicate='contains').tolist()

            # Storage for the roads covering the holes of the face
            holeroads = []

            # Loop through each faces inside this hole
            for hole in iholes[fid]:
                # Retrieve the roads inside and on the borders of the hole
                covering = netree.query(faces[hole], predicate='covers').tolist()
                # Add them to the list of hole roads
                holeroads.extend(covering)

            # Add all hole roads to the dead end groups
            contained.extend(holeroads)
        
            # Retrieve the groups of roads
            connected, nb_entries, unconnected = __topological_grouping(contained, network, face, hull)

            # Add the dead ends to the result
            # Loop through all dead end groups
            for gid, group in enumerate(connected + unconnected):
                connect, root = True, False
                # Loop through each road of the group
                for gcid, cid in enumerate(group):
                    # If the loop goes through unconnected groups, set to unconnected
                    if gid >= len(connected):
                        connect = False
                    else:
                        # Set the road as the entrance of the group if its id is below the entrance number of the group
                        root = True if gcid < nb_entries[gid] else False
                    # Set the road as being a hole boundary if it is
                    ishole = True if cid in holeroads else False
                    # Create the entry
                    deadend = roads[cid]
                    deadend['deadend'] = True
                    deadend['rid'] = cid
                    deadend['face'] = fid
                    deadend['deid'] = gid
                    deadend['connected'] = connect
                    deadend['root'] = root
                    deadend['hole'] = ishole
                    # Append the entry
                    deadends.append(deadend)
                    deadindex.append(cid)

    # Add all roads that are not dead-ends back to the result
    for cid, r in enumerate(roads):
        if cid not in deadindex:
            r['deadend'] = False
            deadends.append(r)

    return gpd.GeoDataFrame(deadends, crs=crs)


def __topological_grouping(indexes, geometries, face, hull):
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
    # Storage for the number of entries in the dead end group
    nb_entries = []

    # Find the starting road(s) for each group
    entries = []
    for i in indexes:
        # Check if the geometry touches the exterior ring of the face
        if shapely.touches(face.exterior, geometries[i]):
            # If outside faces has been included
            if hull is not None:
                # If the geometry touches the boundary of the convex hull, do not add it
                if shapely.touches(hull, geometries[i]):
                    continue
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
            
            # Set the number of entries for this group
            nb_entries.append(len(group))
            
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

    return connected, nb_entries, unconnected
