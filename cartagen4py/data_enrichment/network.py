import shapely as sh
import geopandas as gpd
import numpy as np
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

def is_roundabout(polygon, area_threshold, miller_index):
    """
    Return True or False whether the shape is a roundabout or not depending on the given parameters.
    """
    
    area = polygon.area
    # Check if the area of the polygon is larger than the threshold
    if area < area_threshold:
        # Calculate the Miller index for that polygon
        perimeter = polygon.length
        index = 4 * np.pi * area / (perimeter * perimeter)
        # Return True if this index is above the provided Miller Index
        if index > miller_index:
            return True
    
    return False

def detect_roundabouts(roads, area_threshold=40000, miller_index=0.97):
    """
    This function detects roundabouts inside a road network and returns polygons representing the roundabout extent.
    Returns None if no roundabouts where found.
    Parameters
    ----------
    roads : shapely geometry sequence of LineStrings
        The road network to analyze.
    area_threshold : int optional
        The area (in square meters) above which the object is not considered a roundabout.
        The default value is set to 40000.
    miller_index : float optional
        Index of compactess that determines if the shape is round or not.
        The default value is set to 0.97.
    """

    faces = calculate_network_faces(roads)

    roundabouts = []
    for face in faces:
        if is_roundabout(face, area_threshold, miller_index):
            roundabouts.append(face)

    if len(roundabouts) > 0:
        return roundabouts
    else:
        return None

def is_branching_crossroad(polygon, roads, area_threshold, maximum_surface_difference, allow_middle_node=False, allow_single_4degree_node=False):
    """
    Return True or False whether the given polygon is a branching crossroad or not depending on the given parameters.
    """

    if allow_middle_node:
        if allow_single_4degree_node:
            raise Exception('Parameters set incorrectly.')

    tree = sh.STRtree(roads)

    area = polygon.area
    # Check if the area of the polygon is larger than the threshold
    if area < area_threshold:
        # Retrieve objects that intersects the considered network face using strtree
        intersects = tree.query(polygon)

        lines = []
        for i in intersects:
            l = roads[i]
            # Make an other test to really keep only intersecting roads, spatial index strtree using bbox
            if shapely.intersects(polygon, l):
                lines.append(roads[i])

        # Continue only if more than two lines are selected, otherwise it can't be a branching crossroad
        if len(lines) > 2:
            # Fully dissolve and node the network
            unioned = unary_union(lines)
            # Merge all contiguous lines
            merged = linemerge(unioned)

            network = []
            if merged.geom_type == 'LineString':
                network.append(merged)
            elif merged.geom_type == 'MultiLineString':
                for line in merged.geoms:
                    network.append(line)

            nodes, links = create_graph_from_face(polygon, network)
            print(nodes)

            # If there is not 3 or 4 nodes, this is not a branching crossroads
            if len(nodes) not in [3, 4]:
                return False
            
            # Case of a four nodes crossroads
            if len(nodes) == 4:
                # Check if parameter is set to True
                if allow_middle_node:
                    for node in nodes:
                        degree = node[0]
                        # If a degree different than 3 is found, return False
                        if degree != 3:
                            return False
                else:
                    return False

            # Case of a three nodes crossroad
            if len(nodes) == 3:
                # Check if a 4 degree node is allowed
                if allow_single_4degree_node:
                    nb_four = 0
                    for node in nodes:
                        degree = node[0]
                        if degree not in [3, 4]:
                            return False
                        if degree == 4:
                            nb_four += 1
                    # If more than one node has a degree of 4, return False
                    if nb_four > 1:
                        return False

            return True

def detect_branching_crossroads(roads, area_threshold=2500, maximum_surface_difference=0.5, allow_middle_node=False, allow_single_4degree_node=False):
    """
    This function detects brandching crossroads inside a road network and returns polygons representing their extents.
    Parameters
    ----------
    roads : shapely geometry sequence of LineString.
        The road network to analyze.
    area_threshold : int optional.
        The area (in square meters) above which the object is not considered a branching crossroads.
        The default value is set to 2500.
    maximum_surface_difference : float optional.
        The maximum surface difference between the actual polygon and the triangle formed by the 3 nodes connecting the junction to the rest of the network.
        The default value is set to 0.5.
    allow_middle_node : boolean optional
        If set to True, allow 4 nodes to form the crossroads, but each must have a degree of 3 and the 'middle' node must have an angle of 180°.
        Default value set to False.
    allow_single_4degree_node : boolean optional
        If set to True, allow one and only one node to have a degree of 4.
        Default value set to False.
    Cannot have both the allow_3degree_4nodes and allow_single_4degree_node parameters set to True.
    """

    faces = calculate_network_faces(roads, convex_hull=False)

    crossroads = []
    for i, face in enumerate(faces):
        if is_branching_crossroad(face, roads, area_threshold, maximum_surface_difference, allow_middle_node, allow_single_4degree_node):
            crossroads.append(face)

    if len(crossroads) > 0:
        return crossroads
    else:
        return None