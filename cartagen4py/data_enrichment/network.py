import shapely
import geopandas as gpd
import numpy as np
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *
from cartagen4py.utils.geometry.distances import *

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

def detect_roundabouts(roads, area_threshold=40000, miller_index=0.95):
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

def is_branching_crossroad(polygon, roads, area_threshold,
        maximum_distance_area, roundabouts=None,
        allow_middle_node=False, allow_single_4degree_node=False
    ):
    """
    Return True or False whether the given polygon is a branching crossroad or not depending on the given parameters.
    """

    def is_triangular(dist_area, polygon, nodes, middle=None):
        """
        Check if the network face is triangular given a maximum allowed distance area.
        """
        tnodes = []
        for i, node in enumerate(nodes):
            if middle != i:
                tnodes.append(node[1])

        tnodes.append(tnodes[0])
        triangle = shapely.Polygon(tnodes)

        da = distance_area(polygon, triangle)

        if da <= dist_area:
            return True
        
        return False

    if allow_middle_node:
        if allow_single_4degree_node:
            raise Exception('Parameters set incorrectly.')

    tree = shapely.STRtree(roads)

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

            # If there is not 3 or 4 nodes, this is not a branching crossroads
            if len(nodes) not in [3, 4]:
                return False
            
            # Case of a four nodes crossroads
            if len(nodes) == 4:
                # Check if middle node is allowed
                if allow_middle_node == False:
                    return False

                middle = None
                # Looping through each node
                for i, node in enumerate(nodes):
                    degree = node[0]
                    # If the degree is different than 3, return False
                    if degree != 3:
                        return False
                    # If the angle has been calculated
                    if 0 <= 2 < len(node):
                        angle = node[2]
                        # If the angle is 180, it is the middle node
                        if angle == 180:
                            middle = i
                    else:
                        # If an angle is missing, it means the polygon has a weird shape
                        return False
                # If there is no middle node, this is not a branching crossroad
                if middle is None:
                    return False
                else:
                    if is_triangular(maximum_distance_area, polygon, nodes, middle=middle) == False:
                        return False
                
                # Here, the 4 nodes crossroad is considered a branching crossroad
                return True

            # Case of a 3 nodes crossroad
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
                
                if roundabouts is not None:
                    rtree = shapely.STRtree(roundabouts)
                    # Retrieve roundabouts that intersects the considered network face using strtree
                    rintersects = rtree.query(polygon)
                    rb = []
                    for i in rintersects:
                        r = roundabouts[i]
                        # Make an other test to really keep only intersecting roundabouts, spatial index strtree using bbox
                        if shapely.intersects(polygon, r):
                            # Make sure, it is not the same object
                            if shapely.equals_exact(polygon, r):
                                continue
                            else:
                                rb.append(r)
                            
                    if len(rb) != 0:                        
                        # If more than one roundabout intersects, return False
                        if len(rb) > 1:
                            return False

                        # If the area of the face is twice larger than the area of the roundabout, return False
                        if area > (2 * rb[0].area):
                            return False

                        # There, this is a branching crossroad
                        return True

                if is_triangular(maximum_distance_area, polygon, nodes) == False:
                    return False

                # There, this is a branching crossroad
                return True

def detect_branching_crossroads(roads, area_threshold=2500,
        maximum_distance_area=0.5, roundabouts=None,
        allow_middle_node=False, allow_single_4degree_node=False
    ):
    """
    This function detects brandching crossroads inside a road network and returns polygons representing their extents.
    Parameters
    ----------
    roads : shapely geometry sequence of LineString.
        The road network to analyze.
    area_threshold : int optional.
        The area (in square meters) above which the object is not considered a branching crossroads.
        The default value is set to 2500.
    maximum_distance_area : float optional.
        The maximum distance area between the actual polygon and the triangle formed by the 3 nodes connecting the junction to the rest of the network.
        The default value is set to 0.5.
    roundabouts : shapely geometry sequence of Polygons optional.
        The polygons representing the network faces considered as roundabouts.
        If provided, it offers a better detection of branching crossroads.
        The default value is set to None.
    allow_middle_node : boolean optional
        If set to True, allow 4 nodes to form the crossroads, but each must have a degree of 3 and the 'middle' node must have an angle of 180°.
        Default value set to False.
    allow_single_4degree_node : boolean optional
        If set to True, allow one and only one node to have a degree of 4.
        Default value set to False.
    Cannot have both the allow_middle_node and allow_single_4degree_node parameters set to True.
    """

    faces = calculate_network_faces(roads, convex_hull=False)

    crossroads = []
    for i, face in enumerate(faces):
        if is_branching_crossroad(face, roads, area_threshold, maximum_distance_area, roundabouts, allow_middle_node, allow_single_4degree_node):
            crossroads.append(face)

    if len(crossroads) > 0:
        return crossroads
    else:
        return None