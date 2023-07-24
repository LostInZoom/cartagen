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

def is_branching_crossroad(polygon, roads, area_threshold, maximum_surface_difference):
    """
    Return True or False whether the given polygon is a branching crossroad or not depending on the given parameters.
    """

    tree = sh.STRtree(roads)

    area = face.area
    # Check if the area of the polygon is larger than the threshold
    if area < area_threshold:
        # Retrieve objects that intersects the considered network face using strtree
        intersects = tree.query(face)

        lines = []
        for i in intersects:
            l = roads[i]
            # Make an other test to really keep only intersecting roads, spatial index strtree using bbox
            if shapely.intersects(face, l):
                lines.append(roads[i])

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

        

def detect_branching_crossroads(roads, area_threshold=2500, maximum_surface_difference=0.5):
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
    """

    faces = calculate_network_faces(roads)

    crossroads = []
    for face in faces:
        if is_branching_crossroad(face, roads, area_threshold, maximum_surface_difference):
            roundabouts.append(face)

    if len(crossroads) > 0:
        return crossroads
    else:
        return None