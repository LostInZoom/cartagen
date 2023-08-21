import shapely
import geopandas as gpd
import numpy as np
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

def detect_roundabouts(network, area_threshold=40000, miller_index=0.95):
    """
    This function detects roundabouts inside a road network and returns polygons representing the roundabout extent.
    Returns None if no roundabouts where found.
    Parameters
    ----------
    network : geopandas GeoDataFrame of LineStrings
        The road network to analyze.
    area_threshold : int optional
        The area (in square meters) above which the object is not considered a roundabout.
        The default value is set to 40000.
    miller_index : float optional
        Index of compactess that determines if the shape is round or not.
        The default value is set to 0.97.
    """
    crs = network.crs

    roads = []
    for road in network.geometry:
        roads.append(road)

    faces = calculate_network_faces(roads, convex_hull=False)

    roundabouts = []
    index = 0
    for face in faces:
        add, infos = is_roundabout(face, area_threshold, miller_index)
        if add:
            infos['cid'] = index
            roundabouts.append(infos)
            index += 1

    if len(roundabouts) > 0:
        return gpd.GeoDataFrame(roundabouts, crs=crs)
    else:
        return None

def is_roundabout(polygon, area_threshold, miller_index):
    """
    Return True or False whether the shape is a roundabout or not depending on the given parameters.
    """

    infos = { 'geometry': polygon }
    
    area = polygon.area
    # Check if the area of the polygon is larger than the threshold
    if area < area_threshold:
        # Calculate the Miller index for that polygon
        perimeter = polygon.length
        index = 4 * np.pi * area / (perimeter * perimeter)
        # Return True if this index is above the provided Miller Index
        if index > miller_index:
            infos['index'] = index
            return True, infos
    
    return False, infos