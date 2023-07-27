import shapely
import geopandas as gpd
import numpy as np
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

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

def detect_roundabouts(roads, area_threshold=40000, miller_index=0.95):
    """
    This function detects roundabouts inside a road network and returns polygons representing the roundabout extent.
    Returns None if no roundabouts where found.
    Parameters
    ----------
    roads : geopandas GeoDataFrame of LineStrings
        The road network to analyze.
    area_threshold : int optional
        The area (in square meters) above which the object is not considered a roundabout.
        The default value is set to 40000.
    miller_index : float optional
        Index of compactess that determines if the shape is round or not.
        The default value is set to 0.97.
    """

    network = []
    for road in roads.geometry:
        network.append(road)

    faces = calculate_network_faces(network)

    roundabouts = []
    for i, face in enumerate(faces):
        add, infos = is_roundabout(face, area_threshold, miller_index)
        if add:
            infos['cid'] = i
            roundabouts.append(infos)

    if len(roundabouts) > 0:
        return gpd.GeoDataFrame(roundabouts, crs=roads.crs)
    else:
        return None