import shapely
import geopandas as gpd
import numpy as np
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

def detect_roundabouts(roads, area_threshold=40000, miller_index=0.95):
    """
    Detect roundabouts based on geometric properties (Touya, 2010).

    Detects roundabouts inside a road network.

    Parameters
    ----------
    network : GeoPandas.GeoDataFrame with LineString geometries
        Road network to analyze.
    area_threshold : int, Default=40000
        The area (in square meters) above which the object is not considered a roundabout.
    miller_index : float, Default=0.97
        Index of compactess that determines if the shape is round or not.

    Returns
    -------
    GeoPandas.GeoDataFrame of Polygon geometries

    See Also
    --------
    detect_branching_crossroads
    collapse_roundabouts
    collapse_branching_crossroads
    """
    crs = roads.crs
    roads = roads.to_dict('records')

    network = []
    for road in roads:
        network.append(road['geometry'])

    faces = calculate_network_faces(network, convex_hull=False)

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