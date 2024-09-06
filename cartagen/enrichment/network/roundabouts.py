import geopandas as gpd
import numpy as np
from cartagen.utils.partitioning.network import network_faces

def detect_roundabouts(roads, area_threshold=40000, miller_index=0.95):
    """
    Detect roundabouts based on geometric properties.

    This algorithm proposed by Touya :footcite:p:`touya:2010` detects
    roundabouts inside a road network.

    Parameters
    ----------
    network : GeoDataFrame of LineString
        Road network to analyze.
    area_threshold : int, optional
        The area (in square meters) above which the object is not considered a roundabout.
    miller_index : float, optional
        Index of compactess that determines if the shape is round or not.

    Returns
    -------
    GeoDataFrame of Polygon

    See Also
    --------
    detect_branching_crossroads :
        Detect branching crossroads inside the road network.
    collapse_branching_crossroads :
        Collapse branching crossroads to a point.
    collapse_roundabouts :
        Collapse roundabouts to a point.

    References
    ----------
    .. footbibliography::
    """
    crs = roads.crs
    roads = roads.to_dict('records')

    network = []
    for road in roads:
        network.append(road['geometry'])

    faces = network_faces(network, convex_hull=False)

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
        return gpd.GeoDataFrame()

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