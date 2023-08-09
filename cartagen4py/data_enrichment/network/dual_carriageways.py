import geopandas as gpd
import numpy as np
import shapely

from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

def detect_dual_carriageways(
        network, importance=None, value=None,
        concavity=0.85, elongation=6.0, compactness=0.12,
        area=60000.0, width=20.0, huber=16
    ):
    """
    Detect dual carriageways and return road separators.
    Return None if none were found.
    Parameters
    ----------
    network : geopandas GeoDataFrame of LineStrings
        The road network to analyze.
    importance : str optional.
        The attribute name of the data on which road importance is based.
        Default value is set to None which means every road is taken for the network face calculation.
    value : int optional.
        The maximum value of the importance attribute. Roads with an importance higher than this value will not be taken.
        Default value is set to None.
    """
    crs = network.crs

    network = network.to_dict('records')

    if (importance is not None and value is None) or (importance is None and value is not None):
        raise Exception("Provide both arguments (importance, value) or none.")

    attribute = False
    if importance is not None and value is not None:
        attribute = True

    roads = []
    for road in network:
        if attribute:
            if int(road[importance]) <= value:
                roads.append(road['geometry'])
        else:
            roads.append(road['geometry'])

    faces = calculate_network_faces(roads, convex_hull=False)
    tree = shapely.STRtree(roads)

    separators = []
    index = 0
    for face in faces:
        add, infos = is_dual_carriageway(
            face, roads, tree,
            concavity=concavity,
            elongation=elongation,
            compactness=compactness,
            area=area, width=width, huber=huber
        )
        if add:
            infos['cid'] = index
            separators.append(infos)
            index += 1

    if len(separators) > 0:
        return gpd.GeoDataFrame(separators, crs=crs)
    else:
        return None

def is_dual_carriageway(
        polygon, roads, tree,
        concavity=0.85, elongation=6.0, compactness=0.12,
        area=60000.0, width=20.0, huber=16
    ):
    """
    Return True or False if the network face is detected as a dual carriageway separator along with infos.
    """

    face = NetworkFace(polygon)

    infos = {
        'area': face.area,
        'perimeter': face.perimeter,
        'concavity': face.concavity,
        'elongation': face.elongation,
        'compactness': face.compactness,
        'length': face.length,
        'width': face.width,
        'huber': face.huber,
        'geometry': polygon
    }

    # Check if the face is convex...
    if face.concavity > concavity:
        # False if the width is larger than the threshold
        if face.width > width:
            return False, infos
        # Check the elongation only...
        if (face.elongation > elongation):
            return True, infos
        # ...or the compactness with the elongation.
        if (face.compactness < compactness) and (face.elongation > (elongation / 2)):
            return True, infos
    # ...or not convex.
    else:
        # Check the compactness and the surface area
        if (face.compactness < compactness) and (face.area < area):
            # If the compactness is above half the limit compactness
            if face.compactness > (compactness / 2):
                # Check that the huber width is below the threshold
                if face.huber > huber:
                    return False, infos
            return True, infos 

    # Special case of long motorways, i.e. large area with small compactness
    if (face.compactness < (compactness / 4)) and (face.area < (10 * area)):
        return True, infos

    # Here, it's not a dual carriageway
    return False, infos