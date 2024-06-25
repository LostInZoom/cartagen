import geopandas as gpd
import numpy as np
import shapely

from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *

def detect_dual_carriageways(
        roads, importance=None, value=None,
        concavity=0.85, elongation=6.0, compactness=0.12,
        area=60000.0, width=20.0, huber=16
    ):
    """
    Detect dual carriageways based on geometric properties (Renard, 2009).

    This function detects the network faces as road separator (i.e. separation between
    dual carriageways) when the polygon meets the geometric requirements.
    Those values can be tweaked to fine-tune de detection, but complex interchange will
    nonetheless cause wrong characterization.

    Parameters
    ----------
    roads : GeoPandas.GeoDataFrame with LineString geometries
        Road network to analyze.
    importance : str, Default=None
        The attribute name of the data on which road importance is based.
        Default value is set to None which means every road is taken for the network face calculation.
    value : int, Default=None
        Maximum value of the importance attribute. Roads with an importance higher than this value will not be taken.
    concavity : float, Default=0.85
        Minimum concavity above which the face is a dual carriageway.
        It represents the factor between the polygon surface and its convex hull surface.
    elongation : float, Default=6.0
        Minimum elongation above which the face is a dual carriageway.
        It represents the ratio between the length and the
        width of the minimum rotated rectangle containing the polygon.
    compactness : float, Default=0.12
        Maximum compactness below which the face is a dual carriageway.
        (4*pi*area/perimeter^2)
    area : float, Default=60000.0
        Area factor to detect very long motorways.
    width : float, Default=20.0
        Minimum width above which the face is a dual carriageway.
        It represents the width of the minimum rotated rectangle containing the polygon.
    huber : int, Default=16
        Huber width for long motorways.

    See Also
    --------
    collapse_dual_carriageways
    """
    crs = roads.crs

    roads = roads.to_dict('records')

    if (importance is not None and value is None) or (importance is None and value is not None):
        raise Exception("Provide both arguments (importance, value) or none.")

    attribute = False
    if importance is not None and value is not None:
        attribute = True

    network = []
    for road in roads:
        if attribute:
            if int(road[importance]) <= value:
                network.append(road['geometry'])
        else:
            network.append(road['geometry'])

    faces = calculate_network_faces(network, convex_hull=False)
    tree = shapely.STRtree(network)

    separators = []
    index = 0
    for face in faces:
        add, infos = is_dual_carriageway(
            face, network, tree,
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