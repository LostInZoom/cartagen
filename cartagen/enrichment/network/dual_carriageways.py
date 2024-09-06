import geopandas as gpd
import shapely

from cartagen.utils.partitioning.network import network_faces
from cartagen.utils.network.faces import NetworkFace

def detect_dual_carriageways(
        roads, importance=None, value=None,
        concavity=0.85, elongation=6.0, compactness=0.12,
        area=60000.0, width=20.0, huber=16
    ):
    """
    Detect dual carriageways based on geometric properties.

    This algorithm proposed by Touya :footcite:p:`touya:2010`
    detects the network faces as road separator (*i.e.* separation between
    dual carriageways) when the polygon meets the geometric requirements.
    Those values can be tweaked to fine-tune the detection, but complex interchange will
    nonetheless cause wrong characterization.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        Road network to analyze.
    importance : str, optional
        The attribute name of the data on which road importance is based.
        Default value is set to None which means every road is taken for the network face calculation.
    value : int, optional
        Maximum value of the importance attribute.
        Roads with an importance higher than this value will not be taken.
    concavity : float, optional
        Maximum concavity.
    elongation : float, optional
        Minimum elongation.
    compactness : float, optional
        Maximum compactness.
    area : float, optional
        Area factor to detect very long motorways.
    width : float, optional
        Maximum width of the the :func:`minimum_rotated_rectangle <shapely.minimum_rotated_rectangle>`.
    huber : int, optional
        Huber width for long motorways.

    See Also
    --------
    collapse_dual_carriageways :
        Collapse dual carriageways using a TIN skeleton.

    Notes
    -----
    - **concavity** is the area of the polygon divided by the area of its convex hull.
    - **elongation** is the length of the :func:`minimum_rotated_rectangle <shapely.minimum_rotated_rectangle>`
      divided by its width.
    - **compactness** is calculated using :math:`(4·pi·area)/(perimeter^2)`

    References
    ----------
    .. footbibliography::
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

    faces = network_faces(network, convex_hull=False)
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
        return gpd.GeoDataFrame()
    
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