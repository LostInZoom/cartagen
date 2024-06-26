from shapely.geometry import Polygon
from cartagen4py.utils.math import morphology
from shapely.ops import unary_union

def boffet_areas(buildings, buffer, erosion, simplification_distance=2.0):
    """
    Calculate urban areas from buildings (Boffet, 2003).

    Buffer each provided buildings, simplify and erode the unioned result
    to generate a representation of the urbanized area.

    Parameters
    ----------
    polygons : list of shapely.Polygon
        Buildings to generate the urban area from.
    buffer : float
        The buffer size used to merge buildings that are close from each other.
    erosion : float
        The erosion size to avoid the urban area to expand
        too far from the buildings located on the edge.
    simplification_distance : float, Default=2.0
        The distance threshold used by the
        Douglas-Peucker simplification on the edge.

    Returns
    -------
    list of shapely.Polygon

    See Also
    --------
    morphological_amalgamation :
        Amalgamate buildings using dilation and erosion. Useful for larger scale maps.
    """

    buflist=[]
    # Stores buffer for each building using the dilation size
    for building in buildings:
        buflist += [building.buffer(buffer)]
    dilation = unary_union(buflist)
    simplified = dilation.simplify(simplification_distance)
    closed = morphology.closing_multi_polygon(simplified, erosion)
    if(closed is None):
        return []
    simplified = closed.simplify(simplification_distance)
    final = simplified.buffer(0.0)
    
    urban_areas = []
    if (final.geom_type == "Polygon"):
        urban_areas.append(Polygon(final.exterior))
        return urban_areas

    for geom in final.geoms:
        urban_areas.append(Polygon(geom.exterior))
    
    return urban_areas


