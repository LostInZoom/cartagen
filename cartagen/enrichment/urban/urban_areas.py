from shapely.geometry import Polygon
from cartagen.utils.math.morphology import close_multipolygon
from shapely.ops import unary_union

def boffet_areas(buildings, buffer, erosion, simplification_distance=2.0):
    """
    Calculate urban areas from buildings.

    This algorithm proposed by Boffet :footcite:p:`boffet:2003` uses
    buffer around the buildings, then simplify and erode the unioned result
    to characterize urban areas.

    Parameters
    ----------
    polygons : list of Polygon
        Buildings to generate the urban area from.
    buffer : float
        The buffer size used to merge buildings that are close from each other.
    erosion : float
        The erosion size to avoid the urban area to expand
        too far from the buildings located on the edge.
    simplification_distance : float, optional
        The distance threshold used by the
        Douglas-Peucker simplification on the edge.

    Returns
    -------
    list of Polygon

    See Also
    --------
    morphological_amalgamation :
        Amalgamate buildings using dilation and erosion. Useful for larger scale maps.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> buildings = [Polygon([(1, 0), (9, 0), (9, 6), (1, 6), (1, 0)]), Polygon([(10, 0), (17, 0), (17, 6), (10, 6), (10, 0)])]
    >>> boffet_areas(buildings, 2.0, 1.0)
    <[<POLYGON ((-0.909 5.627, 17.078 7.784, 18.66 -0.02, 0.874 -1.782, -0.909 5.6...>]>
    """

    buflist=[]
    # Stores buffer for each building using the dilation size
    for building in buildings:
        buflist += [building.buffer(buffer)]
    dilation = unary_union(buflist)
    simplified = dilation.simplify(simplification_distance)
    closed = close_multipolygon(simplified, erosion)
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


