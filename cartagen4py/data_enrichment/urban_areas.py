from shapely.geometry import Polygon
from cartagen4py.utils.math import morphology
from shapely.ops import unary_union

# computes the urban areas from a set of buildings, using a method from Boffet
def compute_boffet_urban_areas(buildings, dilation_size, erosion_size, simplification_distance = 2.0):
    buflist=[]
    for building in buildings:
        buflist += [building.buffer(dilation_size)]
    dilation = unary_union(buflist)
    simplified = dilation.simplify(simplification_distance)
    closed = morphology.closing_multi_polygon(simplified, erosion_size)
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


