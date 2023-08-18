import shapely

def distance_area(geom1, geom2):
    """
    Calculcate the distance area between two geometries.
    distance area = area(intersection) / area(union)
    """
    if geom1.geom_type != 'Polygon' or geom2.geom_type != 'Polygon':
        raise Exception('Distance area calculation can be calculated only on Polygon geometries.')

    intersection = shapely.intersection(geom1, geom2)
    union = shapely.union(geom1, geom2)

    if intersection.is_empty == False and union.is_empty == False:
        return 1 - (intersection.area / union.area)
    else:
        return None