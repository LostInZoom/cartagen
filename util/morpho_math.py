# this file contains algorithms of morphological maths, i.e. derivative of dilation and erosion of lines and polygons

import shapely
from shapely import Polygon

def closing_multi_polygon(multipolygon, buffer_size, cap_style="round"):
    buffered = shapely.buffer(multipolygon, buffer_size, cap_style=cap_style)
    # then filter the buffer a little bit to avoid geometry problems
    filtered = shapely.simplify(buffered,1.0)

    # check if it is still a multi polygon after the buffer and the erosion
    if (filtered.geom_type == 'MultiPolygon'):
        return erosion_multipolygon(filtered, buffer_size, cap_style=cap_style)
    
    if (filtered.geom_type == 'Polygon'):
        return erosion(filtered, buffer_size, cap_style=cap_style)

    # if we arrive here, something went wrong during the closing, return the initial geometry
    return multipolygon

# Applies a morphological closing on a polygon (dilation+erosion)
def closing(polygon, buffer_size, cap_style="round"):
    buffered = shapely.buffer(polygon, buffer_size, cap_style=cap_style)
    # then filter the buffer a little bit to avoid geometry problems
    filtered = shapely.simplify(buffered,1.0)
    closed = erosion(filtered, cap_style=cap_style)
    return closed

def opening(polygon, buffer_size, cap_style="round"):
    eroded = erosion(polygon, buffer_size, cap_style=cap_style)

    if (eroded is None):
        return None
    
    dilated = shapely.buffer(eroded, buffer_size, cap_style=cap_style)
    filtered = shapely.simplify(dilated,1.0)

    return filtered

def erosion_multipolygon(multipolygon, buffer_size, cap_style="round"):
    polygons = []

    for simple in multipolygon.geoms:
        eroded = erosion(simple, cap_style=cap_style)
        if eroded is None:
            continue
        if(eroded.geom_type == 'Polygon'):
            polygons.append(eroded)
            continue
        if(eroded.geom_type == 'MultiPolygon'):
            for simple_eroded in eroded.geoms:
                polygons.append(simple_eroded)

    if(len(polygons)==0):
        return None

    return MultiPolygon(polygons)

def erosion(polygon, buffer_size, cap_style="round"):
    # get the outer ring of the polygon
    outer_ring = shapely.get_exterior_ring(polygon)

    # first erode the outer ring
    eroded_outer = erosion_no_hole(Polygon(outer_ring), cap_style=cap_style)
    if(eroded_outer is None):
        return None
    
    # now handle the holes of the polygon
    if (shapely.get_num_interior_rings(polygon) > 0):
        for i in range(0,shapely.get_num_interior_rings(polygon)):
            ring = shapely.get_interior_ring(buffer, i)
            poly = Polygon(ring)
            buffered = shapely.buffer(poly, buffer_size, cap_style=cap_style)
            eroded = shapely.difference(eroded,buffered)
    
    if(shapely.is_empty(eroded)):
        return None
    
    return eroded



# Applies an erosion to a polygon with no hole (or inner ring)
def erosion_no_hole(polygon, buffer_size, cap_style="round"):
    # get the outer ring of the polygon
    outer_ring = shapely.get_exterior_ring(polygon)

    # then buffer the outer ring of the polygon
    buffer = shapely.buffer(outer_ring, buffer_size, cap_style=cap_style)

    # check if the buffer is a valid polygon
    if(buffer.geom_type != 'Polygon'):
        return None

    # at this point, we are interested in the inner rings of buffer
    num_rings = shapely.get_num_interior_rings(buffer)

    # first, check if there is no ring
    if num_rings == 0:
        return None

    # then, if there is only one inner ring, return it as the output polygon
    if num_rings == 1:
        ring = shapely.get_interior_ring(buffer, 0)
        point = ring.coords[0]
        if(shapely.contains(polygon,point)):
            return Polygon(ring)

    # then, if there are several inner rings, we have to join them as a multi-polygon
    rings = []
    for i in range(0,num_rings-1):
        ring = shapely.get_interior_ring(buffer, i)
        poly = Polygon(ring)
        if(shapely.contains(polygon,poly)):
            rings.append(poly)

    return MultiPolygon(rings)



if __name__ == '__main__':
    from shapely.wkt import loads
    from shapely.geometry import Polygon, LineString

    polygon1 = Polygon([(0,0),(0,100),(50,150),(100,100),(100,0),(0,0)])

    eroded = erosion_no_hole(polygon1,10)
    print(eroded)