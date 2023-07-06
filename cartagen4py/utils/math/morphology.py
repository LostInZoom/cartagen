# this file contains algorithms of morphological maths, i.e. derivative of dilation and erosion of lines and polygons

from shapely.geometry import Point, Polygon, MultiPolygon

def closing_multi_polygon(multipolygon, buffer_size, cap_style=1):
    buffered = multipolygon.buffer(buffer_size, cap_style=cap_style)
    # then filter the buffer a little bit to avoid geometry problems
    filtered = buffered.simplify(1.0)

    # check if it is still a multi polygon after the buffer and the erosion
    if (filtered.geom_type == 'MultiPolygon'):
        return erosion_multipolygon(filtered, buffer_size, cap_style)
    
    if (filtered.geom_type == 'Polygon'):
        return erosion(filtered, buffer_size, cap_style)

    # if we arrive here, something went wrong during the closing, return the initial geometry
    return multipolygon

# Applies a morphological closing on a polygon (dilation+erosion)
def closing(polygon, buffer_size, cap_style=1):
    buffered = polygon.buffer(buffer_size, cap_style)
    # then filter the buffer a little bit to avoid geometry problems
    filtered = buffered.simplify(1.0)
    closed = erosion(filtered, cap_style)
    return closed

def opening(polygon, buffer_size, cap_style=1):
    if (polygon.geom_type == 'Polygon'):
        return opening_simple(polygon, buffer_size, cap_style)

    polygons = []
    if (polygon.geom_type == 'MultiPolygon'):
        for simple in polygon.geoms:
            polygons.append(opening_simple(simple, buffer_size, cap_style))
    return MultiPolygon(polygons)

def opening_simple(polygon, buffer_size, cap_style=1):
    eroded = erosion(polygon, buffer_size, cap_style)

    if (eroded is None):
        return None
    
    dilated = eroded.buffer(buffer_size, cap_style)
    filtered = dilated.simplify(1.0)

    return filtered

def erosion_multipolygon(multipolygon, buffer_size, cap_style=1):
    polygons = []

    for simple in multipolygon.geoms:
        eroded = erosion(simple, buffer_size, cap_style)
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

def erosion(polygon, buffer_size, cap_style=1):
    # get the outer ring of the polygon
    outer_ring = polygon.exterior

    # first erode the outer ring
    eroded_outer = erosion_no_hole(Polygon(outer_ring), buffer_size, cap_style)
    if(eroded_outer is None):
        return None
    
    eroded = eroded_outer
    # now handle the holes of the polygon
    if (len(polygon.interiors) > 0):
        for i in range(0,len(polygon.interiors)):
            ring = polygon.interiors[i]
            poly = Polygon(ring)
            buffered = poly.buffer(buffer_size, cap_style)
            eroded = eroded.difference(buffered)
    
    if(eroded.is_empty):
        return None
    
    return eroded



# Applies an erosion to a polygon with no hole (or inner ring)
def erosion_no_hole(polygon, buffer_size, cap_style=1):
    # get the outer ring of the polygon
    outer_ring = polygon.exterior

    # then buffer the outer ring of the polygon
    buffered_ring = outer_ring.buffer(buffer_size, cap_style)

    # check if the buffer is a valid polygon
    if(buffered_ring.geom_type != 'Polygon'):
        return None

    # at this point, we are interested in the inner rings of buffer
    num_rings = len(buffered_ring.interiors)
    
    # first, check if there is no ring
    if num_rings == 0:
        return None

    # then, if there is only one inner ring, return it as the output polygon
    if num_rings == 1:
        ring = buffered_ring.interiors[0]
        point = ring.coords[0]
        if(polygon.contains(Point(point))):
            return Polygon(ring)

    # then, if there are several inner rings, we have to join them as a multi-polygon
    rings = []
    for i in range(0,num_rings):
        ring = buffered_ring.interiors[i]
        poly = Polygon(ring)
        if(polygon.contains(poly)):
            rings.append(poly)
    
    return MultiPolygon(rings)



if __name__ == '__main__':
    from shapely.wkt import loads
    from shapely.geometry import Polygon, LineString
    import matplotlib.pyplot as plt
    import geopandas as gpd

    polygon1 = Polygon([(0,0),(0,100),(50,150),(100,100),(100,0),(0,0)])

    eroded = erosion_no_hole(polygon1,10)
    p1 = gpd.GeoSeries(polygon1)
    p2 = gpd.GeoSeries(eroded)
    base = p1.plot()
    p2.plot(ax=base, color='red')
    plt.show()
    print(eroded)