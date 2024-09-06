from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union

def close_polygon(polygon, size, quad_segs=1):
    """
    Close a polygon using dilation and erosion.

    This algorithm relies on the successive dilation and erosion
    of polygon to merge close polygons together and simplify
    their complexity.

    Parameters
    ----------
    polygon : Polygon or MultiPolygon
        The polygon to close.
    size : float
        The size of the dilation and erosion.
    quad_segs : int
        The number of linear segments in a quarter circle
        when performing the buffer. If above 1, the result
        may have round corners unsuitable for buildings.

    Returns
    -------
    Polygon or MultiPolygon

    See Also
    --------
    open_polygon :
        Open a polygon using erosion and dilation.

    Examples
    --------
    >>> polygon = Polygon([(1, 1), (1, 2), (2, 2), (2, 0), (1, 1)])
    >>> close_polygon(polygon, 1)
    <POLYGON ((1.1508574748015223 1.2059024350708356, 1.5586467790767369 1.1617820026014096, 1.4159303795947187 1.5899312010474647, 1.1111173768873637 1.4375246996937872, 1.1508574748015223 1.2059024350708356))>
    """
    buffered = polygon.buffer(size, quad_segs)
    # then filter the buffer a little bit to avoid geometry problems
    filtered = buffered.simplify(1.0)
    # check if it is still a multi polygon after the buffer and the erosion
    if (polygon.geom_type == 'MultiPolygon'):
        return close_multipolygon(filtered, size, quad_segs)
    if (polygon.geom_type == 'Polygon'):
        closed = erosion(filtered, quad_segs)
        return closed

def open_polygon(polygon, size, quad_segs=1):
    """
    Open a polygon using erosion and dilation.

    This algorithm relies on the successive erosion and dilation
    of a polygon to get rid of thin sections.

    Parameters
    ----------
    polygon : Polygon or MultiPolygon
        The polygon to close.
    size : float
        The size of the erosion and dilation.
    quad_segs : int
        The number of linear segments in a quarter circle
        when performing the buffer. If above 1, the result
        may have round corners unsuitable for buildings.

    Returns
    -------
    Polygon or MultiPolygon

    See Also
    --------
    close_polygon :
        Close a polygon using dilation and erosion.

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 4), (4, 4), (4, 0), (0, 0)])
    >>> open_polygon(polygon, 1)
    <POLYGON ((0 1, 1 4, 4 3, 3 0, 0 1))>
    """
    
    polygons = []
    if (polygon.geom_type == 'Polygon'):
        polygons = open_simple(polygon, size, quad_segs)
    elif (polygon.geom_type == 'MultiPolygon'):
        for simple in polygon.geoms:
            polygons.extend(open_simple(simple, size, quad_segs))
            
    return MultiPolygon(polygons)

def close_multipolygon(multipolygon, buffer_size, quad_segs=1):
    buffered = multipolygon.buffer(buffer_size, quad_segs=quad_segs)
    # then filter the buffer a little bit to avoid geometry problems
    filtered = buffered.simplify(1.0)

    # check if it is still a multi polygon after the buffer and the erosion
    if (filtered.geom_type == 'MultiPolygon'):
        return erode_multipolygon(filtered, buffer_size, quad_segs)
    
    if (filtered.geom_type == 'Polygon'):
        return erosion(filtered, buffer_size, quad_segs)

    # if we arrive here, something went wrong during the closing, return the initial geometry
    return multipolygon

def open_simple(polygon, buffer_size, quad_segs=1):
    eroded = erosion(polygon, buffer_size, quad_segs)
    if (eroded is None):
        return []
    dilated = eroded.buffer(buffer_size, quad_segs)
    filtered = dilated.simplify(1.0)

    opened = []
    if filtered.geom_type == 'MultiPolygon':
        return [ s for s in filtered.geoms ]
    else:
        return [ filtered ]

def erode_multipolygon(multipolygon, buffer_size, quad_segs=1):
    erodedlist=[]
    for simple in multipolygon.geoms:
        eroded = erosion(simple, buffer_size, quad_segs)
        if eroded is None:
            continue
        erodedlist+=[eroded]
    return unary_union(erodedlist)

def erosion(polygon, buffer_size, quad_segs=1):
    # get the outer ring of the polygon
    outer_ring = polygon.exterior

    # first erode the outer ring
    eroded_outer = erosion_no_hole(Polygon(outer_ring), buffer_size, quad_segs)
    if(eroded_outer is None):
        return None
    
    eroded = eroded_outer
    # now handle the holes of the polygon
    if (len(polygon.interiors) > 0):
        for i in range(0,len(polygon.interiors)):
            ring = polygon.interiors[i]
            poly = Polygon(ring)
            buffered = poly.buffer(buffer_size, quad_segs)
            eroded = eroded.difference(buffered)
    
    if(eroded.is_empty):
        return None
    
    return eroded

# Applies an erosion to a polygon with no hole (or inner ring)
def erosion_no_hole(polygon, buffer_size, quad_segs=1):
    # get the outer ring of the polygon
    outer_ring = polygon.exterior

    # then buffer the outer ring of the polygon
    buffered_ring = outer_ring.buffer(buffer_size, quad_segs)

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