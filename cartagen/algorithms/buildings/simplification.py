from math import atan2, pi, sqrt
from shapely.geometry import Polygon, Point, LinearRing
from cartagen.utils.geometry.segment import get_segment_list_polygon

def simplify_building(building, edge_threshold, parallel_limit=20*pi/180, orthogonal_limit=20*pi/180):
    """
    Simplify buildings by removing edges.
    
    This algorithm proposed by Ruas :footcite:p:`ruas:1999` analyses
    the edges of the polygon to find the ones that should be removed and how they can be replaced.
    It was integrated in the AGENT project. Port of the CartAGen implementation of the algorithm.

    Parameters
    ----------
    building : Polygon
        The shapely building to be simplified.
    edge_threshold : float
        Minimum length of an edge to be considered by the simplification algorithm.
    parallel_limit : float, optional
        Limit angle to consider an edge into the parallel case of the simplification algorithm.
        The default value is set to :math:`20·pi/180`
    orthogonal_limit : float, optional
        Limit angle to consider an edge into the orthogonal case of the simplification algorithm.
        The default value is set to :math:`20·pi/180`
    
    Returns
    -------
    Polygon

    See Also
    --------
    square_polygon_ls :
        Squares polygon using the least squares method.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> building = Polygon([(0, 0), (0, 10), (2, 10), (2, 9), (10, 9), (10, 0), (0, 0)])
    >>> simplify_building(building, 2.5)
    <POLYGON ((0 0, 0 9.5, 10 9.5, 10 0, 0 0))>
    """

    # first get all the small segments in the polygon
    small_segments = __get_small_segments(building, edge_threshold)

    # loop on these small segments to remove them
    current_polygon = building
    while(len(small_segments)>0):
        # get the smallest in the list
        smallest = small_segments[0]
        smallest_index = 0
        for i in range (1,len(small_segments)):
            segment = small_segments[i]
            if(segment[0].length() < smallest[0].length()):
                smallest = segment
                smallest_index = i
        
        # remove the smallest segment from the list
        del small_segments[smallest_index]

        # try to remove the segment from the polygon
        ok, new_polygon = __delete_side_polygon(current_polygon, smallest[0], smallest[1], parallel_limit, orthogonal_limit)

        if (ok):
            current_polygon = new_polygon
            small_segments = __get_small_segments(current_polygon, edge_threshold)

    # sometimes the algorithm bugs and produces very large polygon. This verification prevents such problems
    if (current_polygon.area > 2 * building.area):
        return building
    
    return current_polygon

def __get_small_segments(polygon, edge_threshold):
    segment_list = get_segment_list_polygon(polygon)
    small_segments = []
    # retain only the segments smaller than the threshold
    for segment, id in segment_list:
        if(segment.length() < edge_threshold):
            small_segments.append((segment,id))
    return small_segments

def __delete_side_polygon(polygon, segment, ring_index, parallel_limit, orthogonal_limit):
    # first get the ring the segment to remove belongs to

    ring = polygon.exterior
    if(ring_index>=0):
        ring = polygon.interiors[ring_index]
    
    # if there are 4 vertices or less, no simplification possible
    if(len(ring.coords) <= 4):
        return False, polygon
    
    # get the vertices around the segment
    a = segment.point1
    b = segment.point2

    # get the vertices just before and after the segment
    index_a = list(ring.coords).index(a)
    index_a_ = index_a-1
    if(index_a == 0):
        index_a_ = len(ring.coords)-1
    index_b = list(ring.coords).index(b)
    index_b_ = index_b + 1
    if(index_b_ == len(ring.coords)):
        index_b_ = 0
    a_ = ring.coords[index_a_]
    b_ = ring.coords[index_b_]

    # compute the angle between (a, a_) and (b, b_) in the interval ]-pi, pi]
    angle = atan2(b_[1] - b[1], b_[0] - b[0]) - atan2(a_[1] - a[1], a_[0] - a[0])
    if(angle <= -pi):
        angle += 2 * pi
    elif angle > pi:
        angle -= 2 * pi

    # case where the segments around the one to remove are almost orthogonal. In this case, the segments are extended
    if(abs(angle) <= pi / 2 + orthogonal_limit and abs(angle) >= pi / 2 - orthogonal_limit):
    
        # compute the intersection between the segments
        xa = a_[0] - a[0]
        ya = a_[1] - a[1]
        xb = b_[0] - b[0]
        yb = b_[1] - b[1]
        if (xa * yb - ya * xb == 0):
            return False, polygon
        t = (xb * (a[1] - b[1]) - yb * (a[0] - b[0])) / (xa * yb - ya * xb)
        c_ = Point(a[0] + t * xa, a[1] + t * ya)

        # create the new sequence of coordinates where segment vertices are removed and point c_ is added
        new_coords = []
        for point in ring.coords:
            if (point != a and point != b):
                new_coords.append(point)
            elif (point == a):
                new_coords.append(c_.coords[0])
        
        new_ring = LinearRing(new_coords)
        return True, __update_ring_polygon(polygon, new_ring, ring_index)

    # case where the segments around the one to remove are almost parallel. In this case, it is a small protrusion that we want
    # to remove entirely. To do that, a_ is projected on (b, b_) and b_ is projected on (a, a_)
    if(abs(angle) < parallel_limit):
        if(a_ == b_):
            return False, polygon
        # this case removes two vertices so we need at least 6 in the ring
        if(len(ring.coords) == 5):
            return False, polygon

        if ((b[0] - b_[0]) * (b[0] - b_[0]) + (b[1] - b_[1]) * (b[1] - b_[1])) == 0:
            return False, polygon

        if ((a[0] - a_[0]) * (a[0] - a_[0]) + (a[1] - a_[1]) * (a[1] - a_[1])) == 0:
            return False, polygon

        # compute the projections
        aux_a = ((b[0] - b_[0]) * (a_[0] - b_[0])
          + (b[1] - b_[1]) * (a_[1] - b_[1]))  / ((b[0] - b_[0]) * (b[0] - b_[0]) + (b[1] - b_[1]) * (b[1] - b_[1]))
        ca = Point(b_[0] + aux_a * (b[0] - b_[0]),b_[1] + aux_a * (b[1] - b_[1]))
        app_a = (b_[0] - ca.x) * (b[0] - ca.x) + (b_[1] - ca.y) * (b[1] - ca.y) < 0
        aux_b = ((a[0] - a_[0]) * (b_[0] - a_[0])
          + (a[1] - a_[1]) * (b_[1] - a_[1])) / ((a[0] - a_[0]) * (a[0] - a_[0]) + (a[1] - a_[1]) * (a[1] - a_[1]))
        cb = Point(a_[0] + aux_b * (a[0] - a_[0]), a_[1] + aux_b * (a[1] - a_[1]))
        app_b = (a_[0] - cb.x) * (a[0] - cb.x) + (a_[1] - cb.y) * (a[1] - cb.y) < 0

        c1 = None
        c2 = None
        if(app_a == False and app_b == False):
            # none of the projected points belong to the segment. Very rare case.
            return False, polygon
        
        if(app_a == True and app_b == False):
            c1 = ca.coords[0]
            c2 = b_
        elif(app_b == True and app_a == False):
            c1 = a_
            c2 = cb.coords[0]
        else:
            da = Point(a_).distance(Point(a))
            db = Point(b_).distance(Point(b))
            if(da < db):
                c1 = ca.coords[0]
                c2 = b_
            else:
                c1 = a_
                c2 = cb.coords[0]

        # then create the new linear ring where a_, a, b, b_ are removed and replaced by c1 and c2
        new_coords = []
        for point in ring.coords:
            if (point != a and point != b and point != a_ and point != b_):
                new_coords.append(point)
            elif (point == a):
                new_coords.append(c1)
            elif (point == b):
                new_coords.append(c2)
        new_ring = LinearRing(new_coords)
        return True, __update_ring_polygon(polygon, new_ring, ring_index)
    
    # last case where the three segments are replaced by a mean segment
    if(abs(angle) > pi - parallel_limit):
        # first compute the orientation of this new segment
        dx = a[0] - a_[0] + b_[0] - b[0]
        dy = a[1] - a_[1] + b_[1] - b[1]
        n = sqrt(dx * dx + dy * dy)
        dx /= n
        dy /= n
        # then compute the middle of segment (a_,b_)
        xm = (a_[0] + b_[0]) * 0.5
        ym = (a_[1] + b_[1]) * 0.5

        # then project a_ and b_ on this new straight line
        aux1 = (a_[0] - xm) * dx + (a_[1] - ym) * dy
        ca = Point(xm + aux1 * dx, ym + aux1 * dy)
        aux2 = (b_[0] - xm) * dx + (b_[1] - ym) * dy
        cb = Point(xm + aux2 * dx, ym + aux2 * dy)

        # then create the new linear ring where a_, a, b, b_ are removed and replaced by ca and cb
        new_coords = []
        for point in ring.coords:
            if (point != a and point != b and point != a_ and point != b_):
                new_coords.append(point)
            elif (point == a):
                new_coords.append(ca.coords[0])
            elif (point == b):
                new_coords.append(cb.coords[0])
        new_ring = LinearRing(new_coords)
        return True, __update_ring_polygon(polygon, new_ring, ring_index)
        
    return False, polygon

def __update_ring_polygon(polygon, new_ring, ring_index):
    if ring_index == -1:
        return Polygon(new_ring,polygon.interiors)
    else:
        interiors=[]
        for i in range(0,len(polygon.interiors)):
            if i == ring_index:
                interiors.append(new_ring)
            else:
                interiors.append(polygon.interiors[i])
        return Polygon(polygon.exterior, interiors)