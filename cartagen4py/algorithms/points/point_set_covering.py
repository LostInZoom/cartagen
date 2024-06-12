from shapely.ops import triangulate, unary_union
from shapely.geometry import Polygon, Point, MultiPoint
from cartagen4py.utils.geometry.segment import *

def delaunay_concave_hull(points, min_length):
    """
    This algorithm computes a concave hull from a set of points. The algorithm first computes the Delaunay triangulation of the set of points,
    and then removes iteratively the boundary edges that are longer than a parameter.
    The algorithm was proposed by Duckham et al. (https://doi.org/10.1016/j.patcog.2008.03.023)

    Parameters
    ----------
    points : the list of points we want to cover with a concave hull
    min_length : Minimal length of a boundary edge of the triangulation to stop the algorithm.

    Examples
    --------
    >>> points = [Point(1,1),Point(6,6),Point(1,6),Point(6,5),Point(2,4),Point(2,1), Point(1,4)]
    >>> delaunay_concave_hull(points, 2.0)
    <POLYGON ((1 1, 1 4, 1 6, 6 6, 6 5, 2 4, 2 1, 1 1))>
    """
    multipoint = MultiPoint(points)
    triangles = triangulate(multipoint)
    hull = unary_union(triangles)

    removed = []
    shrinkable = True
    while(shrinkable):
        hull_ext = hull.exterior
        longest = -1
        max_length = min_length
        i = 0
        for triangle in triangles:
            if(i in removed):
                i += 1
                continue

            # first check: if the triangle is inside the hull, pass to next
            if(triangle.intersects(hull_ext)==False):
                i += 1
                continue

            # regularity check: if the three points of the triangle are on the boundary, removing this triangle would create an invalid geometry
            if(hull_ext.intersects(Point(triangle.exterior.coords[0])) and hull_ext.intersects(Point(triangle.exterior.coords[1])) 
               and hull_ext.intersects(Point(triangle.exterior.coords[2]))):
                i += 1
                continue
            
            # check that we do not create a multipolygon
            if(hull.difference(triangle).geom_type == 'MultiPolygon'):
                i += 1
                continue

            segments = get_segment_list(triangle)
            for segment in segments:
                segment_length = segment.length()
                if segment_length > max_length:
                    max_length = segment_length
                    longest = i
            i += 1
        
        if longest == -1:
            break
        
        hull = hull.difference(triangles[longest])
        removed.append(longest)
    
    return hull

