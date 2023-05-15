# This file contains measures on buildings that can be used in constraints or inside algorithms.
from cartagen4py.utils.geometry.segment import *
from cartagen4py.utils.geometry.line import to_2d
from shapely.ops import nearest_points, transform

def building_min_width(building):
    min_width = building.length
    segments = get_segment_list(building)
    # loop on the segments in the list
    for i in range (0,len(segments)):
        segment = segments[i]

        other_segments = get_segment_list(building)
        for seg in other_segments:
            if(seg.point1==segment.point1 or seg.point2==segment.point1 or seg.point2==segment.point2 or seg.point1==segment.point2):
                other_segments.remove(seg)

        for other_segment in other_segments:
            nearest_pts = nearest_points(segment.get_geom().centroid,other_segment.get_geom().centroid)
            dist = Point(nearest_pts[0]).distance(Point(nearest_pts[1]))
            point1_2d = nearest_pts[0]
            point2_2d = nearest_pts[1]
            if(point1_2d.has_z):
                point1_2d = transform(to_2d, point1_2d)
            if(point2_2d.has_z):
                point2_2d = transform(to_2d, point2_2d)
            nearest_connection = LineString([point1_2d,point2_2d])
            if(dist >= min_width):
                continue
            if (building.contains(nearest_connection) == False):
                continue
            if(building.exterior.contains(nearest_connection)):
                continue
            min_width = dist
        
    return min_width