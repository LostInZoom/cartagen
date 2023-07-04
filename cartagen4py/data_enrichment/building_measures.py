# This file contains measures on buildings that can be used in constraints or inside algorithms.
from cartagen4py.utils.geometry.segment import *
from cartagen4py.utils.geometry.line import to_2d, get_nearest_vertex
from cartagen4py.utils.geometry.angle import angle_3_pts, angle_to_zero_pi
from shapely.ops import nearest_points, transform, triangulate
from shapely.geometry import MultiPolygon, MultiPoint
from math import pi

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

def block_triangulation(buildings, networks, max_distance):
    points = []
    objects = []
    # first, add the building centroids in the triangulation
    for building in buildings:
        points.append(building.centroid)
        objects.append(building)

    # then add the closest point to a building for each network section
    for section in networks:
        min_dist = max_distance
        nearest = None
        # find the nearest point on section to a building from buildings
        for building in buildings:
            nearest_point = nearest_points(section, building)
            dist = nearest_point[0].distance(nearest_point[1])
            if(dist < min_dist):
                nearest = nearest_point[0]
                min_dist = dist

        # arrived here, add nearest_point to the list of triangulation points
        if min_dist < max_distance:
            points.append(nearest)
            objects.append(section)
    
    # compute the delaunay triangulation
    raw_tri = triangulate(MultiPoint(points), edges=True)

    # filter the edges of the triangulation to remove undesired edges
    triangulation = []
    multi_buildings = MultiPolygon(buildings)
    for edge in raw_tri:
        point1 = Point(edge.coords[0])
        point2 = Point(edge.coords[1])
        is_building1 = multi_buildings.contains(point1)
        is_building2 = multi_buildings.contains(point2)

        # edges between two buildings are kept
        if is_building1 and is_building2:
            building1 = objects[points.index(point1)]
            building2 = objects[points.index(point2)]
            triangulation.append([edge, ((building1,0),(building2,0))])
            continue

        # edges between two roads are removed
        if is_building1 == False and is_building2==False:
            continue

        # here we have an edge between a road and a building
        # we want to keep only the edges that form an almost right angle with the road (the others are considered as redundant)
        if is_building1:
            # first we get the geometry of the section
            building = objects[points.index(point1)]
            section = objects[points.index(point2)]
            vertex_on_section = get_nearest_vertex(section,point2)
            angle = angle_3_pts(vertex_on_section,point2,point1)
            if abs(pi/2 - angle_to_zero_pi(angle)) < 12.5:
                triangulation.append([edge, ((building,0),(section,1))])
                continue
        if is_building2:
            # first we get the geometry of the section
            building = objects[points.index(point2)]
            section = objects[points.index(point1)]
            vertex_on_section = get_nearest_vertex(section,point1)
            angle = angle_3_pts(vertex_on_section,point1,point2)
            if abs(pi/2 - angle_to_zero_pi(angle)) < 12.5:
                triangulation.append([edge, ((section,1),(building,0))])
                continue

    return triangulation

# for a triangulation of a city block, this retrieves the edges of the triangulation that are connected to a given building
def triangulation_edges_around_building(building, triangulation):
    
    edges = []
    for edge in triangulation:
        objects = edge[1]
        if objects[0][1]==0 and objects[0][0].equals(building):
            edges.append((edge[0], objects[1]))
        if objects[1][1]==0 and objects[1][0].equals(building):
            edges.append((edge[0], objects[0]))
    return edges

def building_congestion(building, triangulation, max_distance, nb_orientations=16):
    congested = []
    nb_edges_2buildings = 0
    for i in range(0,nb_orientations):
        congested.append(0)
    step = 2* pi / nb_orientations

    # get the edges of the triangulation connected to building
    edges = triangulation_edges_around_building(building, triangulation)

    # case where all edges have been removed around building
    if len(edges) == 0:
        return congested
    
    for edge in edges:
        if(edge[0].length > max_distance):
            continue

        # count edges with buildings
        if edge[1][1] == 0:
            nb_edges_2buildings += 1
        
        # compute the congestion for this edge
        congestion = (max_distance-edge[0].length )/max_distance

        # now compute the direction of the congestion
        point1 = Point(edge[0].coords[0])
        point2 = Point(edge[0].coords[1])
        angle = angle_2_pts(point1, point2)
        if angle < 0:
            angle = angle + 2 * pi
        for i in range(0,nb_orientations):
            orientation = i * step
            # compute the difference between angle and orientation
            diff = abs(angle_to_zero_pi(angle)-orientation)
            if diff > pi/4:
                continue
            projected_congestion = congestion * (1 - 4*diff/pi)
            congested[i] = projected_congestion
        
    return congested

# computes the compactness of a polygon using the Miller index
def polygon_compactness(polygon):
    return 4 * pi * polygon.area / (polygon.length * polygon.length)

# computes the elongation of a polygon
def polygon_elongation(polygon):
    mbr = polygon.minimum_rotated_rectangle
    coords = mbr.exterior.coords
    length1 = Point(coords[0]).distance(Point(coords[1]))
    length2 = Point(coords[1]).distance(Point(coords[2]))
    if(length1 < length2):
        return length2 / length1
    else:
        return length1 / length2