# This file contains measures on buildings that can be used in constraints or inside algorithms.
from cartagen.utils.geometry.segment import get_segment_list
from cartagen.utils.geometry.line import to_2d, get_nearest_vertex, resample_line
from cartagen.utils.geometry.angle import angle_3_pts, angle_2_pts, angle_to_zero_pi, angle_between_2lines
from shapely.ops import nearest_points, transform, triangulate, substring, split
from shapely.geometry import MultiPolygon, MultiPoint, Polygon
from math import pi

def building_min_width(building):
    """
    Returns the minimum width inside a building. The minimum width is the minimum distance between
    two edges of the buildings that are not adjacent. 
    The measure was proposed during the AGENT project. 'building' should be a shapely ''Polygon''.
    """
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
        is_building1 = False
        is_building2 = False
        if multi_buildings.is_valid:
            is_building1 = multi_buildings.contains(point1)
            is_building2 = multi_buildings.contains(point2)
        else:
            for b in buildings:
                if b.contains(point1):
                    is_building1 = True
                if b.contains(point2):
                    is_building2 = True

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

# this function computes the congestion around a building following the principles proposed in Anne Ruas PhD.
# It returns an array of congested direction around the building, with a zero value when that direction is not congested.
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
        return (congested, 0)
    
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

    return (congested, nb_edges_2buildings)
    
# Compute the corner buildings in a block.
# The list contains the ids of the buildings
# (i.e. their number in the list of buildings),
# not the buildings themselves.
def corner_buildings(buildings, networks, angle_tolerance=22.5, triangle_edge = 20.0):
    corner = []
    radians_tol = angle_tolerance * pi / 180.0
    #road_parts = __cut_roads_at_corners(networks, radians_tol)
    # first create the corner areas
    corner_areas = []
    road_part_id = 0
    processed_roads = []
    for road_part in networks:
        processed_roads.append(road_part_id)
        # get the connected roads
        connected_roads = __get_corner_connected_roads(road_part_id, networks)
        print(len(connected_roads))
        for connected_road in connected_roads:
            if connected_road in processed_roads:
                continue
            # compute the angle between road_part and connected_road at their intersection point
            angle = angle_to_zero_pi(angle_between_2lines(road_part, networks[connected_road]))
            if angle > (pi / 2 - radians_tol) and angle < (pi / 2 + radians_tol):
                print("pass the angle test")
                # build a triangle as a new corner area
                # first get the three vertices of the triangle
                other_road = networks[connected_road]
                a = road_part.intersection(other_road)
                distab = triangle_edge * min(angle, pi/2) / max(angle, pi/2)
                b = road_part.interpolate(distab)
                if a.equals(Point(road_part.coords[len(road_part.coords)-1])):
                    reversed_road = substring(road_part, 1, 0, normalized=True)
                    b = reversed_road.interpolate(distab)
                c = other_road.interpolate(distab)
                if a.equals(Point(other_road.coords[len(other_road.coords)-1])):
                    reversed_road = substring(other_road, 1, 0, normalized=True)
                    c = reversed_road.interpolate(distab)
                corner_areas.append(Polygon([a,b,c,a]))
        road_part_id += 1

    # now get the buildings intersecting one of the corner areas
    print("corner areas:")
    print(len(corner_areas))
    id = 0
    for building in buildings:
        for triangle in corner_areas:
            if building.intersects(triangle):
                corner.append(id)
                break
        id += 1
    return corner, corner_areas

def __cut_roads_at_corners(networks, angle_tolerance):
    road_parts = []

    for road in networks:
        cumulated_angle = 0.0
        ante_pt = None
        prev_pt = None
        initial_pt = road.coords[0]
        nb_pts_itr = 0
        # simplify the line to get sharper angles
        simp_road = road.simplify(1.0)
        # densify the simplified road
        dense_road = resample_line(simp_road, 10.0)
        for pt in dense_road.coords:
            if nb_pts_itr > 3:
                # add the angle between [0,Pi/2] to the cumulated angle
                angle = angle_to_zero_pi(angle_3_pts(Point(ante_pt), Point(prev_pt), Point(pt)))
                if angle > pi / 2:
                    angle = abs(angle - pi)
                cumulated_angle += angle
                if(cumulated_angle > pi/2 - angle_tolerance):
                    # cut the road at prev_pt
                    last_portion = dense_road
                    if initial_pt != road.coords[0]:
                        last_portion = split(dense_road, Point(initial_pt)).geoms[1]
                    split_road = split(last_portion,Point(prev_pt)).geoms[0]  
                    road_parts.append(split_road)                
                    cumulated_angle = 0
                    nb_pts_itr = 0
                    initial_pt = prev_pt

            nb_pts_itr += 1
            ante_pt = prev_pt
            prev_pt = pt
        
        # append the final part of the road
        last_portion = dense_road
        if initial_pt != road.coords[0]:
            last_portion = split(dense_road, Point(initial_pt)).geoms[1]
        road_parts.append(last_portion)

    return road_parts

def __get_corner_connected_roads(road_id, networks):
    connected_roads = []
    road_geom = networks[road_id]
    id = 0
    for road in networks:
        if id == road_id:
            id += 1
            continue
        if road_geom.touches(road):
            connected_roads.append(id)
        id += 1

    return connected_roads