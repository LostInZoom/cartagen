import shapely
import geopandas as gpd
import numpy as np
from cartagen4py.utils.partitioning import *
from cartagen4py.utils.network import *
from cartagen4py.utils.geometry.distances import *


def is_branching_crossroad(polygon, roads, area_threshold,
        maximum_distance_area, roundabouts=None,
        allow_middle_node=False, allow_single_4degree_node=True
    ):
    """
    Return True or False whether the given polygon is a branching crossroad or not depending on the given parameters.
    """

    def is_triangular(dist_area, polygon, nodes, middle=None):
        """
        Check if the network face is triangular given a maximum allowed distance area.
        """
        tnodes = []
        for i, node in enumerate(nodes):
            if middle != i:
                tnodes.append(node[1])

        tnodes.append(tnodes[0])
        triangle = shapely.Polygon(tnodes)

        da = distance_area(polygon, triangle)

        if da <= dist_area:
            return True, da
        
        return False, None

    tree = shapely.STRtree(roads)

    infos = { 'geometry': polygon }

    area = polygon.area
    # Check if the area of the polygon is larger than the threshold
    if area < area_threshold:
        # Retrieve objects that intersects the considered network face using strtree
        intersects = tree.query(polygon)
        lines = []
        for i in intersects:
            l = roads[i]
            # Make an other test to really keep only intersecting roads, spatial index strtree using bbox
            if shapely.intersects(polygon, l):
                lines.append(roads[i])

        # Continue only if more than two lines are selected, otherwise it can't be a branching crossroad
        if len(lines) > 2:
            # Fully dissolve and node the subnetwork
            unioned = unary_union(lines)
            # Merge all contiguous lines
            merged = linemerge(unioned)

            network = []
            if merged.geom_type == 'LineString':
                network.append(merged)
            elif merged.geom_type == 'MultiLineString':
                for line in merged.geoms:
                    network.append(line)

            nodes, links = create_graph_from_face(polygon, network)

            # If there is not 3 or 4 nodes, this is not a branching crossroads
            if len(nodes) not in [3, 4]:
                return False, infos
            
            # Case of a four nodes crossroads
            if len(nodes) == 4:
                # Check if middle node is allowed
                if allow_middle_node == False:
                    return False, infos

                middle = None
                # Looping through each node
                for i, node in enumerate(nodes):
                    degree = node[0]
                    # If the degree is different than 3, return False
                    if degree != 3:
                        return False, infos
                    # If the angle has been calculated
                    if 0 <= 2 < len(node):
                        angle = node[2]
                        # If the angle is 180, it is the middle node
                        if angle == 180:
                            middle = i
                    else:
                        # If an angle is missing, it means the polygon has a weird shape
                        return False, infos
                # If there is no middle node, this is not a branching crossroad
                if middle is None:
                    return False, infos
                else:
                    is_triangle, dist_area = is_triangular(maximum_distance_area, polygon, nodes, middle=middle)
                    if is_triangle == False:
                        return False, infos
                    else:
                        infos['distance_area'] = dist_area
                
                # Here, the 4 nodes crossroad is considered a branching crossroad
                infos['type'] = '4 nodes'
                return True, infos

            # Case of a 3 nodes crossroad
            if len(nodes) == 3:
                nb_four = 0
                for node in nodes:
                    degree = node[0]
                    if degree not in [3, 4]:
                        return False, infos
                    if degree == 4:
                        nb_four += 1
                
                if nb_four > 0:
                    if allow_single_4degree_node == False:
                        return False, infos

                    # If more than one node has a degree of 4, return False
                    if nb_four > 1:
                        return False, infos

                    if nb_four == 1:
                        infos['type'] = '3 nodes single 4 degree'
                
                if roundabouts is not None:
                    roundabouts = roundabouts.to_dict('records')

                    rlist = []
                    for roundabout in roundabouts:
                        rlist.append(roundabout['geometry'])

                    rtree = shapely.STRtree(rlist)
                    # Retrieve roundabouts that intersects the considered network face using strtree
                    rintersects = rtree.query(polygon)
                    rb = []
                    for i in rintersects:
                        r = roundabouts[i]['geometry']
                        ri = roundabouts[i]['cid']
                        # Make an other test to really keep only intersecting roundabouts, spatial index strtree using bbox
                        if shapely.intersects(polygon, r):
                            # Make sure, it is not the same object
                            if shapely.equals_exact(polygon, r):
                                continue
                            else:
                                rb.append([ri, r])
                            
                    if len(rb) != 0:                        
                        # If more than one roundabout intersects, return False
                        if len(rb) > 1:
                            return False, infos

                        # If the area of the face is twice larger than the area of the roundabout, return False
                        if area > (2 * rb[0][1].area):
                            return False, infos
                        else:
                            infos['roundabout'] = rb[0][0]

                is_triangle, dist_area = is_triangular(maximum_distance_area, polygon, nodes)
                if is_triangle == False:
                    return False, infos
                else:
                    infos['distance_area'] = dist_area

                if 'type' not in infos.keys():
                    infos['type'] = '3 nodes'

                # There, this is a branching crossroad
                return True, infos

    return False, infos


def detect_branching_crossroads(roads, area_threshold=2500,
        maximum_distance_area=0.5, roundabouts=None,
        allow_middle_node=False, allow_single_4degree_node=True
    ):
    """
    This function detects brandching crossroads inside a road network and returns polygons representing their extents.
    Parameters
    ----------
    roads : geopandas GeoDataFrame of LineStrings.
        The road network to analyze.
    area_threshold : int optional.
        The area (in square meters) above which the object is not considered a branching crossroads.
        The default value is set to 2500.
    maximum_distance_area : float optional.
        The maximum distance area between the actual polygon and the triangle formed by the 3 nodes connecting the junction to the rest of the network.
        The default value is set to 0.5.
    roundabouts : geopandas GeoDataFrame of Polygons optional
        The polygons representing the network faces considered as roundabouts.
        If provided, it offers a better detection of branching crossroads.
        The default value is set to None.
    allow_middle_node : boolean optional
        If set to True, allow 4 nodes to form the crossroads, but each must have a degree of 3 and the 'middle' node must have an angle of 180°.
        Default value set to False.
    allow_single_4degree_node : boolean optional
        If set to True, allow one and only one node to have a degree of 4.
        Default value set to False.
    Cannot have both the allow_middle_node and allow_single_4degree_node parameters set to True.
    """

    network = []
    for road in roads.geometry:
        network.append(road)

    faces = calculate_network_faces(network, convex_hull=False)

    crossroads = []
    for i, face in enumerate(faces):
        add, infos = is_branching_crossroad(
            face, network, area_threshold,
            maximum_distance_area, roundabouts,
            allow_middle_node, allow_single_4degree_node
        )
        if add:
            infos['cid'] = i
            crossroads.append(infos)

    if len(crossroads) > 0:
        return gpd.GeoDataFrame(crossroads, crs=roads.crs)
    else:
        return None