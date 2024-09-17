import shapely
import geopandas as gpd
import numpy as np
from cartagen.utils.partitioning.network import network_faces
from cartagen.utils.network.roads import Crossroad
from cartagen.utils.geometry.distances import distance_area

def detect_branching_crossroads(roads, roundabouts=None,
        area_threshold=2500, maximum_distance_area=0.5, 
        allow_middle_node=True, middle_angle_tolerance=10.0,
        allow_single_4degree_node=True
    ):
    """
    Detect branching crossroads based on geometric properties.

    This algorithm proposed by Touya :footcite:p:`touya:2010` detects
    branching crossroads inside a road network based on the proximity between
    the geometry of the network face and a triangle.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        Road network to analyze.
    roundabouts : GeoDataFrame of Polygon, optional
        The polygons representing the network faces considered as roundabouts.
        If provided, links the branching crossroad to a roundabout for collapsing.
    area_threshold : int, optional
        The area (in square meters) above which the object is not considered a branching crossroads.
    maximum_distance_area : float, optional
        The maximum distance area between the actual polygon
        and the triangle formed by the 3 nodes connecting
        the junction to the rest of the network.
    allow_middle_node : bool, optional
        If set to True, allow 4 nodes to form the crossroads,
        but each must have a degree of 3 and the 'middle'
        node must have an angle of 180°.
    middle_angle_tolerance : float, optional
        If allow_middle_node is set to True,
        indicate an angle tolerance in degree
        for the fourth node of the crossroad to be considered the middle node.
    allow_single_4degree_node : bool, optional
        If set to True, allow one and only one node to have a degree of 4.

    Returns
    -------
    GeoDataFrame of Polygon

    Warning
    -------
    Detecting roundabouts beforehand is important as a branching crossroad
    may be an entrance to a roundabout. This algorithm will link branching
    crossroads to a roundabout when applicable, and this will help collapsing
    both objects.

    See Also
    --------
    detect_roundabouts : 
        Detect roundabouts inside the road network.
    collapse_roundabouts :
        Collapse roundabouts to a point.
    collapse_branching_crossroads :
        Collapse branching crossroads to a point.
    
    References
    ----------
    .. footbibliography::
    """

    crs = roads.crs
    roads = roads.to_dict('records')

    network = []
    for road in roads:
        network.append(road['geometry'])

    faces = network_faces(network, convex_hull=False)
    tree = shapely.STRtree(network)

    crossroads = []
    index = 0
    for face in faces:
        add, infos = is_branching_crossroad(
            face, network, tree, area_threshold,
            maximum_distance_area, roundabouts,
            allow_middle_node, middle_angle_tolerance,
            allow_single_4degree_node
        )
        if add:
            infos['cid'] = index
            crossroads.append(infos)
            index += 1

    if len(crossroads) > 0:
        return gpd.GeoDataFrame(crossroads, crs=crs)
    else:
        return gpd.GeoDataFrame()


def is_branching_crossroad(polygon, roads, tree,
        area_threshold, maximum_distance_area, roundabouts=None,
        allow_middle_node=True, middle_angle_tolerance=10,
        allow_single_4degree_node=True
    ):
    """
    Detect if the provided network face is a branching crossroad.

    Return True or False whether the given polygon is a
    branching crossroad or not depending on the given parameters.

    Parameters
    ----------
    polygon : Polygon
        The geometry of the network face to check.
    roads : GeoDataFrame of LineString
        The road network to analyze.
    tree : STRTree
        A shapely spatial index calculated on all provided roads.
    area_threshold : int, optional
        The area (in square meters) above which the object is not considered a branching crossroads.
        The default value is set to 2500.
    maximum_distance_area : float, optional
        The maximum distance area between the actual polygon and the triangle formed by the 3 nodes connecting the junction to the rest of the network.
        The default value is set to 0.5.
    roundabouts : GeoDataFrame of Polygon
        The polygons representing the network faces considered as roundabouts.
        If provided, it offers a better detection of branching crossroads.
        The default value is set to None.
    allow_middle_node : bool, optional
        If set to True, allow 4 nodes to form the crossroads, but each must have a degree of 3 and the 'middle' node must have an angle of 180°.
        Default value set to False.
    allow_single_4degree_node : bool, optional
        If set to True, allow one and only one node to have a degree of 4.
        Default value set to False.
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

    infos = {
        'geometry': polygon,
        'roundabout': -1,
        'middle': -1
    }

    area = polygon.area
    # Check if the area of the polygon is larger than the threshold
    if area < area_threshold:
        crossroad = Crossroad(roads, tree, polygon)

        if crossroad is not None:
            nodes = crossroad.nodes
            links = crossroad.links

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
                        # If the angle is below pi/8, it is considered flat
                        if angle < (np.pi / 8):
                            middle = i
                    else:
                        # If an angle is missing, it means the polygon has a weird shape
                        return False, infos
                # If there is no middle node, this is not a branching crossroad
                if middle is None:
                    return False, infos
                else:
                    # Test if the shape is triangular
                    is_triangle, dist_area = is_triangular(maximum_distance_area, polygon, nodes, middle=middle)
                    if is_triangle == False:
                        return False, infos
                    else:
                        # Retrieve the main road id
                        # Add the id of the middle node and the distance area to the dict
                        infos['middle'] = middle
                        infos['distance_area'] = dist_area
                
                if roundabouts is not None:
                    proceed, rb = __find_roundabout_index(roundabouts, polygon, area)
                    if proceed == False:
                        return False, infos
                    else:
                        if rb is not None:
                            infos['roundabout'] = rb

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
                    proceed, rb = __find_roundabout_index(roundabouts, polygon, area)
                    if proceed == False:
                        return False, infos
                    else:
                        if rb is not None:
                            infos['roundabout'] = rb

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

def __find_roundabout_index(roundabouts, face, area):
    """
    Return a list of indexes of intersecting roundabouts
    """
    roundabouts = roundabouts.to_dict('records')
    rlist = []
    for roundabout in roundabouts:
        rlist.append(roundabout['geometry'])

    rtree = shapely.STRtree(rlist)
    # Retrieve roundabouts that intersects the considered network face using strtree
    rintersects = rtree.query(face)
    rb = []
    for i in rintersects:
        r = roundabouts[i]['geometry']
        ri = roundabouts[i]['cid']
        # Make an other test to really keep only intersecting roundabouts, spatial index strtree using bbox
        if shapely.intersects(face, r):
            # Make sure, it is not the same object
            if shapely.equals_exact(face, r):
                continue
            else:
                rb.append([ri, r])

    if len(rb) != 0:
        # If more than one roundabout intersects, return None
        if len(rb) > 1:
            return False, None

        # If the area of the face is twice larger than the area of the roundabout, return None
        if area > (2 * rb[0][1].area):
            return False, None
        else:
            return True, rb[0][0]
    else:
        return True, None