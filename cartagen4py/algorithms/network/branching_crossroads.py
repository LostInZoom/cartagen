import geopandas as gpd
import shapely
from shapely.ops import unary_union, linemerge

from cartagen4py.utils.network import *
from cartagen4py.utils.geometry.line import *

def collapse_branching_crossroads(roads, crossroads, roundabouts=None, maximum_area=None,):
    """
    Collapse detected branching crossroads below the provided area to a point and return the new network.
    Parameters
    ----------
    roads : geopandas GeoDataFrame of LineStrings.
        The road network where branching crossroads will be collapsed.
    crossroads : geopandas GeoDataFrame of Polygons.
        The polygons representing the faces of the network detected as branching crossroads.
    maximum_area : float, optional.
        The area, in square meter, below which branching crossroads are collapsed.
        Default value is set to None. 
    roundabouts : geopandas GeoDataFrame of Polygons, optional.
        The polygons representing the faces of the network detected as roundabouts.
        Provide a better collapsing when provided.
        Default value is set to None.
    """

    # Retrieve crs for output
    crs = roads.crs

    # Convert geodataframe to list of dicts
    roads = roads.to_dict('records')
    crossroads = crossroads.to_dict('records')
    if roundabouts is not None:
        roundabouts = roundabouts.to_dict('records')
    
    # Store road geometries
    network = []
    for road in roads:
        network.append(road['geometry'])

    # Build the spatial index
    tree = shapely.STRtree(network)

    # This list will store the indexes of roads to throw away
    originals = []
    # This list will store the geometries of internal and external roads.
    # It is important when crossroads are too close to each other and external and internal roads might overlap.
    internals = []
    externals = []
    # This will store the new collapsed geometries
    collapsed = []

    # Looping through crossroads
    for index, crossroad in enumerate(crossroads):
        middle = crossroad['middle']

        # If it's connected to a roundabout, don't collapse
        # The collapsing will take place during the roundabout collapsing function
        if crossroad['roundabout'] == -1:
            # Retrieve geometry, linear ring and area of the crossroad
            polygon = crossroad['geometry']
            linear_ring = polygon.exterior
            area = polygon.area
            
            # Check if area is below the given parameter before proceeding
            if maximum_area is not None:
                if area >= maximum_area:
                    continue

            intersection = Crossroad(network, tree, polygon)

            if intersection is not None:
                main_roads = __find_main_roads(intersection, middle)

                # Test if main road(s) were found
                if main_roads is not None:
                    # Finding the minor external road
                    minor_road = __find_minor_road(intersection, main_roads)
                    # Test if the minor road was found
                    if minor_road is not None:
                        c = __collapse(intersection, main_roads, minor_road)
                        i = []
                        for internal in intersection.internals:
                            i.append(intersection.network[internal])
                        e = []
                        for external in intersection.externals:
                            e.append(intersection.network[external])
                        originals.extend(intersection.original)
                        internals.extend(i)
                        externals.extend(e)
                        collapsed.extend(c)
    
    result = []
    for rid, road in enumerate(roads):
        if rid not in originals:
            result.append(road)

    alreadyc = []
    # Add the collapsed road sections
    for c in collapsed:
        add = True
        # Check if the collapsed section hasn't been added yet
        for ac in alreadyc:
            if shapely.equals(ac, c):
                add = False
        if add:
            result.append({'geometry': c})
            alreadyc.append(c)

    alreadye = []
    # Now add all the external road sections
    for e in externals:
        add = True
        # Check if this external road section has not been added yet
        for ae in alreadye:
            if shapely.equals(ae, e):
                add = False
        # Check conflicts with other internal road sections
        for i in internals:
            if shapely.intersects(i, e):
                # If the external road equals the internal road section of an other crossroad
                if shapely.equals(i, e):
                    add = False
                # If the external road overlaps an other internal road
                elif shapely.overlaps(i, e):
                    add = False
                # If it is contained inside an other internal road section
                elif shapely.contains(i, e):
                    add = False
        if add:
            result.append({'geometry': e})
            alreadye.append(e)
                
    return gpd.GeoDataFrame(result, crs=crs)


def __collapse(branching, mains, minor):
    """
    Collapse the main road(s) and the minor road to a point.
    Returns the indexes of roads to get rid of, and the collapsed crossroad.
    TODO: Propagate the attributes when collapsing the minor road on to the main one.
    """

    # Create the new minor road extension inside the branching crossroad
    def create_minor_extension(point, main_road):
        # Calculate the distance between the starting point of the main road and the projection of the point on the road.
        dist = main_road.project(point)
        # Create the point projected on the main road.
        intersection = shapely.Point(list(main_road.interpolate(dist).coords))

        new_minor = None
        new_main = []

        # Case if the projected point if already the start or end of the main road linestring
        if dist == 0 or dist == main_road.length:
            # The main road is the one provided, no changes required
            new_main.append(main_road)
            # The minor extension is the line between the provided point and the intersection point
            new_minor = shapely.LineString([point, intersection])
        # Case if the projected point is along the main road
        else:
            # Split the main road at the given point
            main1, main2 = split_line_at_point(main_road, point)
            new_main.extend([main1, main2])
            new_minor = shapely.LineString([point, main2.coords[0]])

        return new_minor, new_main

    oroad = None
    # If the are two main roads
    if len(mains) > 1:
        distance = None
        # Loop through main roads to keep only the closest to the minor road
        # Stores in oroad the index of the other main road to keep
        for index, main in enumerate(mains):
            dist = shapely.distance(branching.network[minor], branching.network[main])
            if distance is not None:
                if dist > distance:
                    oroad = index
                    continue
                else:
                    oroad = index - 1
            distance = dist
            main_road = main
    else:
        main_road = mains[0]
    
    # Retrieve start and end point of the minor road
    start = shapely.Point(branching.network[minor].coords[0])
    end = shapely.Point(branching.network[minor].coords[-1])

    # Retrieve the geometry of the concerned main road
    main_geom = branching.network[main_road]

    startdist = shapely.distance(start, main_geom)
    enddist = shapely.distance(end, main_geom)
    # Check if start or end point is the one to be extended
    # There it's the start of the line
    if startdist < enddist:
        new_minor, new_main = create_minor_extension(start, main_geom)
    # There it's the end point
    else:
        new_minor, new_main = create_minor_extension(end, main_geom)

    # If oroad is not None, it means that there are two main roads (with a middle node)
    if oroad is not None:
        new_main.append(branching.network[mains[oroad]])

    collapsed = [new_minor]

    for nm in new_main:
        collapsed.append(nm)

    return collapsed


def __find_main_roads(branching, middle):
    """
    Find the main interior road(s) of the given branching crossroad.
    """
    # Case of a four nodes crossroads, it's the road containing the middle node
    if len(branching.nodes) == 4:
        # Safety, but 4 nodes crossroads should have a middle node
        if middle != -1:
            # Retrieve coordinates
            mcoords = branching.nodes[middle][1]
            main = []
            # Loop through internal roads
            for internal in branching.internals:
                coords = branching.network[internal].coords
                start, end = coords[0], coords[-1]
                # Check if start or end point intersects with the middle node
                if start == mcoords or end == mcoords:
                    main.append(internal)
            return main
    
    # Case of a three nodes crossroads, it's the internal road with two flat angles with external roads
    if len(branching.nodes) == 3:
        main = []
        # Loop through internal roads
        for internal in branching.internals:
            # Set the number of flat angle with external road to 0
            flat = 0
            # Loop through external roads
            for external in branching.externals:
                # If it touches the internal road
                if shapely.touches(branching.network[internal], branching.network[external]):
                    # Calculate the angle between the two lines
                    angle = abs(angle_between_2lines(branching.network[internal], branching.network[external]))
                    if angle > (np.pi / 2):
                        angle = abs(np.pi - angle)
                    if angle < (np.pi / 8):
                        flat += 1
            
            if flat == 2:
                main.append(internal)

        if len(main) == 1:
            return main
            
    return None


def __find_minor_road(branching, main_roads):
    """
    Find the minor road, i.e. the road usually perpendicular to the main road.
    It it the only external road that doesn't intersect in any way the main roads(s).
    """
    minor = []

    # If there are more than one main road, union all the lines
    if len(main_roads) > 1:
        fullmain = []
        for main in main_roads:
            fullmain.extend(branching.network[main].coords)
        main_road = shapely.LineString(fullmain)
    else:
        main_road = branching.network[main_roads[0]]

    for external in branching.externals:
        if shapely.intersects(main_road, branching.network[external]) == False:
            minor.append(external)
    
    if len(minor) == 1:
        return minor[0]
    else:
        return None
