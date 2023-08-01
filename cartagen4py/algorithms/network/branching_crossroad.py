import geopandas as gpd
import shapely
from shapely.ops import unary_union, linemerge

from cartagen4py.utils.network import *

def collapse_branching_crossroads(roads, crossroads, roundabouts=None, maximum_area=None):
    """
    Collapse detected branching crossroads below the provided area to a point and return the new network.
    Parameters
    ----------
    roads : geopandas GeoDataFrame of LineStrings.
        The road network where branching crossroads will be collapsed.
    crossroads : geopandas GeoDataFrame of Polygons.
        The polygons representing the faces of the network detected as branching crossroads.
    roundabouts : geopandas GeoDataFrame of Polygons, optional.
        The polygons representing the faces of the network detected as roundabouts.
        Provide a better collapsing when provided.
        Default value is set to None.
    maximum_area : float, optional.
        The area, in square meter, below which branching crossroads are collapsed.
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

    # This list will store the indexes of interior minor roads to throw away
    nokeep = []

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

            branching = Crossroad(network, polygon, tree)

            if branching is not None:
                print(crossroad['cid'])
                main_roads = __find_main_roads(branching, middle)

                # Test if main road(s) were found
                if main_roads is not None:
                    # Finding the minor external road
                    minor_road = __find_minor_road(branching, main_roads)
                    # Test if the minor road was found
                    if minor_road is not None:
                        __collapse(branching, main_roads, minor_road)

                        major_export = []
                        for net in main_roads:
                            major_export.append({'geometry': branching.network[net]})

                        major_export_gdf = gpd.GeoDataFrame(major_export, crs=crs)
                        major_export_gdf.to_file("cartagen4py/data/mains.geojson", driver="GeoJSON")

                        minor_export = [{'geometry': branching.network[minor_road]}]
                        minor_export_gdf = gpd.GeoDataFrame(minor_export, crs=crs)
                        minor_export_gdf.to_file("cartagen4py/data/minor.geojson", driver="GeoJSON")

                        if len(main_roads) > 1:
                            break

                # if crossroad["cid"] == 50:
                #     break


def __collapse(branching, mains, minor):
    """
    Collapse the main road(s) and the minor road to a point.
    """
    if len(mains) > 1:
        unioned = []
        for main in mains:
            unioned.extend(branching.network[main].coords)
        main_road = shapely.LineString(unioned)
    else:
        main_road = branching.network[mains[0]]

    start = shapely.Point(branching.network[minor].coords[0])
    end = shapely.Point(branching.network[minor].coords[-1])

    if shapely.distance(start, main_road) < shapely.distance(end, main_road):
        dist = main_road.project(start)
        intersection = list(main_road.interpolate(dist).coords)
    else:
        dist = main_road.project(start)
        intersection = list(main_road.interpolate(dist).coords)

    print(dist, main_road.length, intersection)


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
