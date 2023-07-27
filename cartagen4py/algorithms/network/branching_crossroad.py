import geopandas as gpd
import shapely

def collapse_branching_crossroads(network, crossroads, roundabouts=None, maximum_area=None):
    """
    Collapse detected branching crossroads below the provided area to a point and return the new network.
    Parameters
    ----------
    network : geopandas GeoDataFrame of LineStrings.
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
    crs = network.crs

    # Convert geodataframe to list of dicts
    crossroads = crossroads.to_dict('records')
    network = network.to_dict('records')

     # Create a list of all the roads geometry of the network
    roads = []
    for n in network:
        roads.append(n['geometry'])

    # Calculate the spatial index on roads
    rtree = shapely.STRtree(roads)

    # This list will store the indexes of interior roads to throw away
    nokeep = []

    # Looping through roundabouts
    for crossroad in crossroads:
        # Retrieve geometry, linear ring and area of the crossroad
        geometry = crossroad['geometry']
        linear_ring = geometry.exterior
        area = geometry.area
        
        # Check if area is below the given parameter before proceeding
        if maximum_area is not None:
            if area >= maximum_area:
                continue