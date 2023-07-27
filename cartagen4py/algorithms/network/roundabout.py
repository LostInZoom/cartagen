import geopandas as gpd
import shapely

from cartagen4py.utils.geometry.line import extend_line_with_point

def collapse_roundabouts(network, roundabouts, maximum_diameter=None):
    """
    Collapse detected roundabouts below the provided diameter to a point and return the new network.
    Parameters
    ----------
    network : geopandas GeoDataFrame of LineStrings.
        The road network where roundabouts will be collapsed.
    roundabouts : geopandas GeoDataFrame of Polygons.
        The polygons representing the faces of the network detected as roundabouts.
    maximum_diameter : float, optional.
        The diameter, in meter, below which roundabouts are collapsed.
        Default value is set to None. 
    """

    # Retrieve crs for output
    crs = network.crs

    # Convert geodataframe to list of dicts
    roundabouts = roundabouts.to_dict('records')
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
    for roundabout in roundabouts:
        # Retrieve geometry, linear ring and centroid of the roundabout
        geometry = roundabout['geometry']
        linear_ring = geometry.exterior
        centroid = geometry.centroid
        # Calculate the diameter of the roundabout, which is the closest distance between its centroid and border
        diameter = shapely.distance(centroid, linear_ring)

        # Check if diameter is below the given parameter before proceeding
        if maximum_diameter is not None:
            if diameter >= maximum_diameter:
                continue

        # Selecting intersecting roads index
        intersecting = rtree.query(geometry)

        # Future storage of internal (the roundabout roads themselves)
        # and external (the entrances/exits of the roundabout)
        internal = []
        external = []

        # Looping through intersecting roads index
        for i in intersecting:
            # Get the road geometry
            r = roads[i]

            # If the road is contained within the linear ring of the roundabout, it is interior
            if shapely.contains(linear_ring, r):
                internal.append(i)
            # If the road is only touching the linear ring of the roundabout, it is exterior
            elif shapely.touches(linear_ring, r):
                external.append(i)

        # Loop through internal roads
        for i in internal:
            # Add the internal road for removal
            if i not in nokeep:
                nokeep.append(i)

        # Future storage for multiple external roads entering the roundabout at the same point
        noextend = []
        # Loop through external roads index
        for e in external:
            # Retrieve external road geometry
            ext_road = roads[e]
            # Retrieve coordinates of start and end of the road line
            start, end = ext_road.coords[0], ext_road.coords[-1]
            # Creating points for the start and end
            first = shapely.Point(start)
            last = shapely.Point(end)

            # If the first point is closer to the centroid of the roundabout than the last
            if shapely.distance(first, centroid) < shapely.distance(last, centroid):
                # If the point as not already been extended to the centroid
                if start not in noextend:
                    # Extend the line from the start to the centroid
                    network[e]['geometry'] = extend_line_with_point(ext_road, centroid, position='start')
                    # Add the point if an other road passes through it to avoid overlapping
                    noextend.append(start)
            # Else, do the same for the last point
            else:
                if end not in noextend:
                    network[e]['geometry'] = extend_line_with_point(ext_road, centroid, position='end')
                    noextend.append(end)

    # Storing the collapsed network into a new list
    collapsed = []
    for i, n in enumerate(network):
        # Check if the index is to be kept
        if i not in nokeep:
            collapsed.append(n)
    
    return gpd.GeoDataFrame(collapsed, crs=crs)