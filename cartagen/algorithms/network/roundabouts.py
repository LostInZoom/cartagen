import geopandas as gpd
import shapely

from cartagen.utils.geometry.line import extend_line_with_point

def collapse_roundabouts(roads, roundabouts, crossroads=None, maximum_diameter=None):
    """
    Collapse roundabouts to a point.

    This algorithm proposed by Touya :footcite:p:`touya:2010` collapses roundabouts to a point
    if their diameter is below the given value.
    
    Parameters
    ----------
    roads : GeoDataFrame of LineString
        Road network where roundabouts will be collapsed.
    roundabouts : GeoDataFrame of Polygon
        Polygons representing the faces of the network detected as roundabouts.
    crossroads : GeoDataFrame of Polygon, optional
        Polygons representing the faces of the network detected as branching crossroads. This
        allows incoming branching crossroads on roundabouts to be collapsed as well. 
    maximum_diameter : float, optional
        Diameter, in meter, below which roundabouts are collapsed. Collpase all roundabouts if left to None.

    Returns
    -------
    GeoDataFrame of LineString

    Warning
    -------
    Detecting branching crossroads beforehand is important as a branching crossroad
    may be an entrance to a roundabout. This algorithm will collapse the roundabout
    as well as all its connected branching crossroads.

    See Also
    --------
    detect_roundabouts : 
        Detect roundabouts inside the road network.
    detect_branching_crossroads : 
        Detect branching crossroads inside the road network.
    collapse_branching_crossroads :
        Collapse branching crossroads to a point.

    References
    ----------
    .. footbibliography::
    """
    # Retrieve crs for output
    crs = roads.crs

    # Convert geodataframe to list of dicts
    roundabouts = roundabouts.to_dict('records')

    if len(roundabouts) == 0:
        return roads

    roads = roads.to_dict('records')
    if crossroads is not None:
        crossroads = crossroads.to_dict('records')

    # Create a list of all the roads geometry of the network
    network = []
    for n in roads:
        network.append(n['geometry'])

    # Calculate the spatial index on roads
    tree = shapely.STRtree(network)

    # This list will store the indexes of interior roads to throw away
    nokeep = []

    # Looping through roundabouts
    for roundabout in roundabouts:
        # Retrieve geometry and centroid of the roundabout
        polygon = roundabout['geometry']
        centroid = polygon.centroid
        linear_ring = polygon.exterior
        # Calculate the diameter of the roundabout, which is the closest distance between its centroid and border
        diameter = shapely.distance(centroid, linear_ring)

        # Check if diameter is below the given parameter before proceeding
        if maximum_diameter is not None:
            if diameter >= maximum_diameter:
                continue

        connected_crossroads = []
        if crossroads is not None:
            for crossroad in crossroads:
                if crossroad['roundabout'] == roundabout['cid']:
                    connected_crossroads.append(crossroad)

        if len(connected_crossroads) > 0:
            cgeom = []
            for c in connected_crossroads:
                cgeom.append(c['geometry'])
            face = shapely.ops.unary_union([polygon, *cgeom])
        else:
            face = polygon

        # Selecting intersecting roads index
        intersecting = tree.query(face)

        full_ring = face.exterior

        # Future storage of internal (the roundabout roads themselves)
        # and external (the entrances/exits of the roundabout)
        internals = []
        externals = []

        # Looping through intersecting roads index
        for i in intersecting:
            # Get the road geometry
            r = network[i]
            if shapely.intersects(full_ring, r):
                # If the road is contained within the linear ring of the roundabout, it is interior
                if shapely.contains(full_ring, r):
                    internals.append(i)
                # If the road is only touching the linear ring of the roundabout but is not contained within its polygon, it is exterior
                elif shapely.touches(full_ring, r):
                    if shapely.contains(linear_ring, r):
                        # Add the road for removal
                        if i not in nokeep:
                            nokeep.append(i)
                    else:
                        externals.append(i)

        # Loop through internal roads
        for i in internals:
            # Add the internal road for removal
            if i not in nokeep:
                nokeep.append(i)

        # Future storage for multiple external roads entering the roundabout at the same point
        noextend = []
        # Loop through external roads index
        for external in externals:
            # Retrieve external road geometry
            e = roads[external]['geometry']
            # Retrieve coordinates of start and end of the road line
            start, end = e.coords[0], e.coords[-1]
            # Creating points for the start and end
            first = shapely.Point(start)
            last = shapely.Point(end)

            # If the first point is closer to the centroid of the roundabout than the last
            if shapely.distance(first, centroid) < shapely.distance(last, centroid):
                # If the point as not already been extended to the centroid
                if start not in noextend:
                    # Extend the line from the start to the centroid
                    roads[external]['geometry'] = extend_line_with_point(e, centroid, position='start')
                    # Add the point if an other road passes through it to avoid overlapping
                    noextend.append(start)
            # Else, do the same for the last point
            else:
                if end not in noextend:
                    roads[external]['geometry'] = extend_line_with_point(e, centroid, position='end')
                    noextend.append(end)

    # Storing the collapsed network into a new list
    collapsed = []
    for i, n in enumerate(roads):
        # Check if the index is to be kept
        if i not in nokeep:
            collapsed.append(n)

    return gpd.GeoDataFrame(collapsed, crs=crs)