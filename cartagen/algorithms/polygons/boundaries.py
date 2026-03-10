from cartagen.utils.geometry.topo_map import TopoMap
from cartagen.utils.lines import visvalingam_whyatt, douglas_peucker, raposo, li_openshaw
from cartagen.utils.geometry.polygon import surfacic_distance

def boundaries_visvalingam_whyatt(polygons, threshold):
    """
    Applies the Visvalingam-Whyatt :footcite:p:`visvalingam:1993` algorithm to the boundaries of the given polygons.
    As most polygons share a boundary with another polygon, the simplification is only applied to the common line,
    so that no topological disconnection is created between adjacent polygons.

    Parameters
    ----------
    polygons : GeoDataFrame of Polygons
        Polygons forming the set of boundaries to simplify.
    threshold : float
        The minimum triangle area to keep a vertex in the line.
        Higher values = more points kept (less aggressive simplification).

    Returns
    -------
    GeoDataFrame of Polygons

    Warning
    -------
    This algorithm cannot create MultiPolygon, unlike :func:`hull_swinging_arm`.
    Using a length too low can produce an invalid geometry.

    See Also
    --------
    visvalingam_whyatt :
        Area-based line simplification.
    boundaries_raposo :
        Apply the Raposo line simplification to the boundaries of polygons.
    boundaries_douglas_peucker :
        Apply the Douglas-Peucker line simplification to the boundaries of polygons.
    boundaries_li_openshaw :
        Apply the Li-Openshaw line simplification to the boundaries of polygons.

    References
    ----------
    .. footbibliography::
    """

    polygons_simplified = polygons.copy()
    polygons_to_search = []

    for index, feature in polygons.iterrows():
        polygons_to_search.append(feature['geometry'])

    # first create a TopoMap
    tm = TopoMap(snap_tolerance=1e-9, build_infinite_face=True)
    tm.build_from_polygons(polygons_to_search)


    # then loop on all arcs of the TopoMap
    for aid, arc in tm.arcs.items():
        geom = arc.geom
        simplified = visvalingam_whyatt(geom, threshold=threshold)
        tm.arcs[aid].geom = simplified

    # update the TopoMap with the simplified arcs
    tm.clear_faces()
    tm.build_faces()# we need to rebuild the faces from scratch, as some arcs may have been removed

    # then match the simplified faces of the TopoMap with the initial polygons
    for face in tm.faces:
        poly = tm.faces[face].polygon
        # check if it is the infinite face
        if poly == None:
            continue
        intersection = polygons_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = polygons_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            polygons_simplified.at[match, 'geometry'] = poly
    
    return polygons_simplified


def boundaries_douglas_peucker(polygons, threshold, preserve_topology=True):
    """
    Applies the Ramer-Douglas-Peucker :footcite:p:`ramer:1972` :footcite:p:`douglas:1973` algorithm to the boundaries of polygons.
    As most polygons share a boundary with another polygon, the simplification is only applied to the common line,
    so that no topological disconnection is created between adjacent polygons.

    Parameters
    ----------
    polygons : GeoDataFrame of Polygons
        The polygon forming the set of boundaries to simplifiy.
    threshold : float
        The distance threshold to remove the vertex from the line.
    preserve_topology : bool, optional
        If set to True, the algorithm will prevent invalid geometries
        from being created (checking for collapses, ring-intersections, etc).
        The trade-off is computational expensivity.

    Returns
    -------
    GeoDataFrame of Polygons

    Warning
    -------
    This algorithm cannot create MultiPolygon, so having input MultiPolygon can cause matching problems.

    See Also
    --------
    douglas_peucker :
        Distance-based line simplification.
    boundaries_visvalingam :
        Apply the Visvalingam-Whyatt line simplification to the boundaries of polygons.
    boundaries_raposo :
        Apply the Raposo line simplification to the boundaries of polygons.
    boundaries_li_openshaw :
        Apply the Li-Openshaw line simplification to the boundaries of polygons.

    References
    ----------
    .. footbibliography::
    """

    polygons_simplified = polygons.copy()
    polygons_to_search = []

    for index, feature in polygons.iterrows():
        polygons_to_search.append(feature['geometry'])

    # first create a TopoMap
    tm = TopoMap(snap_tolerance=1e-9, build_infinite_face=True)
    tm.build_from_polygons(polygons_to_search)


    # then loop on all arcs of the TopoMap
    for aid, arc in tm.arcs.items():
        geom = arc.geom
        simplified = douglas_peucker(geom, threshold, preserve_topology)
        tm.arcs[aid].geom = simplified

    # update the TopoMap with the simplified arcs
    tm.clear_faces()
    tm.build_faces()# we need to rebuild the faces from scratch, as some arcs may have been removed

    # then match the simplified faces of the TopoMap with the initial polygons
    for face in tm.faces:
        poly = tm.faces[face].polygon
        # check if it is the infinite face
        if poly == None:
            continue
        intersection = polygons_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = polygons_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            polygons_simplified.at[match, 'geometry'] = poly
    
    return polygons_simplified



def boundaries_raposo(polygons, initial_scale, final_scale, centroid=True, tobler=False):
    """
    Applies the Raposo :footcite:p:`raposo:2013` algorithm to the boundaries of polygons.
    As most polygons share their boundaries with another polygon, the simplification is only applied to the common line,
    so that no topological disconnection is created between adjacent polygons.

    Parameters
    ----------
    polygons : GeoDataFrame of Polygon
        The polygon forming the set of boundaries to simplifiy.
    initial_scale : float
        Initial scale of the provided line (25000.0 for 1:25000 scale).
    final_scale : float
        Final scale of the simplified line.
    centroid : bool, optional
        If true, uses the center of the hexagonal cells as the new vertex,
        if false, the center is projected on the nearest point in the initial line.
    tobler : bool, optional
        If True, compute cell resolution based on Tobler’s formula, else uses Raposo's formula

    Returns
    -------
    GeoDataFrame of Polygon

    Warning
    -------
    This algorithm cannot create MultiPolygon, so having input MultiPolygon can cause matching problems.

    See Also
    --------
    raposo :
        Hexagon-based line simplification.
    boundaries_douglas_peucker :
        Apply the Douglas-Peucker line simplification to the boundaries of polygons.
    boundaries_visvalingam :
        Apply the Visvalingam-Whyatt line simplification to the boundaries of polygons.
    boundaries_li_openshaw :
        Apply the Li-Openshaw line simplification to the boundaries of polygons.

    References
    ----------
    .. footbibliography::
    """

    polygons_simplified = polygons.copy()
    polygons_to_search = []

    for index, feature in polygons.iterrows():
        polygons_to_search.append(feature['geometry'])

    # first create a TopoMap
    tm = TopoMap(snap_tolerance=1e-9, build_infinite_face=True)
    tm.build_from_polygons(polygons_to_search)


    # then loop on all arcs of the TopoMap
    for aid, arc in tm.arcs.items():
        geom = arc.geom
        simplified = raposo(geom, initial_scale, final_scale, centroid, tobler)
        tm.arcs[aid].geom = simplified

    # update the TopoMap with the simplified arcs
    tm.clear_faces()
    tm.build_faces()# we need to rebuild the faces from scratch, as some arcs may have been removed

    # then match the simplified faces of the TopoMap with the initial polygons
    for face in tm.faces:
        poly = tm.faces[face].polygon
        # check if it is the infinite face
        if poly == None:
            continue
        intersection = polygons_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = polygons_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            polygons_simplified.at[match, 'geometry'] = poly
    
    return polygons_simplified


def boundaries_li_openshaw(polygons, cell_size):
    """
    Applies the Li-Openshaw :footcite:p:`li:1993` algorithm to the boundaries of the polygons. As most polygons share a boundary
    with another polygon, the simplification is only applied to the common line, so that no topological disconnection is
    created between adjacent polygons.

    Parameters
    ----------
    polygons : GeoDataFrame of Polygon
        The polygon forming the set of boundaries to simplifiy.
    cell_size : float
        The size of the regular grid used to divide the line.

    Returns
    -------
    GeoDataFrame of Polygon

    Warning
    -------
    This algorithm cannot create MultiPolygon, so having input MultiPolygon can cause matching problems.

    See Also
    --------
    li_openshaw :
        Square grid-based line simplification.
    boundaries_douglas_peucker :
        Applies the Douglas-Peucker line simplification to the boundaries of polygons.
    boundaries_raposo :
        Applies the Raposo line simplification to the boundaries of polygons.
    boundaries_visvalingam :
        Applies the Visvalingam-Whyatt line simplification to the boundaries of polygons.

    References
    ----------
    .. footbibliography::
    """

    polygons_simplified = polygons.copy()
    polygons_to_search = []

    for index, feature in polygons.iterrows():
        polygons_to_search.append(feature['geometry'])

    # first create a TopoMap
    tm = TopoMap(snap_tolerance=1e-9, build_infinite_face=True)
    tm.build_from_polygons(polygons_to_search)


    # then loop on all arcs of the TopoMap
    for aid, arc in tm.arcs.items():
        geom = arc.geom
        simplified = li_openshaw(geom, cell_size)
        tm.arcs[aid].geom = simplified

    # update the TopoMap with the simplified arcs
    tm.clear_faces()
    tm.build_faces()# we need to rebuild the faces from scratch, as some arcs may have been removed

    # then match the simplified faces of the TopoMap with the initial polygons
    for face in tm.faces:
        poly = tm.faces[face].polygon
        # check if it is the infinite face
        if poly == None:
            continue
        intersection = polygons_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = polygons_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            polygons_simplified.at[match, 'geometry'] = poly
    
    return polygons_simplified