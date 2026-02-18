from cartagen.utils.geometry.topo_map import TopoMap
from cartagen.utils.geometry.line import visvalingam_whyatt, douglas_peucker, raposo, li_openshaw
from cartagen.utils.geometry.polygon import surfacic_distance

def boundaries_visvalingam(boundaries, area_tolerance):
    """
    Applies the Visvialingam-Whyatt algorithm to the boundaries of the polygons. As most polygons share their boundaries 
    with another polygon, the simplification is only applied to the common line, so that no topological disconnection is
    created between adjacent polygons.

    The algorithm proposed by Visvalingam and Whyatt :footcite:p:`visvalingam:1993` is used to simplify the boundaries.

    Parameters
    ----------
    boundaries : geodataframe of Polygons
        The polygon forming the set of boundaries to simplifiy.
    area_tolerance : float
        The minimum triangle area to keep a vertex in the line.

    Returns
    -------
    geodataframe of Polygons

    Warning
    -------
    This algorithm cannot create multiple polygons, unlike :func:`hull_swinging_arm`.
    Using a length too low can produce an invalid geometry.

    See Also
    --------
    boundaries_raposo :
        Applies the Raposo line simplification to the boundaries of the polygons.
    boundaries_douglas_peucker :
        Applies the Douglas-Peucker line simplification to the boundaries of the polygons.
    boundaries_li_openshaw :
        Applies the Li-Openshaw line simplification to the boundaries of the polygons.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    ...
    """

    boundaries_simplified = boundaries.copy()
    polygons_to_search = []

    for index, feature in boundaries.iterrows():
        polygons_to_search.append(feature['geometry'])

    # first create a TopoMap
    tm = TopoMap(snap_tolerance=1e-9, build_infinite_face=True)
    tm.build_from_polygons(polygons_to_search)


    # then loop on all arcs of the TopoMap
    for aid, arc in tm.arcs.items():
        geom = arc.geom
        simplified = visvalingam_whyatt(geom, area_tolerance)
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
        intersection = boundaries_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = boundaries_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            boundaries_simplified.at[match, 'geometry'] = poly
    
    return boundaries_simplified


def boundaries_douglas_peucker(boundaries, threshold, preserve_topology=True):
    """
    Applies the Douglas-Peucker-Ramer algorithm to the boundaries of the polygons. As most polygons share their boundaries 
    with another polygon, the simplification is only applied to the common line, so that no topological disconnection is
    created between adjacent polygons.

    The algorithm proposed by Ramer :footcite:p:`ramer:1972` and by Douglas and Peucker :footcite:p:`douglas:1973` is used to simplify the boundaries.

    Parameters
    ----------
    boundaries : geodataframe of Polygons
        The polygon forming the set of boundaries to simplifiy.
    threshold : float
        The distance threshold to remove the vertex from the line.

    Returns
    -------
    geodataframe of Polygons

    Warning
    -------
    This algorithm cannot create multiple polygons, so having input multiple polygons can cause matching problems.

    See Also
    --------
    boundaries_visvalingam :
        Applies the Visvalingam-Whyatt line simplification to the boundaries of the polygons.
    boundaries_raposo :
        Applies the Raposo line simplification to the boundaries of the polygons.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    ...
    """

    boundaries_simplified = boundaries.copy()
    polygons_to_search = []

    for index, feature in boundaries.iterrows():
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
        intersection = boundaries_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = boundaries_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            boundaries_simplified.at[match, 'geometry'] = poly
    
    return boundaries_simplified



def boundaries_raposo(boundaries, initial_scale, final_scale, centroid=True, tobler=False):
    """
    Applies the Raposo algorithm to the boundaries of the polygons. As most polygons share their boundaries 
    with another polygon, the simplification is only applied to the common line, so that no topological disconnection is
    created between adjacent polygons.

    This algorithm proposed by Raposo :footcite:p:`raposo:2013` simplifies lines based on a
    hexagonal tessellation. The algorithm also works for the simplification of the border of a polygon object.
    The idea of the algorithm is to put a hexagonal tessallation on top of the line to simplify,
    the size of the cells depending on the targeted granularity of the line.

    Parameters
    ----------
    boundaries : geodataframe of Polygons
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
    geodataframe of Polygons

    Warning
    -------
    This algorithm cannot create multiple polygons, so having input multiple polygons can cause matching problems.

    See Also
    --------
    boundaries_douglas_peucker :
        Applies the Douglas-Peucker line simplification to the boundaries of the polygons.
    boundaries_visvalingam :
        Applies the Visvalingam-Whyatt line simplification to the boundaries of the polygons.

    References
    ----------
    .. footbibliography::

    """

    boundaries_simplified = boundaries.copy()
    polygons_to_search = []

    for index, feature in boundaries.iterrows():
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
        intersection = boundaries_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = boundaries_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            boundaries_simplified.at[match, 'geometry'] = poly
    
    return boundaries_simplified


def boundaries_li_openshaw(boundaries, cell_size):
    """
    Applies the Li-Openshaw algorithm to the boundaries of the polygons. As most polygons share their boundaries 
    with another polygon, the simplification is only applied to the common line, so that no topological disconnection is
    created between adjacent polygons.

    This algorithm proposed by Li & Openshaw :footcite:p:`li:1993` simplifies lines based on a
    regular square grid. It first divide the line vertexes into groups partionned by a regular
    grid, then each group of vertexes is replaced by their centroid.

    Parameters
    ----------
    boundaries : geodataframe of Polygons
        The polygon forming the set of boundaries to simplifiy.
    cell_size : float
        The size of the regular grid used to divide the line.

    Returns
    -------
    geodataframe of Polygons

    Warning
    -------
    This algorithm cannot create multiple polygons, so having input multiple polygons can cause matching problems.

    See Also
    --------
    boundaries_douglas_peucker :
        Applies the Douglas-Peucker line simplification to the boundaries of the polygons.
    boundaries_raposo :
        Applies the Raposo line simplification to the boundaries of the polygons.
    boundaries_visvalingam :
        Applies the Visvalingam-Whyatt line simplification to the boundaries of the polygons.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    ...
    """

    boundaries_simplified = boundaries.copy()
    polygons_to_search = []

    for index, feature in boundaries.iterrows():
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
        intersection = boundaries_simplified.sindex.query(poly)
        # find the intersecting polygon with the smallest surfacic distance
        dist_min = 1.5
        match = None
        for idx in intersection:
            candidate = boundaries_simplified.iloc[idx].geometry
            dist = surfacic_distance(poly, candidate)
            if dist < dist_min:
                dist_min = dist
                match = idx
        
        if match is not None:
            boundaries_simplified.at[match, 'geometry'] = poly
    
    return boundaries_simplified