from shapely.ops import shared_paths
from shapely.strtree import STRTree

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

    References
    ----------
    .. footbibliography::

    Examples
    --------
    ...
    """

    polygon_arcs = []
    polygons_to_search = []
    arcs = []
    arc_id = 0
    # first initialise the list of arcs for each polygon
    for index, feature in boundaries.iterrows():
        polygon_arcs.append([index, []])
        polygons_to_search.append(feature['geometry'])

    # then loop on all polygons
    for index, feature in boundaries.iterrows():
        geom = polygons_to_search.pop()
        tree = STRTree(polygons_to_search)
        for neighbour in tree.query(geom,predicate="intersects"):
            onward,backward = shared_paths(geom,shared_paths)

