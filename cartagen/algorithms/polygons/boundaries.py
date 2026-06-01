from cartagen.utils.geometry.topo_map import TopoMap
from cartagen.utils.geometry.polygon import surfacic_distance

def generalize_boundaries(polygons, algorithm, *args, **kwargs):
    """
    Generalize shared polygon boundaries while strictly preserving topology.

    This function applies a simplification or smoothing algorithm to a set of
    adjacent polygons. By isolating and processing only the shared boundary segments
    (lines) before reconstructing the polygons, it prevents gaps, overlaps, or 
    topological disconnections from forming between neighboring features.

    Parameters
    ----------
    polygons : GeoDataFrame of Polygon
        The input layer containing adjacent polygons to be generalized.
    algorithm : callable
        The simplification or smoothing function to apply to the boundaries.
    *args : tuple
        Additional positional arguments forwarded to the `algorithm`.
    **kwargs : dict
        Additional keyword arguments forwarded to the `algorithm`.

    Returns
    -------
    GeoDataFrame of Polygon
    """
    polygons_simplified = polygons.copy()
    polygons_to_search = []

    # 1. EXTRACTION : Flatten MultiPolygons into simple Polygons for TopoMap
    for index, feature in polygons.iterrows():
        geom = feature['geometry']
        if geom.geom_type == 'MultiPolygon':
            polygons_to_search.extend(list(geom.geoms))
        else:
            polygons_to_search.append(geom)

    # First, create the TopoMap
    tm = TopoMap(snap_tolerance=1e-9, build_infinite_face=True)
    tm.build_from_polygons(polygons_to_search)

    # 2. GENERALIZATION : Apply the custom algorithm to each shared arc
    for aid, arc in tm.arcs.items():
        geom = arc.geom
        # Forwarding the geometry along with all extra arguments
        simplified = algorithm(geom, *args, **kwargs)
        tm.arcs[aid].geom = simplified

    # Update the TopoMap with the generalized arcs
    tm.clear_faces()
    # Rebuild the faces from scratch as some arcs may have changed or been removed
    tm.build_faces() 

    # Dictionary to collect all resulting faces belonging to each original row index
    matched_geometries = {idx: [] for idx in polygons_simplified.index}

    # 3. MATCHING : Match the new faces back to the original polygons
    for face in tm.faces:
        poly = tm.faces[face].polygon
        if poly is None:
            continue
            
        intersection = polygons_simplified.sindex.query(poly)
        dist_min = 1.5
        match = None
        
        for idx in intersection:
            candidate_geom = polygons_simplified.iloc[idx].geometry
            
            if candidate_geom.geom_type == 'MultiPolygon':
                for part in candidate_geom.geoms:
                    dist = surfacic_distance(poly, part)
                    if dist < dist_min:
                        dist_min = dist
                        match = idx
            else:
                dist = surfacic_distance(poly, candidate_geom)
                if dist < dist_min:
                    dist_min = dist
                    match = idx
        
        if match is not None:
            matched_geometries[match].append(poly)
    
    # 4. RECONSTRUCTION : Rebuild Polygons or MultiPolygons for the GeoDataFrame
    for idx, geoms in matched_geometries.items():
        if not geoms:
            continue 
        elif len(geoms) == 1:
            polygons_simplified.at[idx, 'geometry'] = geoms[0]
        else:
            polygons_simplified.at[idx, 'geometry'] = MultiPolygon(geoms)
    
    return polygons_simplified