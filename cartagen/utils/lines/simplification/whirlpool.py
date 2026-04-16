import numpy as np
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

def simplify_whirlpool(geometry, threshold):
    """
    Simplify a line or polygon using an epsilon-circle based selection.

    This algorithm proposed by Dougenik and Chrisman :footcite:p:`dougenik:1979` performs a
    line simplification that removes spiky vertices while preserving the overall shape of the line. 
    It works by iterating through the vertices of the line and removing those that are within a specified
    distance (epsilon) from the last kept vertex.
    
    This method is particularly effective at simplifying lines with many small,
    sharp angles, such as rivers or coastlines, while maintaining the general form of the line.

    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to simplify.
    threshold : float
        The minimum epsilon-distance to consider a vertex to be removed.
        Higher values = fewer points kept (more aggressive simplification).

    Returns
    -------
    LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

    See Also
    --------
    simplify_angular :
        Simplify a line or polygon by removing vertexes with small angles.
    simplify_douglas_peucker :
        Simplify a line or polygon using a distance-based selection.
    simplify_lang :
        Simplify a line or polygon using a look-ahead distance-based selection.
    simplify_li_openshaw :
        Simplify a line or a polygon using a regular grid.
    simplify_raposo :
        Simplify a line or a polygon using an hexagonal tessellation.
    simplify_reumann_witkam :
        Simplify a line or polygon using a directional distance-based selection.
    simplify_topographic :
        Simplify a line or polygon and mimic hand-made cartographic generalization.
    simplify_visvalingam_whyatt :
        Simplify a line or polygon using an area-based selection.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> simplify_whirlpool(line, threshold=2.0)
    <LINESTRING (0 0, 5 3)>
    """
    # --- 1. Recursive handling for Multi-geometries ---
    if geometry.geom_type == 'MultiLineString':
        geoms = [simplify_whirlpool(geometry, threshold) for g in geometry.geoms]
        return MultiLineString(geoms)

    if geometry.geom_type == 'MultiPolygon':
        geoms = [simplify_whirlpool(geometry, threshold) for g in geometry.geoms]
        return MultiPolygon(geoms)

    # --- 2. Handling Polygons ---
    if geometry.geom_type == 'Polygon':
        # Simplify the exterior ring
        simplified_exterior = simplify_whirlpool(geometry.exterior, threshold)
        
        # VALIDITY CHECK: A LinearRing must have at least 4 coordinates
        if len(simplified_exterior.coords) < 4:
            # Option A: Return an empty geometry (it will be filtered out in the main loop)
            return Polygon() 
            # Option B: return geom (if you want to keep the original instead of deleting it)
        
        # Ensure it's a LinearRing (closed)
        exterior_ring = LinearRing(simplified_exterior.coords)

        # Simplify interior rings (holes)
        simplified_interiors = []
        for interior in geometry.interiors:
            s_int = simplify_whirlpool(interior, threshold)
            # Only keep holes that are still large enough to be rings
            if len(s_int.coords) >= 4:
                simplified_interiors.append(LinearRing(s_int.coords))
        
        poly = Polygon(exterior_ring, simplified_interiors)
        
        # Final topological repair
        if not poly.is_valid:
            poly = poly.buffer(0)
        return poly

    # --- 3. Core Linear Logic (for LineString, LinearRing, etc.) ---
    if geometry.geom_type not in ['LineString', 'LinearRing']:
        raise ValueError(f'{geometry.geom_type} geometry type cannot be simplified.')

    coords = list(geometry.coords)
    if len(coords) < 2:
        return geometry

    simplified_coords = [coords[0]]
    last_kept_pt = np.array(coords[0])

    for i in range(1, len(coords)):
        current_pt = np.array(coords[i])
        # Calcul de la distance euclidienne entre le dernier point gardé et le point actuel
        dist = np.linalg.norm(last_kept_pt - current_pt)
        
        # Selon WHIRLPOOL, on ne garde le point que s'il est à une distance >= d
        if dist >= threshold:
            simplified_coords.append(tuple(current_pt))
            last_kept_pt = current_pt
            
    # On s'assure que le dernier point de la ligne originale est conservé pour la forme
    if tuple(coords[-1]) not in simplified_coords:
        simplified_coords.append(coords[-1])

    return LineString(simplified_coords)