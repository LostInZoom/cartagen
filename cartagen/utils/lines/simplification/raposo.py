from shapely.ops import transform, nearest_points
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

from cartagen.utils.partitioning.tessellation import HexagonalTessellation

def simplify_raposo(geometry, initial_scale, final_scale, centroid=True, tobler=False):
    """
    Simplify a line or a polygon using an hexagonal tessellation.
    
    This algorithm proposed by Raposo :footcite:p:`raposo:2013` simplifies lines based on a
    hexagonal tessellation. The algorithm also works for the simplification of the border of a polygon object.
    The idea of the algorithm is to put a hexagonal tessellation on top of the line to simplify,
    the size of the cells depending on the targeted granularity of the line.
    Similarly to the Li-Openshaw algorithm, only one vertex is kept inside each cell.
    This point can be the centroid of the removed vertices, or a projection on the initial line of this centroid.
    The shapes obtained with this algorithm are less sharp than the ones obtained with other algorithms such as Douglas-Peucker.

    The algorithm is dedicated to the smooth simplification of natural features such as rivers, forests, coastlines, lakes.

    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The line to simplify.
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
    simplify_reumann_witkam :
        Simplify a line or polygon using a directional distance-based selection.
    simplify_topographic :
        Simplify a line or polygon and mimic hand-made cartographic generalization.
    simplify_visvalingam_whyatt :
        Simplify a line or polygon using an area-based selection.
    simplify_whirlpool :
        Simplify a line or polygon using an epsilon-circle based selection.

    Notes
    -----
    The Tobler based formula to compute hexagonal cell size is :math:`cellsize=5·l·s`
    where :math:`l` is the width of the line in the map in meters
    (*e.g.* 0.0005 for 0.5 mm), and :math:`s` is the target scale denominator.

    Raposo’s formula to compute hexagonal cell size is :math:`cellsize={l/n}·{t/d}`
    where :math:`l` is the length of the line, :math:`n` the number of vertices of the line,
    :math:`t` the denominator of the target scale, and :math:`d` the denominator of the initial scale.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.simplify_raposo(line, 5000, 10000)
    <LINESTRING (0 0, 1 0.3333333333333333, 5 3)>
    """

    # --- 1. Recursive handling for Multi-geometries ---
    if geometry.geom_type == 'MultiLineString':
        geoms = [simplify_raposo(g, initial_scale, final_scale, centroid, tobler) for g in geometry.geoms]
        return MultiLineString(geoms)

    if geometry.geom_type == 'MultiPolygon':
        geoms = [simplify_raposo(g, initial_scale, final_scale, centroid, tobler) for g in geometry.geoms]
        return MultiPolygon(geoms)

    # --- 2. Handling Polygons ---
    if geometry.geom_type == 'Polygon':
        # Simplify the exterior ring
        simplified_exterior = simplify_raposo(geometry.exterior, initial_scale, final_scale, centroid, tobler)
        
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
            s_int = simplify_raposo(interior, initial_scale, final_scale, centroid, tobler)
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

    # Calculate hexagons width
    if tobler:
        width = final_scale * 5 / 4000 
    else:
        first_factor = geometry.length / (len(geometry.coords) - 1)
        second_factor = final_scale / initial_scale
        width = first_factor * second_factor

    # Create the hex tesselation
    tessellation = HexagonalTessellation(geometry.envelope, width)
    
    # Convert coords to 2d
    coords_2d = [transform(lambda *args: args[:2], Point(c)).coords[0] for c in geometry.coords]
    
    # Initialize with first point
    final_coords = [coords_2d[0]]
    
    i = 1
    previous_cell = None
    
    while i < len(coords_2d):
        # Find the cell containing the current point
        containing_cells = tessellation.get_containing_cells(coords_2d[i])
        
        # Ignore previous cell if present
        if previous_cell and previous_cell in containing_cells:
            containing_cells.remove(previous_cell)
        
        if not containing_cells:
            i += 1
            continue
            
        current_cell = containing_cells[0]
        
        # Collect points inside the cell
        point_cloud = []
        j = i
        
        while j < len(coords_2d):
            cells = tessellation.get_containing_cells(coords_2d[j])
            
            if current_cell in cells:
                point_cloud.append(coords_2d[j])
                j += 1
            else:
                break
        
        # Add the point representative of the cell
        if point_cloud:
            multipoint = MultiPoint(point_cloud)
            
            if centroid:
                final_coords.append(multipoint.centroid.coords[0])
            else:
                # Find the point closest to centroid
                nearest = nearest_points(multipoint, multipoint.centroid)
                final_coords.append(nearest[0].coords[0])
        
        previous_cell = current_cell
        i = j
    
    # Logic for closing the ring or ensuring the last point
    is_ring = geometry.geom_type == 'LinearRing' or (coords_2d[0] == coords_2d[-1])
    
    if is_ring:
        # If it was a ring, ensure the last point is exactly the same as the first
        if final_coords[-1] != final_coords[0]:
            final_coords.append(final_coords[0])
    elif coords_2d[-1] != final_coords[-1]:
        final_coords.append(coords_2d[-1])
        
    return LineString(final_coords)