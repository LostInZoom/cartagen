import shapely
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

from cartagen.utils.partitioning import partition_grid

def simplify_li_openshaw(geometry, cell_size, preserve_extremities=True):
    """
    Simplify a line or a polygon using a regular grid.

    This algorithm proposed by Li & Openshaw :footcite:p:`li:1993` simplifies lines based on a
    regular square grid. It first divide the line vertexes into groups partionned by a regular
    grid, then each group of vertexes is replaced by their centroid.

    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to simplify.
    cell_size : float
        The size of the regular grid used to divide the line.
    preserve_extremities : bool, optional
        Whether the algorithm should preserve the first and last vertex
        of the input geometry.

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
    simplify_raposo :
        Simplify a line or a polygon using an hexagonal tessellation.
    simplify_reumann_witkam :
        Simplify a line or polygon using a directional distance-based selection.
    simplify_topographic :
        Simplify a line or polygon and mimic hand-made cartographic generalization.
    simplify_visvalingam_whyatt :
        Simplify a line or polygon using an area-based selection.
    simplify_whirlpool :
        Simplify a line or polygon using an epsilon-circle based selection.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.simplify_li_openshaw(line, 1)
    <LINESTRING (0 0, 0.5 0.5, 2 0, 5 3)>
    """
    # --- 1. Recursive handling for Multi-geometries ---
    if geometry.geom_type == 'MultiLineString':
        geoms = [simplify_li_openshaw(geometry, cell_size, preserve_extremities) for g in geometry.geoms]
        return MultiLineString(geoms)

    if geometry.geom_type == 'MultiPolygon':
        geoms = [simplify_li_openshaw(geometry, cell_size, preserve_extremities) for g in geometry.geoms]
        return MultiPolygon(geoms)

    # --- 2. Handling Polygons ---
    if geometry.geom_type == 'Polygon':
        # Simplify the exterior ring
        simplified_exterior = simplify_li_openshaw(geometry.exterior, cell_size, preserve_extremities=False)
        
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
            s_int = simplify_li_openshaw(interior, cell_size, preserve_extremities=False)
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
    
    coords = np.array(geometry.coords)
    n = len(coords)
    
    if n <= 2:
        return geometry
    
    # Calculate cell index for each vertex
    cell_indices = (coords // cell_size).astype(int)

    # Build a map of cell_id -> list of vertices, preserving traversal order.
    # Crucially, if the line re-enters a cell it already visited, that second
    # run is treated as a NEW group (matching the Java reference logic).
    # This prevents vertices from non-adjacent parts of the line being averaged
    # together, which is the root cause of self-intersections.
    simplified = []
    seen_cells = set()
    current_cell_id = None
    current_group = []

    for i, (ci, coord) in enumerate(zip(cell_indices, coords)):
        cell_id = (ci[0], ci[1])

        if cell_id == current_cell_id:
            # Still in the same cell: accumulate
            current_group.append(coord)
        else:
            # Flush the previous group
            if current_group:
                centroid = np.mean(current_group, axis=0)
                simplified.append(tuple(centroid))

            # Start a new group for this cell.
            # If we have already visited this cell earlier in the line
            # (re-entry), we still start a fresh group — do NOT merge with
            # the earlier visit, as that would connect distant parts of the
            # line and risk self-intersections.
            current_cell_id = cell_id
            current_group = [coord]

    # Flush the last group
    if current_group:
        centroid = np.mean(current_group, axis=0)
        simplified.append(tuple(centroid))

    if preserve_extremities:
        if simplified[0] != tuple(coords[0]):
            simplified.insert(0, tuple(coords[0]))
        if simplified[-1] != tuple(coords[-1]):
            simplified.append(tuple(coords[-1]))
    
    return LineString(simplified)