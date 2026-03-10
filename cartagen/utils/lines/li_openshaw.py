import shapely
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString

from cartagen.utils.partitioning import partition_grid

def li_openshaw(line, cell_size):
    """
    Regular grid-based line simplification.

    This algorithm proposed by Li & Openshaw :footcite:p:`li:1993` simplifies lines based on a
    regular square grid. It first divide the line vertexes into groups partionned by a regular
    grid, then each group of vertexes is replaced by their centroid.

    Parameters
    ----------
    line : LineString, MultiLineString
        The line to simplify.
    cell_size : float
        The size of the regular grid used to divide the line.

    Returns
    -------
    LineString, MultiLineString

    See Also
    --------
    douglas_peucker :
        Distance-based line simplification.
    visvalingam_whyatt :
        Area-based line simplification.
    raposo :
        Hexagon-based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.li_openshaw(line, 1)
    <LINESTRING (0 0, 0.5 0.5, 2 0, 5 3)>
    """
    if line.geom_type not in ['LineString', 'MultiLineString']:
        raise ValueError(f'{line.geom_type} geometry type cannot be simplified.')
    
    if line.geom_type == 'MultiLineString':
        geoms = [li_openshaw(geom, cell_size) for geom in line.geoms]
        return MultiLineString(geoms)
    
    vertexes = [Point(x) for x in line.coords]
    gdf = gpd.GeoDataFrame(geometry=vertexes)
    groups, squares = partition_grid(gdf, cell_size)
    
    # Créer un mapping vertex_index -> group_index pour O(1) lookup
    vertex_to_group = {}
    for gi, group in enumerate(groups):
        for vertex_idx in group:
            vertex_to_group[vertex_idx] = gi
    
    simplified = []
    processed_groups = set()
    
    for i in range(len(vertexes)):
        group_idx = vertex_to_group.get(i)
        if group_idx is not None and group_idx not in processed_groups:
            group = groups[group_idx]
            centroid = shapely.centroid(MultiPoint([vertexes[v] for v in group]))
            simplified.append(centroid)
            processed_groups.add(group_idx)
    
    # Assurer la préservation des extrémités
    if not simplified[0].equals(vertexes[0]):
        simplified.insert(0, vertexes[0])
    if not simplified[-1].equals(vertexes[-1]):
        simplified.append(vertexes[-1])
    
    return LineString(simplified)