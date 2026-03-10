import shapely
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString

from cartagen.utils.partitioning import partition_grid

def li_openshaw(line, cell_size, preserve_extremities=True):
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
    preserve_extremities : bool, optional
        Whether the algorithm should preserve the first and last vertex
        of the input line.

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
    
    coords = np.array(line.coords)
    n = len(coords)
    
    if n <= 2:
        return line
    
    # Calculate cell indexes
    cell_indices = (coords // cell_size).astype(int)
    
    # Create unique id for each cell
    cell_ids = cell_indices[:, 0] * 1000000 + cell_indices[:, 1]
    
    # Find unique cell while preserving order
    _, unique_idx = np.unique(cell_ids, return_index=True)
    unique_idx = np.sort(unique_idx)
    
    # Calculate centroid for each cell
    simplified = []
    for idx in unique_idx:
        mask = cell_ids == cell_ids[idx]
        centroid = coords[mask].mean(axis=0)
        simplified.append(tuple(centroid))
    
    if preserve_extremities:
        # Keep extremities
        if simplified[0] != tuple(coords[0]):
            simplified.insert(0, tuple(coords[0]))
        if simplified[-1] != tuple(coords[-1]):
            simplified.append(tuple(coords[-1]))
    
    return LineString(simplified)