from shapely.ops import transform, nearest_points
from shapely import Point, MultiPoint, LineString, MultiLineString

from cartagen.utils.partitioning.tessellation import HexagonalTessellation

def raposo(line, initial_scale, final_scale, centroid=True, tobler=False):
    """
    Hexagon-based line simplification.
    
    This algorithm proposed by Raposo :footcite:p:`raposo:2013` simplifies lines based on a
    hexagonal tessellation. The algorithm also works for the simplification of the border of a polygon object.
    The idea of the algorithm is to put a hexagonal tessallation on top of the line to simplify,
    the size of the cells depending on the targeted granularity of the line.
    Similarly to the Li-Openshaw algorithm, only one vertex is kept inside each cell.
    This point can be the centroid of the removed vertices, or a projection on the initial line of this centroid.
    The shapes obtained with this algorithm are less sharp than the ones obtained with other algorithms such as Douglas-Peucker.

    The algorithm is dedicated to the smooth simplification of natural features such as rivers, forests, coastlines, lakes.

    Parameters
    ----------
    line : LineString, MultiLineString
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
    LineString, MultiLineString

    See Also
    --------
    douglas_peucker :
        Distance-based line simplification.
    visvalingam_whyatt :
        Area-based line simplification.
    li_openshaw :
        Square grid-based line simplification.

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
    >>> c4.raposo(line, 5000, 10000)
    <LINESTRING (0 0, 1 0.3333333333333333, 5 3)>
    """

    if line.geom_type not in ['LineString', 'MultiLineString']:
        raise ValueError(f'{line.geom_type} geometry type cannot be simplified.')

    if line.geom_type == 'MultiLineString':
        geoms = [raposo(geom, initial_scale, final_scale, centroid, tobler) for geom in line.geoms]
        return MultiLineString(geoms)

    # Calculate hexagons width
    if tobler:
        width = final_scale * 5 / 4000 
    else:
        first_factor = line.length / (len(line.coords) - 1)
        second_factor = final_scale / initial_scale
        width = first_factor * second_factor

    # Create the hex tesselation
    tessellation = HexagonalTessellation(line.envelope, width)
    
    # Convert coords to 2d
    coords_2d = [transform(lambda *args: args[:2], Point(c)).coords[0] for c in line.coords]
    
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
    
    # Add the last point if not present
    if coords_2d[-1] != final_coords[-1]:
        final_coords.append(coords_2d[-1])
    
    return LineString(final_coords)