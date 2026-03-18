import numpy as np
from shapely.geometry import LineString

def simplify_whirlpool(line, threshold):
    """
    Epsilon-circle based line simplification.

    This algorithm proposed by Dougenik and Chrisman :footcite:p:`dougenik:1979` performs a
    line simplification that removes spiky vertices while preserving the overall shape of the line. 
    It works by iterating through the vertices of the line and removing those that are within a specified
    distance (epsilon) from the last kept vertex.
    
    This method is particularly effective at simplifying lines with many small,
    sharp angles, such as rivers or coastlines, while maintaining the general form of the line.

    Parameters
    ----------
    line : LineString, MultiLineString
        The line to simplify.
    threshold : float
        The minimum epsilon-distance to consider a vertex to be removed.
        Higher values = fewer points kept (more aggressive simplification).

    Returns
    -------
    LineString, MultiLineString

    See Also
    --------
    simplify_douglas_peucker :
        Distance-based line simplification.
    simplify_lang :
        Look-ahead distance-based line simplification.
    simplify_li_openshaw :
        Square grid-based line simplification.
    simplify_raposo :
        Hexagon-based line simplification.
    simplify_reumann_witkam :
        Directional distance-based line simplification.
    simplify_visvalingam_whyatt :
        Area-based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> simplify_whirlpool(line, threshold=2.0)
    <LINESTRING (0 0, 5 3)>
    """
    if line.is_empty:
        return line

    if line.geom_type not in ['LineString', 'MultiLineString', 'LinearRing']:
        raise ValueError(f'{line.geom_type} geometry type cannot be simplified.')

    if line.geom_type == 'MultiLineString':
        geoms = [simplify_whirlpool(geom, threshold) for geom in line.geoms]
        return MultiLineString(geoms)

    coords = list(line.coords)
    if len(coords) < 2:
        return line

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