import numpy as np
from shapely.geometry import LineString, MultiLineString

def simplify_lang(line, tolerance: float, look_ahead: int=5) -> LineString:
    """
    Look-ahead distance-based line simplification.

    This algorithm, proposed by Lang :footcite:p:`lang:1969`, performs a 
    simplification by defining a search region of a fixed number of vertices 
    (look-ahead). It serves as a middle ground between local sequential filters 
    and global algorithms like Douglas-Peucker.

    The principle of the algorithm is to create a segment between the current 
    vertex and a vertex further down the line. 
    The perpendicular distances from all intermediate vertices to this segment 
    are calculated. If any distance exceeds the tolerance, the search region is 
    shrunk by moving the end vertex one step closer to the start, 
    and the process repeats until all intermediate points fall within the 
    tolerance. Once a valid segment is found, all intermediate points are 
    removed, and the process restarts from the end of that segment.

    Parameters
    ----------
    line : LineString, MultiLineString
        The line to simplify.
    tolerance : float
        The maximum allowed perpendicular distance between the original 
        vertices and the simplified segment.
    look_ahead : int, optional
        The maximum number of vertices to consider in a single search window. 
        Higher values allow for more aggressive simplification but increase 
        computational cost.

    Returns
    -------
    LineString, MultiLineString

    See Also
    --------
    simplify_douglas_peucker :
        Distance-based line simplification.
    simplify_li_openshaw :
        Square grid-based line simplification.
    simplify_raposo :
        Hexagon-based line simplification.
    simplify_reumann_witkam :
        Directional distance-based line simplification.
    simplify_visvalingam_whyatt :
        Area-based line simplification.
    simplify_whirlpool :
        Epsilon-circle based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 0.1), (2, 0.2), (3, 5), (4, 0)])
    >>> simplify_lang(line, tolerance=0.5, look_ahead=4)
    <LINESTRING (0 0, 2 0.2, 3 5, 4 0)>
    """
    if line.is_empty:
        return line
    
    # Gestion des MultiLineString
    if isinstance(line, MultiLineString):
        return MultiLineString([simplify_lang(l, tolerance, look_ahead) for l in line.geoms])

    coords = np.array(line.coords)
    n = len(coords)
    if n <= 2:
        return line

    mask = np.ones(n, dtype=bool)
    i = 0
    
    while i < n - 1:
        # Déterminer la fin de la fenêtre de recherche
        j = min(i + look_ahead, n - 1)
        
        # Réduire la fenêtre jusqu'à ce que tous les points intermédiaires soient valides
        while j > i + 1:
            p_start = coords[i]
            p_end = coords[j]
            
            # Points intermédiaires à tester
            p_mid = coords[i + 1 : j]
            
            # Calcul de distance point-segment (vecteur p_start -> p_end)
            v = p_end - p_start
            v_mag_sq = np.dot(v, v)
            
            if v_mag_sq == 0:
                # Segment nul, la distance est simplement la distance euclidienne
                dists = np.linalg.norm(p_mid - p_start, axis=1)
            else:
                # Distance perpendiculaire à la droite (p_start, p_end)
                # Formule : |(x2-x1)(y1-y3) - (x1-x3)(y2-y1)| / dist(p1,p2)
                areas = np.abs(v[0] * (p_start[1] - p_mid[:, 1]) - 
                               (p_start[0] - p_mid[:, 0]) * v[1])
                dists = areas / np.sqrt(v_mag_sq)
            
            if np.all(dists <= tolerance):
                # Tous les points intermédiaires sont dans la tolérance
                mask[i + 1 : j] = False
                i = j # On saute à la fin de cette fenêtre validée
                break
            else:
                # Au moins un point dépasse, on réduit la fenêtre par la fin
                j -= 1
        else:
            # Si j arrive à i + 1, on ne peut rien simplifier ici
            i += 1

    return LineString(coords[mask])