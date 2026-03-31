import numpy as np
from shapely.geometry import LineString, MultiLineString

def simplify_reumann_witkam(line: LineString, tolerance: float) -> LineString:
    """
    Directional distance-based line simplification.

    This algorithm, proposed by Reumann and Witkam :footcite:p:`reumann:1974`, performs a 
    sequential line simplification by using a "corridor" or "tube" defined by the 
    direction of the first segment. Unlike the Douglas-Peucker algorithm, which 
    considers the line globally, Reumann-Witkam is a local, streaming-friendly 
    filter that processes vertices in order.

    The principle of the algorithm is to define a search pipe using the first two 
    points of a segment. For all subsequent points, the perpendicular distance to 
    the infinite line passing through this initial segment is calculated. 
    As long as the points stay within the tolerance distance, they are marked for 
    deletion. When a point falls outside the pipe, the current point becomes the 
    new starting vertex, and a new pipe direction is established.

    The algorithm is particularly efficient for reducing the density of points in 
    datasets where the direction of the line is relatively constant, making it 
    ideal for real-time thinning of trajectory data or GPS traces.

    Parameters
    ----------
    line : LineString, MultiLineString
        The line to simplify.
    tolerance : float
        The width (radius) of the search corridor. Points within this distance 
        from the segment's trajectory are removed.
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
    simplify_visvalingam_whyatt :
        Area-based line simplification.
    simplify_whirlpool :
        Epsilon-circle based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (10, 0.1), (20, -0.1), (25, 10)])
    >>> simplify_reumann_witkam(line, tolerance=0.5)
    <LINESTRING (0 0, 20 -0.1, 25 10)>
    """
    # Gestion des MultiLineString
    if isinstance(line, MultiLineString):
        return MultiLineString([simplify_reumann_witkam(l, tolerance) for l in line.geoms])

    # Extraction des coordonnées en une seule fois (très rapide)
    coords = np.array(line.coords)
    n = len(coords)
    
    if n <= 2:
        return line

    # Initialisation du masque (True = on garde le point)
    mask = np.ones(n, dtype=bool)
    
    # On travaille sur des sous-vues pour éviter les copies de mémoire
    # i est l'indice du début du segment de référence
    i = 0
    while i < n - 2:
        p1 = coords[i]
        p2 = coords[i+1]
        
        # Vecteur directeur p1 -> p2
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        mag_sq = dx*dx + dy*dy
        
        # Sécurité : si p1 et p2 sont identiques, on passe au suivant
        if mag_sq == 0:
            i += 1
            continue
            
        inv_mag = 1.0 / np.sqrt(mag_sq)
        
        # On teste les points suivants
        j = i + 2
        while j < n:
            p3 = coords[j]
            # Distance point-droite (p1, p2) : |(x2-x1)(y1-y3) - (x1-x3)(y2-y1)| / mag
            dist = abs(dx * (p1[1] - p3[1]) - (p1[0] - p3[0]) * dy) * inv_mag
            
            if dist <= tolerance:
                mask[j-1] = False  # On élimine le point intermédiaire
                j += 1
            else:
                # Le point sort du "tube", on recommence un segment à partir d'ici
                i = j - 1
                break
        else:
            # Fin de la ligne atteinte
            break

    # Création de la nouvelle LineString à partir des points filtrés
    return LineString(coords[mask])