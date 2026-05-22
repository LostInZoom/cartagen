import math
import numpy as np
from shapely.geometry import Polygon

def square_polygon_orientation(polygon):
    """
    Squares a building by forcing it into an orthogonal coordinate system.
    
    This algorithm described by Hakim and Tsai :footcite:p:`hakim:2026` is a heuristic approach
    that leverages the Manhattan World assumption. It identifies the 
    primary structural axis of the building (its longest wall segment) and performs a 
    global rotation to align the footprint with the X/Y grid. Facades are then classified 
    and merged into strictly horizontal or vertical lines. Perfect 90-degree corners are 
    reconstructed via direct line intersections before an inverse rotation is applied to 
    restore the building to its original geospatial location.

    Parameters
    ----------
    polygon : Polygon
        The polygon to square.

    Returns
    -------
    Polygon

    See Also
    --------
    square_polygon_ls :
        Squares a polygon using the method of least squares.
    square_polygon_naive:
        Squares a polygon according to the statistical weighted orientation (SWO).

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 1), (1.1, 1), (1, 0)])
    >>> square_polygon_orientation(polygon)
    <POLYGON ((0 1, 1.05 1, 1.05 0, 0 0, 0 1))>        
    """
    # Extraction des coordonnées uniques (on retire le point de fermeture dupliqué à la fin)
    coords = np.array(polygon.exterior.coords)
    if np.allclose(coords[0], coords[-1]):
        coords = coords[:-1]
        
    N = len(coords)
    if N < 3:
        return polygon
        
    # --- ÉTAPE 1 : Trouver le segment le plus long (Orientation dominante) ---
    max_len = -1
    dominant_angle = 0
    for i in range(N):
        p1 = coords[i]
        p2 = coords[(i + 1) % N]
        v = p2 - p1
        length = np.linalg.norm(v)
        if length > max_len:
            max_len = length
            dominant_angle = math.atan2(v[1], v[0]) # Angle en radians
            
    # --- ÉTAPE 2 : Rotation globale pour aligner le bâtiment sur l'axe X ---
    cos_a, sin_a = math.cos(-dominant_angle), math.sin(-dominant_angle)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    coords_rot = coords @ R.T
    
    # --- ÉTAPE 3 : Classer les segments en Horizontaux (H) ou Verticaux (V) ---
    segments = []
    for i in range(N):
        p1 = coords_rot[i]
        p2 = coords_rot[(i + 1) % N]
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        
        # Si le déplacement horizontal est majeur, c'est un mur horizontal
        if abs(dx) >= abs(dy):
            segments.append({'type': 'H', 'pts': [p1, p2]})
        else:
            segments.append({'type': 'V', 'pts': [p1, p2]})
            
    # --- ÉTAPE 4 : Fusionner les segments consécutifs de même direction ---
    # (Si deux murs horizontaux se suivent, ils ne forment qu'une seule grande façade)
    merged_segments = [segments[0]]
    for i in range(1, len(segments)):
        if segments[i]['type'] == merged_segments[-1]['type']:
            merged_segments[-1]['pts'].extend(segments[i]['pts'])
        else:
            merged_segments.append(segments[i])
            
    # Gestion du rebouclage cyclique (entre le premier et le dernier segment)
    if len(merged_segments) > 1 and merged_segments[0]['type'] == merged_segments[-1]['type']:
        merged_segments[0]['pts'].extend(merged_segments[-1]['pts'])
        merged_segments.pop(-1)
        
    # Un bâtiment orthogonal valide doit avoir au moins 4 façades (un rectangle)
    if len(merged_segments) < 4:
        return polygon 
        
    # --- ÉTAPE 5 : Calculer la position moyenne de chaque façade ---
    for seg in merged_segments:
        pts = np.array(seg['pts'])
        if seg['type'] == 'H':
            seg['val'] = np.mean(pts[:, 1]) # Hauteur Y moyenne de la façade
        else:
            seg['val'] = np.mean(pts[:, 0]) # Position X moyenne de la façade
            
    # --- ÉTAPE 6 : Reconstruction des coins par intersection alternée ---
    new_coords_rot = []
    M = len(merged_segments)
    for i in range(M):
        seg_curr = merged_segments[i]
        seg_next = merged_segments[(i + 1) % M]
        
        # L'intersection d'une ligne Y=constante et X=constante est immédiate
        if seg_curr['type'] == 'H' and seg_next['type'] == 'V':
            new_coords_rot.append([seg_next['val'], seg_curr['val']])
        elif seg_curr['type'] == 'V' and seg_next['type'] == 'H':
            new_coords_rot.append([seg_curr['val'], seg_next['val']])
            
    if not new_coords_rot:
        return polygon
        
    new_coords_rot = np.array(new_coords_rot)
    
    # --- ÉTAPE 7 : Rotation inverse pour remettre le bâtiment à sa place ---
    cos_inv, sin_inv = math.cos(dominant_angle), math.sin(dominant_angle)
    R_inv = np.array([[cos_inv, -sin_inv], [sin_inv, cos_inv]])
    new_coords = new_coords_rot @ R_inv.T
    
    # Fermeture propre du polygone Shapely
    new_coords = np.vstack([new_coords, new_coords[0]])
    return Polygon(new_coords)