import math
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, LineString, Point
from shapely.ops import unary_union

def simplify_building_morphological(building, tolerance, max_iteration=100):
    """
    Simplify buildings using morphological closing and structural generalization.
    
    This algorithm performs a global structural simplification. For MultiPolygon, 
    it uses a morphological closing (buffer expansion and contraction using mitre joins) 
    to merge nearby structures and clear micro-gaps. For individual rings (exterior and interior), 
    it iteratively detects and processes small segments below the threshold, relocating vertices 
    to force parallel extensions or orthogonal alignments.

    This algorithm is a Python translation of the Lithuanian openmap association vector-map repository
    `simplification algorithm <https://github.com/openmaplt/vector-map/blob/master/db/func/stc_simplify_building_line.sql>`_.

    Parameters
    ----------
    building : Polygon or MultiPolygon
        The shapely geometry (building) to be simplified.
    tolerance : float
        The maximum distance tolerance for geometric simplification. Also used as the 
        buffer distance for morphological closing and as a square root threshold (tolerance^2) 
        for filtering out small interior rings.
    max_iteration : int, optional
        The maximum number of iterative vertex adjustments allowed during the line 
        simplification phase. The default value is set to 100.
    
    Returns
    -------
    Polygon or MultiPolygon
        The simplified building geometry of the same type as the input.

    See Also
    --------
    simplify_building_ruas :
        Simplify buildings by removing edges using Ruas algorithm.

    Examples
    --------
    >>> from shapely.geometry import Polygon
    >>> building = Polygon([(0, 0), (0, 10), (2, 10), (2, 9), (10, 9), (10, 0), (0, 0)])
    >>> simplify_building_morphological(building, 2.0)
    <POLYGON ((0 0, 0 9, 10 9, 10 0, 0 0))>
    """
    tp = building.geom_type

    if tp in ('LineString', 'LinearRing'):
        return simplify_boundary_morphological(building, tolerance, max_iteration=max_iteration)

    elif tp == 'Polygon':
        ext_line = LineString(building.exterior.coords)
        n_line = simplify_boundary_morphological(ext_line, tolerance, max_iteration=max_iteration)
        
        if len(n_line.coords) < 4:
            return Polygon()  # Le polygone s'est effondré
        
        simplified_interiors = []
        for interior in building.interiors:
            hole_poly = Polygon(interior)
            if hole_poly.area > (tolerance ** 2):
                e_line = simplify_boundary_morphological(LineString(interior.coords), tolerance, max_iteration=max_iteration)
                if len(e_line.coords) >= 4:
                    e_poly = Polygon(e_line)
                    if e_poly.area > (tolerance ** 2):
                        simplified_interiors.append(e_line)
                        
        return Polygon(n_line, simplified_interiors)

    elif tp == 'MultiPolygon':
        if len(building.geoms) == 1:
            simplified = simplify_building_morphological(building.geoms[0], tolerance, max_iteration)
            return MultiPolygon([simplified]) if simplified.geom_type == 'Polygon' else simplified
        else:
            buffered = building.buffer(tolerance, join_style='mitre')
            buffered_back = buffered.buffer(-tolerance, join_style='mitre')
            n = buffered_back.union(building)
            
            if n.geom_type == 'Polygon':
                n = MultiPolygon([n])
                
            if len(n.geoms) == len(building.geoms):
                n = building
                
            g = []
            for geom in n.geoms:
                if geom.area > max_iteration:
                    simplified_part = simplify_building_morphological(geom, tolerance, max_iteration)
                    if not simplified_part.is_empty:
                        g.append(simplified_part)
                        
            if not g:
                return MultiPolygon()
                
            union_g = unary_union(g)
            if union_g.geom_type == 'Polygon':
                return MultiPolygon([union_g])
            return union_g
    else:
        return building

def simplify_boundary_morphological(line_geom, tolerance, max_iteration=100, debug=False):
    """
    Version Python de simplify_boundary_morphological.
    Généralise une ligne ou un anneau de bâtiment en préservant ses angles.
    """
    is_closed = line_geom.coords[0] == line_geom.coords[-1]
    current_coords = clean_coords(list(line_geom.coords), is_closed)
    excluded_edges = []

    for _ in range(max_iteration):
        num_pts = len(current_coords)
        if num_pts <= 5:  # Un anneau fermé a besoin d'au moins 4-5 points
            break

        smallest = 1000000.0
        sn = -1

        # Recherche du plus petit segment non exclu inférieur à la tolérance `t`
        for i in range(2, num_pts + 1):
            p_start = get_pt_sql(current_coords, i - 1, is_closed)
            p_end = get_pt_sql(current_coords, i, is_closed)
            edge_line = LineString([p_start, p_end])
            length = edge_line.length

            if length < smallest and length <= tolerance:
                is_ex = any(edge_line.distance(ex) < 0.1 for ex in excluded_edges)
                if not is_ex:
                    smallest = length
                    sn = i

        if sn == -1:
            break  # Plus aucun segment à simplifier

        backup_coords = list(current_coords)

        # Récupération des sommets voisins (indexation calquée sur le SQL)
        p_sn_minus_3 = get_pt_sql(current_coords, sn - 3, is_closed)
        p_sn_minus_2 = get_pt_sql(current_coords, sn - 2, is_closed)
        p_sn_minus_1 = get_pt_sql(current_coords, sn - 1, is_closed)
        p_sn         = get_pt_sql(current_coords, sn, is_closed)
        p_sn_plus_1  = get_pt_sql(current_coords, sn + 1, is_closed)
        p_sn_plus_2  = get_pt_sql(current_coords, sn + 2, is_closed)

        aze = get_azimuth(p_sn_minus_1, p_sn)
        azp = get_azimuth(p_sn, p_sn_plus_1)
        lp  = LineString([p_sn, p_sn_plus_1]).length

        azm = get_azimuth(p_sn_minus_2, p_sn_minus_1)
        lm  = LineString([p_sn_minus_2, p_sn_minus_1]).length

        # Calcul de l'angle de déviation pour intrusion/extrusion
        ac = azp - aze
        if ac > 180: ac -= 360
        elif ac < -180: ac += 360
        intrusion = ac < 0

        # Calcul du changement d'angle global entre les voisins
        ac_neighbors = azp - azm
        if ac_neighbors > 180: ac_neighbors -= 360
        elif ac_neighbors < -180: ac_neighbors += 360

        modified = False
        np_pt = None

        # -----------------------------------------------------------------
        # CAS 1 : BRANCHE CHANGE ~0 (Lignes presque parallèles)
        # -----------------------------------------------------------------
        if -40 <= ac_neighbors <= 40:
            if lm > lp:
                l1 = LineString([p_sn_minus_1, project_pt(p_sn_minus_1, 1000, azm)])
                ew = LineString([p_sn_plus_1, p_sn_plus_2])
                inter = l1.intersection(ew)
                if inter.is_empty or inter.geom_type != 'Point':
                    az_next = get_azimuth(p_sn_plus_1, p_sn_plus_2)
                    l2 = LineString([p_sn_plus_1, project_pt(p_sn_plus_1, -1000, az_next)])
                    inter = l1.intersection(l2)

                if not inter.is_empty and inter.geom_type == 'Point':
                    np_pt = (inter.x, inter.y)
                    set_point_sql(current_coords, sn - 1, np_pt, is_closed)
                    set_point_sql(current_coords, sn, np_pt, is_closed)
                    modified = True
            else:
                l1 = LineString([p_sn, project_pt(p_sn, -1000, azp)])
                ew = LineString([p_sn_minus_2, p_sn_minus_3])
                inter = l1.intersection(ew)
                if inter.is_empty or inter.geom_type != 'Point':
                    az_prev = get_azimuth(p_sn_minus_2, p_sn_minus_3)
                    l2 = LineString([p_sn_minus_2, project_pt(p_sn_minus_2, -1000, az_prev)])
                    inter = l1.intersection(l2)

                if not inter.is_empty and inter.geom_type == 'Point':
                    np_pt = (inter.x, inter.y)
                    set_point_sql(current_coords, sn - 1, np_pt, is_closed)
                    set_point_sql(current_coords, sn - 2, np_pt, is_closed)
                    set_point_sql(current_coords, sn - 3, np_pt, is_closed)
                    modified = True

        # -----------------------------------------------------------------
        # CAS 2 : BRANCHE CHANGE 90 (Angles droits - Le plus fréquent sur le bâti)
        # -----------------------------------------------------------------
        elif (50 <= ac_neighbors <= 110) or (-110 <= ac_neighbors <= -50):
            l1 = LineString([p_sn_minus_2, project_pt(p_sn_minus_1, 100, azm)])
            l2 = LineString([p_sn_plus_1, project_pt(p_sn, -1000, azp)])
            inter = l1.intersection(l2)
            if not inter.is_empty and inter.geom_type == 'Point':
                np_pt = (inter.x, inter.y)
                set_point_sql(current_coords, sn - 1, np_pt, is_closed)
                set_point_sql(current_coords, sn - 2, np_pt, is_closed)
                modified = True

        # Note: Les cas complexes de changement à 180° (Int/Ext) du SQL original
        # n'étant pas finalisés ("TODO A/B/C/D"), ils tombent dans la section non-gérée.

        # -----------------------------------------------------------------
        # VALIDATION DE LA TOPOLOGIE ET NETTOYAGE
        # -----------------------------------------------------------------
        if not modified:
            edge_geom = LineString([p_sn_minus_1, p_sn])
            excluded_edges.append(edge_geom)
        else:
            current_coords = clean_coords(current_coords, is_closed)
            
            # Sécurité anti-effondrement : On s'assure que la géométrie reste valide
            valid_topology = True
            if is_closed:
                if len(current_coords) < 4:
                    valid_topology = False
                else:
                    poly_test = Polygon(current_coords)
                    if not poly_test.is_valid or poly_test.area < (tolerance ** 2):
                        valid_topology = False
            else:
                if len(current_coords) < 2:
                    valid_topology = False

            if not valid_topology:
                # Rollback si la topologie est cassée
                current_coords = backup_coords
                edge_geom = LineString([p_sn_minus_1, p_sn])
                excluded_edges.append(edge_geom)

    return LineString(current_coords)

def get_azimuth(p1, p2):
    """Équivalent de ST_Azimuth (0 = Nord, 90 = Est, en degrés)"""
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    az = math.degrees(math.atan2(dx, dy))
    return az % 360

def project_pt(p, dist, az_deg):
    """Équivalent de ST_Project en coordonnées cartésiennes plane"""
    az_rad = math.radians(az_deg)
    dx = dist * math.sin(az_rad)
    dy = dist * math.cos(az_rad)
    return (p[0] + dx, p[1] + dy)

def clean_coords(coords, is_closed):
    """
    Équivalent combiné de stc_simplify_turbo, stc_simplify_angle et stc_remove_spike.
    Nettoie les doublons, les pointes (allers-retours) et les points colinéaires.
    """
    if len(coords) < 3:
        return coords

    # Retrait des doublons consécutifs
    unique_pts = []
    for pt in coords:
        if not unique_pts or pt != unique_pts[-1]:
            unique_pts.append(pt)
            
    if is_closed and len(unique_pts) > 1 and unique_pts[0] == unique_pts[-1]:
        unique_pts.pop()

    changed = True
    while changed:
        changed = False
        n = len(unique_pts)
        if n < 3:
            break
        filtered = []
        for i in range(n):
            p_prev = unique_pts[i - 1]
            p_curr = unique_pts[i]
            p_next = unique_pts[(i + 1) % n]

            # Supprime les "spikes" (allers-retours sur le même point)
            if p_prev == p_next:
                changed = True
                continue

            # Supprime les points colinéaires (angle de ~0 ou ~180 degrés)
            az1 = get_azimuth(p_prev, p_curr)
            az2 = get_azimuth(p_curr, p_next)
            diff = (az2 - az1) % 360
            if diff > 180:
                diff = 360 - diff

            if diff < 1.0 or diff > 179.0:
                changed = True
                continue

            filtered.append(p_curr)
        unique_pts = filtered

    if is_closed and unique_pts:
        unique_pts.append(unique_pts[0])
    return unique_pts

def get_pt_sql(coords_list, sql_idx, is_closed):
    """Permet de requêter la liste Python en utilisant la logique d'index 1-basé de PostGIS"""
    n = len(coords_list)
    idx_0b = sql_idx - 1  # Conversion 1-basé vers 0-basé
    if is_closed:
        u_len = n - 1  # Exclure le point de fermeture pour le modulo
        return coords_list[idx_0b % u_len]
    else:
        return coords_list[max(0, min(idx_0b, n - 1))]

def set_point_sql(coords_list, sql_idx, new_pt, is_closed):
    """Modifie un point dans la liste en utilisant la logique d'index 0-basé de ST_SetPoint"""
    n = len(coords_list)
    if is_closed:
        u_len = n - 1
        actual_idx = sql_idx % u_len
        coords_list[actual_idx] = new_pt
        if actual_idx == 0:
            coords_list[u_len] = new_pt  # Mettre à jour la fermeture
    else:
        if 0 <= sql_idx < n:
            coords_list[sql_idx] = new_pt