import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points
from collections import deque
from typing import Dict, Set, Tuple, List

def network_propagation(
    gdf: gpd.GeoDataFrame,
    modified_index: int,
    new_geometry: LineString,
    tolerance: float = 1e-6,
    max_propagation_distance: float = None,
    damping_factor: float = 0.1,
    inplace: bool = False
) -> gpd.GeoDataFrame:
    
    if not inplace:
        gdf = gdf.copy()

    # 1. Calcul du déplacement initial aux nœuds de la ligne modifiée
    old_geom = gdf.loc[modified_index, 'geometry']
    
    def get_vec(p1, p2):
        return (p2.x - p1.x, p2.y - p1.y)

    disp_start = get_vec(Point(old_geom.coords[0]), Point(new_geometry.coords[0]))
    disp_end = get_vec(Point(old_geom.coords[-1]), Point(new_geometry.coords[-1]))

    # Mise à jour de la ligne centrale
    gdf.at[modified_index, 'geometry'] = new_geometry

    # 2. Index de connectivité et identification des vraies bordures
    conn_map = {}
    for idx, row in gdf.iterrows():
        coords = row.geometry.coords
        for pos, pt_coords in [(0, coords[0]), (-1, coords[-1])]:
            pt_key = tuple(np.round(pt_coords, int(-np.log10(tolerance))))
            if pt_key not in conn_map:
                conn_map[pt_key] = []
            conn_map[pt_key].append((idx, pos == 0))

    # Points terminaux réels du réseau (degré 1)
    true_boundaries = {pt for pt, lines in conn_map.items() if len(lines) == 1}

    # 3. BFS avec gestion des "fausses bordures"
    queue = deque()
    visited = {modified_index}

    # Initialisation avec les voisins directs
    for pos, vec in [(0, disp_start), (-1, disp_end)]:
        if np.linalg.norm(vec) < tolerance: continue
        
        pt_key = tuple(np.round(old_geom.coords[pos], int(-np.log10(tolerance))))
        for neighbor_idx, is_start in conn_map.get(pt_key, []):
            if neighbor_idx not in visited:
                # On commence avec une distance de 0.0 pour les voisins directs
                queue.append((neighbor_idx, vec, is_start, 0.0))

    while queue:
        curr_idx, vec_anchor, is_start, dist = queue.popleft()
        if curr_idx in visited: continue
        
        line = gdf.loc[curr_idx, 'geometry']
        L = line.length
        # La distance cumulée est calculée à l'extrémité opposée de la ligne
        dist_at_end = dist + L
        
        # --- LOGIQUE DES BORDURES ---
        opp_pos = -1 if is_start else 0
        opp_pt_key = tuple(np.round(line.coords[opp_pos], int(-np.log10(tolerance))))
        
        # Une "fausse bordure" est activée si on dépasse la distance max
        is_fake_boundary = max_propagation_distance is not None and dist_at_end > max_propagation_distance
        is_true_boundary = opp_pt_key in true_boundaries
        
        if is_fake_boundary or is_true_boundary:
            # On force l'extrémité à rester immobile (amortissement total)
            target_damping = 0.0
            should_continue = False 
        else:
            # Propagation normale avec amortissement standard
            target_damping = damping_factor
            should_continue = True

        # Application de la transformation
        coords = np.array(line.coords)
        factors = np.linspace(1.0, target_damping, len(coords)) if is_start else np.linspace(target_damping, 1.0, len(coords))
        
        new_coords = coords.copy()
        new_coords[:, 0] += vec_anchor[0] * factors
        new_coords[:, 1] += vec_anchor[1] * factors
        
        gdf.at[curr_idx, 'geometry'] = LineString(new_coords)
        visited.add(curr_idx)

        # On ne propage aux voisins que si on n'a pas atteint une bordure (vraie ou fausse)
        if should_continue:
            vec_at_end = (vec_anchor[0] * target_damping, vec_anchor[1] * target_damping)
            for next_idx, next_is_start in conn_map.get(opp_pt_key, []):
                if next_idx not in visited:
                    queue.append((next_idx, vec_at_end, next_is_start, dist_at_end))

    return gdf

def network_propagation_batch(
    gdf: gpd.GeoDataFrame,
    modifications: List[Tuple[int, LineString]],
    tolerance: float = 1e-6,
    max_propagation_distance: float = None,
    damping_factor: float = 0.1,
    inplace: bool = False
) -> gpd.GeoDataFrame:
    """
    Applique une liste de modifications de manière séquentielle.
    
    Chaque modification est propagée dans le réseau résultant de la modification précédente.
    """
    result_gdf = gdf if inplace else gdf.copy()
    
    for i, (line_idx, new_geom) in enumerate(modifications):
        # On applique la propagation sur l'état actuel du GeoDataFrame
        result_gdf = network_propagation(
            result_gdf,
            modified_index=line_idx,
            new_geometry=new_geom,
            tolerance=tolerance,
            max_propagation_distance=max_propagation_distance,
            damping_factor=damping_factor,
            inplace=True # On travaille sur la copie en cours pour la séquence
        )
        
    return result_gdf