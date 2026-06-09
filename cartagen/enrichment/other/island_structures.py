import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
import pandas as pd
from shapely.geometry import Point
from scipy.spatial import KDTree
import networkx as nx

def detect_polygon_patches_structures(gdf, k=3.0, radius_method='area', min_line_length=3):
    """
    Implementation of the algorithm of Steiniger et al. :footcite:p:`steiniger:2006` for grouping 
    polygon patches (islands/lakes/forests) into meso-structures. The algorithm is based on the concept 
    of "perception horizon" and a proximity graph (Minimum Spanning Tree).
    
    Parameters:
    -----------
    gdf : geopandas.GeoDataFrame
        The GeoDataFrame containing the polygons (indexed from 0 to N-1).
    k : float
        The multiplier simulating the "calculation of the perception horizon" (default 3.0).
    radius_method : str
        'area' (radius of an equivalent circle) or 'mbr' (half of the longest side of the MBR).
    
    Returns
    -------
    geopandas.GeoDataFrame
        The GeoDataFrame with an additional column 'group_id' indicating the group to which each polygon belongs.

    References
    ----------
    .. footbibliography::
    """
    gdf = gdf.reset_index(drop=True).copy()
    n_polygons = len(gdf)
    
    # --- STEP 1 : Calculation of the perception horizons (d_max) ---
    r = np.zeros(n_polygons)
    for idx, poly in enumerate(gdf.geometry):
        if radius_method == 'area':
            r[idx] = np.sqrt(poly.area / np.pi)
        elif radius_method == 'mbr':
            mbr = poly.minimum_rotated_rectangle
            x, y = mbr.exterior.coords.xy
            edge_lengths = [np.hypot(x[i+1]-x[i], y[i+1]-y[i]) for i in range(3)]
            r[idx] = max(edge_lengths) / 2.0
            
    d_max = k * r

    # --- STEP 2 : Construction du Proximity Graph (Two-way edge) ---
    G_prox = nx.Graph()
    G_prox.add_nodes_from(range(n_polygons))
    for i in range(n_polygons):
        for j in range(i + 1, n_polygons):
            dist = gdf.geometry.iloc[i].distance(gdf.geometry.iloc[j])
            if dist <= d_max[i] and dist <= d_max[j]:
                G_prox.add_edge(i, j, weight=dist)

    # --- STEP 3 : Réduction au Minimal Spanning Tree (MST) ---
    MST = nx.Graph()
    MST.add_nodes_from(range(n_polygons))
    for cc in nx.connected_components(G_prox):
        subgraph = G_prox.subgraph(cc)
        if len(subgraph) > 1:
            mst_sub = nx.minimum_spanning_tree(subgraph, weight='weight')
            MST.add_edges_from(mst_sub.edges(data=True))

    # --- STEP 4 : Sélection Robuste des Graines (Seeds) ---
    # 4a. Graines de Cluster (Stricte condition de l'article : degré >= 4)
    cluster_seeds = [node for node, degree in MST.degree() if degree >= 4]
    
    # 4b. AJOUT : Graines d'Alignement (Pour capturer la chaîne linéaire du bas)
    # On cherche les sous-chemins linéaires significatifs (nœuds de degré 2 consécutifs)
    line_seeds = []
    # Création d'un sous-graphe uniquement avec les nœuds de chemin (degré 2) et d'extrémités (degré 1)
    path_nodes = [node for node, degree in MST.degree() if degree in [1, 2]]
    sub_paths = MST.subgraph(path_nodes)
    
    for cc in nx.connected_components(sub_paths):
        # Si la chaîne linéaire contient assez d'îles, on élit le milieu comme graine
        if len(cc) >= min_line_length:
            # Trouver un nœud central dans cette chaîne isolée
            cc_list = list(cc)
            degrees = [MST.degree(n) for n in cc_list]
            # On évite de prendre un nœud qui est en fait directement collé au grand cluster
            valid_nodes = [n for n in cc_list if MST.degree(n) == 2]
            if valid_nodes:
                # On prend le nœud au milieu de la chaîne
                line_seeds.append(valid_nodes[len(valid_nodes) // 2])
            elif cc_list:
                line_seeds.append(cc_list[len(cc_list) // 2])

    # Combinaison de toutes les graines uniques
    seeds = list(set(cluster_seeds + line_seeds))
    
    if not seeds:
        gdf['group_id'] = -1
        return gdf

    # --- STEP 5 : Voronoi regions of the seeds & traversal of the MST ---
    centroids = gdf.geometry.centroid
    seed_coords = np.array([[centroids.iloc[s].x, centroids.iloc[s].y] for s in seeds])
    spatial_tree = KDTree(seed_coords)
    
    island_coords = np.array([[c.x, c.y] for c in centroids])
    _, closest_seed_indices = spatial_tree.query(island_coords)
    voronoi_assignment = {i: seeds[closest_seed_indices[i]] for i in range(n_polygons)}

    initial_groups = {s: {s} for s in seeds}
    for s in seeds:
        queue = [s]
        visited = {s}
        while queue:
            current = queue.pop(0)
            for neighbor in MST.neighbors(current):
                if neighbor not in visited:
                    if voronoi_assignment[neighbor] == s:
                        visited.add(neighbor)
                        initial_groups[s].add(neighbor)
                        queue.append(neighbor)

    # --- STEP 6 : Extension aux îles restantes (un seul groupe limitrophe) ---
    assigned_islands = set().union(*initial_groups.values())
    unassigned_islands = set(range(n_polygons)) - assigned_islands
    
    # Répéter l'assignation tant qu'on peut assigner des membres clairs
    changed = True
    while changed:
        changed = False
        for island in list(unassigned_islands):
            neighbors = list(MST.neighbors(island))
            connected_to_groups = set()
            for n in neighbors:
                for s, group in initial_groups.items():
                    if n in group:
                        connected_to_groups.add(s)
            if len(connected_to_groups) == 1:
                s_target = list(connected_to_groups)[0]
                initial_groups[s_target].add(island)
                unassigned_islands.remove(island)
                changed = True

    # --- STEP 7 : Fusion des méso-structures trop proches ---
    G_seeds = nx.Graph()
    G_seeds.add_nodes_from(seeds)
    for i in range(len(seeds)):
        for j in range(i + 1, len(seeds)):
            s1, s2 = seeds[i], seeds[j]
            dist = gdf.geometry.iloc[s1].distance(gdf.geometry.iloc[s2])
            if dist <= d_max[s1] and dist <= d_max[s2]:
                G_seeds.add_edge(s1, s2)
                
    final_clusters = []
    processed_seeds = set()
    for s in seeds:
        if s not in processed_seeds:
            connected_seeds = nx.node_connected_component(G_seeds, s)
            merged_group = set()
            for cs in connected_seeds:
                merged_group.update(initial_groups[cs])
                processed_seeds.add(cs)
            final_clusters.append(merged_group)

    # --- STEP 8 : Final filling ---
    gdf['group_id'] = -1
    for cluster_idx, cluster in enumerate(final_clusters):
        for island in cluster:
            gdf.loc[island, 'group_id'] = cluster_idx
            
    return gdf