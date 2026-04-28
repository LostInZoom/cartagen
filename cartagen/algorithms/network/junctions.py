import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
from shapely.geometry import LineString
from scipy.cluster.hierarchy import linkage, fcluster

def collapse_junctions(network, threshold):
    """
    Collapse road junctions using hierarchical clustering.

    This algorithm, proposed by Mackaness and Mackechnie :footcite:p:`mackaness:1999`,
    simplifies road intersections by grouping nearby nodes 
    using a single-linkage clustering method. For each cluster, topological 
    connectivity is analyzed to identify independent road groups (connected 
    components). Each component is then collapsed toward its specific centroid, 
    preserving the network's functional topology while removing micro-segments.

    Parameters
    ----------
    network : GeoDataFrame of LineString
        The road network to simplify.
    threshold : float
        The maximum distance (in the network's coordinate system units) 
        to group nodes together into a single cluster.

    Returns
    -------
    GeoDataFrame of LineString
        The simplified road network with collapsed junctions and updated geometries.

    Notes
    -----
    The algorithm first performs a spatial clustering to identify potential 
    junctions. Within each spatial cluster, it builds a graph to ensure that 
    only nodes that are physically connected are collapsed together, avoiding 
    artificial connections between nearby but disconnected layers (e.g., overpasses).

    References
    ----------
    .. footbibliography::
    """

    gdf = network.copy()
    if gdf.index.duplicated().any():
        gdf = gdf.reset_index(drop=True)

    def normalize(coords):
        return (round(float(coords[0]), 8), round(float(coords[1]), 8))

    # --- 1. Capture de la topologie initiale (Anti-KeyError) ---
    node_to_edges = {} 
    edge_to_nodes = {} # Mémorise les IDs de nœuds pour chaque ligne
    unique_coords = []
    coord_to_id = {}

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty: continue
        
        # On extrait et normalise une fois pour toutes
        s_coords = normalize(geom.coords[0])
        e_coords = normalize(geom.coords[-1])
        
        # Gestion des IDs de nœuds
        for pt in [s_coords, e_coords]:
            if pt not in coord_to_id:
                coord_to_id[pt] = len(unique_coords)
                unique_coords.append(pt)
            
            node_id = coord_to_id[pt]
            if node_id not in node_to_edges:
                node_to_edges[node_id] = []
            node_to_edges[node_id].append(idx)
        
        # On stocke la relation edge -> (node_id_start, node_id_end)
        edge_to_nodes[idx] = (coord_to_id[s_coords], coord_to_id[e_coords])

    nodes_df = pd.DataFrame({'coords': unique_coords, 'node_id': range(len(unique_coords))})

    # --- 2. Clustering Spatial ---
    coords_array = np.array([list(c) for c in unique_coords])
    Z = linkage(coords_array, method='single', metric='euclidean')
    nodes_df['cluster_id'] = fcluster(Z, t=threshold, criterion='distance')

    edges_to_delete = set()

    # --- 3. Traitement par Cluster et Sous-Graphes ---
    for cluster_id, cluster_group in nodes_df.groupby('cluster_id'):
        if len(cluster_group) <= 1:
            continue
            
        G_cluster = nx.Graph()
        cluster_node_ids = set(cluster_group['node_id'])
        
        # On construit le graphe à partir de notre mémoire edge_to_nodes
        for n_id in cluster_node_ids:
            G_cluster.add_node(n_id)
            for edge_idx in node_to_edges.get(n_id, []):
                s_id, e_id = edge_to_nodes[edge_idx]
                if s_id in cluster_node_ids and e_id in cluster_node_ids:
                    G_cluster.add_edge(s_id, e_id, edge_index=edge_idx)

        for component in nx.connected_components(G_cluster):
            if len(component) <= 1:
                continue
            
            # Centroïde
            comp_points = np.array([list(unique_coords[nid]) for nid in component])
            centroid = tuple(np.mean(comp_points, axis=0))
            centroid_norm = (float(centroid[0]), float(centroid[1]))
            
            for n_id in component:
                # On utilise node_to_edges (fixe) au lieu de recalculer sur la géométrie
                for edge_idx in node_to_edges.get(n_id, []):
                    s_id, e_id = edge_to_nodes[edge_idx]
                    other_id = e_id if s_id == n_id else s_id
                    
                    if other_id in component:
                        edges_to_delete.add(edge_idx)
                    else:
                        # On ne modifie la géométrie que pour les routes qui sortent
                        line = gdf.loc[edge_idx].geometry
                        new_coords = list(line.coords)
                        # On compare les IDs (plus robuste que les coords)
                        if s_id == n_id:
                            new_coords[0] = centroid_norm
                        else:
                            new_coords[-1] = centroid_norm
                        
                        gdf.at[edge_idx, 'geometry'] = LineString(new_coords)

    return gdf.drop(list(edges_to_delete))