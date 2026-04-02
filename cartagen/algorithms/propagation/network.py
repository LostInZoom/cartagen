import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points
from collections import deque
from typing import Dict, Set, Tuple, List

def propagation_network(network, index, geometry, propagation_distance=None, damping=0.1, tolerance=1e-6, inplace=False):
    """
    Propagate a displacement along the network.

    This algorithm proposed by Lecordix *et al.* :footcite:p:`lecordix:1997` propagate
    a geometry change inside a network along the whole network. The extremities of
    the network is always kept.

    Parameters
    ----------
    network : GeoDataFrame of LineString or List[LineString]
        The road network in which to propagate the change.
    index : int
        The index of the modified geometry.
    geometry : LineString
        The new modified geometry.
    propagation_distance : float, optional
        The maximum distance to limit the propagation.
        Lower value will modify more heavily closer lines.
    damping : float, optional
        The damping factor.
    tolerance : float, optional
        The tolerance to consider touching roads.
    inplace : bool, optional
        Whether to modify the GeoDataFrame in place.

    Returns
    -------
    GeoDataFrame of LineString or List[LineString]

    See Also
    --------
    propagation_crow_flies :
        Propagate a displacement using the method "as the crow flies".
    propagation_network_batch :
        Propagate a displacement along the network by batch.

    References
    ----------
    .. footbibliography::
    """
    if not inplace:
        network = network.copy()

    is_list = isinstance(network, list)

    # 1. Calculate the initial displacement at the nodes of the modified line
    old_geom = network[index] if is_list else network.loc[index, 'geometry']
    
    def get_vec(p1, p2):
        return (p2.x - p1.x, p2.y - p1.y)

    disp_start = get_vec(Point(old_geom.coords[0]), Point(geometry.coords[0]))
    disp_end = get_vec(Point(old_geom.coords[-1]), Point(geometry.coords[-1]))

    # Update the central line
    if is_list:
        network[index] = geometry
    else:
        network.at[index, 'geometry'] = geometry

    # 2. Connectivity index and identification of true boundaries
    conn_map = {}
    
    # We iterate differently depending on whether it's a list or a DataFrame
    if is_list:
        for idx, geom in enumerate(network):
            coords = geom.coords
            for pos, pt_coords in [(0, coords[0]), (-1, coords[-1])]:
                pt_key = tuple(np.round(pt_coords, int(-np.log10(tolerance))))
                if pt_key not in conn_map:
                    conn_map[pt_key] = []
                conn_map[pt_key].append((idx, pos == 0))
    else:
        for idx, row in network.iterrows():
            coords = row.geometry.coords
            for pos, pt_coords in [(0, coords[0]), (-1, coords[-1])]:
                pt_key = tuple(np.round(pt_coords, int(-np.log10(tolerance))))
                if pt_key not in conn_map:
                    conn_map[pt_key] = []
                conn_map[pt_key].append((idx, pos == 0))

    # True network terminal points (degree 1)
    true_boundaries = {pt for pt, lines in conn_map.items() if len(lines) == 1}

    # 3. BFS with "fake boundaries" management
    queue = deque()
    visited = {index}

    # Initialization with direct neighbors
    for pos, vec in [(0, disp_start), (-1, disp_end)]:
        if np.linalg.norm(vec) < tolerance: continue
        
        pt_key = tuple(np.round(old_geom.coords[pos], int(-np.log10(tolerance))))
        for neighbor_idx, is_start in conn_map.get(pt_key, []):
            if neighbor_idx not in visited:
                # We start with a distance of 0.0 for direct neighbors
                queue.append((neighbor_idx, vec, is_start, 0.0))

    while queue:
        curr_idx, vec_anchor, is_start, dist = queue.popleft()
        if curr_idx in visited: continue
        
        line = network[curr_idx] if is_list else network.loc[curr_idx, 'geometry']
        L = line.length
        # The cumulative distance is calculated at the opposite end of the line
        dist_at_end = dist + L
        
        # --- BOUNDARY LOGIC ---
        opp_pos = -1 if is_start else 0
        opp_pt_key = tuple(np.round(line.coords[opp_pos], int(-np.log10(tolerance))))
        
        # A "fake boundary" is triggered if we exceed the maximum distance
        is_fake_boundary = propagation_distance is not None and dist_at_end > propagation_distance
        is_true_boundary = opp_pt_key in true_boundaries
        
        if is_fake_boundary or is_true_boundary:
            # We force the extremity to remain still (total damping)
            target_damping = 0.0
            should_continue = False 
        else:
            # Normal propagation with standard damping
            target_damping = damping
            should_continue = True

        # Apply the transformation
        coords = np.array(line.coords)
        factors = np.linspace(1.0, target_damping, len(coords)) if is_start else np.linspace(target_damping, 1.0, len(coords))
        
        new_coords = coords.copy()
        new_coords[:, 0] += vec_anchor[0] * factors
        new_coords[:, 1] += vec_anchor[1] * factors
        
        if is_list:
            network[curr_idx] = LineString(new_coords)
        else:
            network.at[curr_idx, 'geometry'] = LineString(new_coords)
            
        visited.add(curr_idx)

        # Propagate to neighbors only if a boundary (true or fake) has not been reached
        if should_continue:
            vec_at_end = (vec_anchor[0] * target_damping, vec_anchor[1] * target_damping)
            for next_idx, next_is_start in conn_map.get(opp_pt_key, []):
                if next_idx not in visited:
                    queue.append((next_idx, vec_at_end, next_is_start, dist_at_end))

    return network

def propagation_network_batch(network, modifications, propagation_distance=None, damping=0.1, tolerance=1e-6, inplace=False):
    """
    Propagate a displacement along the network by batch.

    This algorithm proposed by Lecordix *et al.* :footcite:p:`lecordix:1997` propagate
    a geometry change inside a network along the whole network. The extremities of
    the network is always kept.

    Parameters
    ----------
    network : GeoDataFrame of LineString
        The road network in which to propagate the change.
    modifications : List[Tuple[int, LineString]]
        A list of modifications inside the network. A tuple of
        two -> the index of the row to modify, the new geometry.
    propagation_distance : float, optional
        The maximum distance to limit the propagation.
        Lower value will modify more heavily closer lines.
    damping : float, optional
        The damping factor.
    tolerance : float, optional
        The tolerance to consider touching roads.
    inplace : bool, optional
        Whether to modify the GeoDataFrame in place.

    Returns
    -------
    GeoDataFrame of LineString

    See Also
    --------
    propagation_network :
        Propagate a displacement along the network.

    References
    ----------
    .. footbibliography::
    """
    result_network = network if inplace else network.copy()
    
    for i, (line_idx, new_geom) in enumerate(modifications):
        # On applique la propagation sur l'état actuel du GeoDataFrame
        result_network = propagation_network(result_network, line_idx, new_geom, propagation_distance, damping, tolerance, True)
        
    return result_network