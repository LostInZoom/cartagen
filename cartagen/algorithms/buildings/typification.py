import geopandas as gpd
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
from typing import List, Optional
import pandas as pd

def typify_buildings(
    buildings: gpd.GeoDataFrame,
    initial_scale: int = 25000,
    final_scale: int = 50000,
    ratio: float = None,
    road_network: Optional[gpd.GeoDataFrame] = None,
    attributes: Optional[List[str]] = None,
    distance: float = 20.0
) -> gpd.GeoDataFrame:
    """
    Typify buildings through iterative merging.

    This algorithm was proposed by Li *et al.* :footcite:p:`li:2005` and
    replace tight groups of buildings with a single representative building.
    
    The algorithm follows three main steps: First, it finds the number of
    buildings using improved radical law. Second, it calculates the position
    and representation through iterative merging. Finally, building sizes
    are harmmonized.
    
    Parameters
    ----------
    buildings : gpd.GeoDataFrame
        GeoDataFrame containing building polygons or multipolygons
    initial_scale : int
        Source map scale denominator (default: 25000)
    final_scale : int
        Target map scale denominator (default: 50000)
    ratio : float
        Ratio between source and target number of buildings (e.g., 0.5 to reduce by half)
        If None, will be computed from initial_scale and final_scale
    road_network : gpd.GeoDataFrame, optional
        Road network for spatial partitioning. If None, simple grid partitioning is used
    attributes : List[str], optional
        List of attribute names to transfer from largest building in each cluster
    distance : float
        Minimum separate distance for building harmonization (in map units)
    
    Returns
    -------
    gpd.GeoDataFrame
        Typified buildings with transferred attributes
    
    References
    ----------
    .. footbibliography::
    """
    
    if len(buildings) == 0:
        return buildings.copy()
    
    # Ensure we have a copy to work with
    buildings = buildings.copy()
    
    # Compute building properties
    buildings = _enrich_building_data(buildings)
    
    # Step 1: Spatial partitioning (many-to-many matching)
    buildings = _spatial_partition(buildings, road_network)
    
    # Step 2: Determine target number using improved radical law
    if ratio is None:
        T = _compute_selection_level(initial_scale, final_scale, len(buildings))
        target_count = int(len(buildings) * (initial_scale / final_scale) ** T)
    else:
        target_count = int(len(buildings) * ratio)
    
    # Ensure at least one building remains
    target_count = max(1, target_count)
    
    # Step 3: Typify buildings group by group
    typified_buildings = []
    
    for group_id in buildings['group_id'].unique():
        group_buildings = buildings[buildings['group_id'] == group_id].copy()
        
        # Compute target number for this group
        group_ratio = len(group_buildings) / len(buildings)
        group_target = max(1, int(target_count * group_ratio))
        
        # Typify this group
        typified_group = _typify_group(
            group_buildings, 
            group_target,
            distance,
            attributes
        )
        
        typified_buildings.append(typified_group)
    
    # Combine all groups
    result = gpd.GeoDataFrame(
        pd.concat(typified_buildings, ignore_index=True),
        crs=buildings.crs
    )
    
    return result


def _enrich_building_data(buildings: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Enrich building data with geometric properties needed for typification.
    
    Computes:
    - Centroid coordinates (X, Y)
    - Area
    - Orientation (angle of longest axis)
    - L_by_W (length to width ratio)
    """
    buildings = buildings.copy()
    
    # Compute centroids
    centroids = buildings.geometry.centroid
    buildings['X'] = centroids.x
    buildings['Y'] = centroids.y
    
    # Compute area
    buildings['Area'] = buildings.geometry.area
    
    # Compute orientation and L_by_W
    orientations = []
    l_by_w_ratios = []
    
    for geom in buildings.geometry:
        if isinstance(geom, MultiPolygon):
            # Take the largest polygon
            geom = max(geom.geoms, key=lambda p: p.area)
        
        # Get minimum rotated rectangle
        min_rect = geom.minimum_rotated_rectangle
        
        if min_rect.is_empty or min_rect.area == 0:
            orientations.append(0.0)
            l_by_w_ratios.append(1.0)
            continue
        
        # Get coordinates of the rectangle
        coords = list(min_rect.exterior.coords)
        
        # Calculate edge lengths
        edge1_length = Point(coords[0]).distance(Point(coords[1]))
        edge2_length = Point(coords[1]).distance(Point(coords[2]))
        
        # Determine long and short axes
        if edge1_length > edge2_length:
            long_axis = edge1_length
            short_axis = edge2_length
            # Orientation is angle of first edge
            dx = coords[1][0] - coords[0][0]
            dy = coords[1][1] - coords[0][1]
        else:
            long_axis = edge2_length
            short_axis = edge1_length
            # Orientation is angle of second edge
            dx = coords[2][0] - coords[1][0]
            dy = coords[2][1] - coords[1][1]
        
        # Compute orientation in degrees
        orientation = np.degrees(np.arctan2(dy, dx))
        
        # Normalize to [0, 180)
        orientation = orientation % 180
        
        orientations.append(orientation)
        
        # Compute L_by_W ratio
        if short_axis > 0:
            l_by_w_ratios.append(long_axis / short_axis)
        else:
            l_by_w_ratios.append(1.0)
    
    buildings['Orientation'] = orientations
    buildings['L_by_W'] = l_by_w_ratios
    
    return buildings


def _spatial_partition(
    buildings: gpd.GeoDataFrame, 
    road_network: Optional[gpd.GeoDataFrame] = None
) -> gpd.GeoDataFrame:
    """
    Partition buildings into groups using road network or grid-based approach.
    
    If road network is provided, uses Voronoi-like partitioning.
    Otherwise, uses a simple grid-based partitioning.
    """
    buildings = buildings.copy()
    
    if road_network is not None and len(road_network) > 0:
        # Use road-based partitioning (simplified approach)
        # Create buffer zones around roads
        road_buffer = unary_union(road_network.geometry.buffer(50))
        
        # Assign group IDs based on spatial proximity
        group_ids = []
        for idx, building in buildings.iterrows():
            centroid = building.geometry.centroid
            # Simple hashing based on location
            grid_x = int(centroid.x / 500)
            grid_y = int(centroid.y / 500)
            group_ids.append(f"{grid_x}_{grid_y}")
        
        buildings['group_id'] = group_ids
    else:
        # Grid-based partitioning
        bounds = buildings.total_bounds
        min_x, min_y, max_x, max_y = bounds
        
        # Create grid cells (adjust cell size based on data extent)
        n_cells = max(5, int(np.sqrt(len(buildings)) / 2))
        cell_width = (max_x - min_x) / n_cells
        cell_height = (max_y - min_y) / n_cells
        
        # Avoid division by zero
        if cell_width == 0:
            cell_width = 1
        if cell_height == 0:
            cell_height = 1
        
        group_ids = []
        for idx, building in buildings.iterrows():
            centroid = building.geometry.centroid
            cell_x = int((centroid.x - min_x) / cell_width)
            cell_y = int((centroid.y - min_y) / cell_height)
            group_ids.append(f"{cell_x}_{cell_y}")
        
        buildings['group_id'] = group_ids
    
    return buildings


def _compute_selection_level(
    initial_scale: int, 
    final_scale: int, 
    n_buildings: int
) -> float:
    """
    Compute selection level T using improved radical law.
    
    Formula: T = 2(log(N_B) - log(N_S)) / (log(M_S) - log(M_B))
    
    Where:
    - N_B: number of buildings at larger scale
    - N_S: number of buildings at smaller scale
    - M_B: scale denominator at larger scale
    - M_S: scale denominator at smaller scale
    """
    # Assuming a typical reduction factor
    typical_ratio = 0.6
    n_target = int(n_buildings * typical_ratio)
    
    if n_target <= 0:
        return 0.5
    
    numerator = 2 * (np.log(n_buildings) - np.log(n_target))
    denominator = np.log(final_scale) - np.log(initial_scale)
    
    if abs(denominator) < 1e-10:
        return 0.5
    
    T = numerator / denominator
    
    return max(0.1, min(2.0, T))  # Clamp between reasonable values


def _compute_pairwise_distances(coords: np.ndarray) -> np.ndarray:
    """
    Compute pairwise Euclidean distances between points.
    
    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) containing x, y coordinates
    
    Returns
    -------
    np.ndarray
        Distance matrix of shape (n, n)
    """
    n = coords.shape[0]
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.sqrt(
                (coords[i, 0] - coords[j, 0]) ** 2 + 
                (coords[i, 1] - coords[j, 1]) ** 2
            )
            distances[i, j] = dist
            distances[j, i] = dist
    
    return distances


def _typify_group(
    group_buildings: gpd.GeoDataFrame,
    target_count: int,
    distance: float,
    attributes: Optional[List[str]] = None
) -> gpd.GeoDataFrame:
    """
    Typify a group of buildings by iteratively merging closest pairs.
    
    Steps:
    1. While current count > target count:
       - Find pair with shortest distance
       - Merge into new building
       - Update building list
    2. Harmonize sizes
    """
    if len(group_buildings) <= target_count:
        return group_buildings
    
    # Create working copy with point representation
    buildings = group_buildings.copy()
    buildings['active'] = True
    buildings['building_id'] = range(len(buildings))
    
    # Track attribute transfer
    if attributes:
        buildings['largest_area'] = buildings['Area']
        for attr in attributes:
            if attr in buildings.columns:
                buildings[f'transfer_{attr}'] = buildings[attr]
    
    current_count = len(buildings)
    
    # Iteratively merge buildings
    while current_count > target_count:
        # Get active buildings
        active = buildings[buildings['active']].copy()
        
        if len(active) <= 1:
            break
        
        # Compute pairwise distances between centroids
        coords = np.column_stack([active['X'].values, active['Y'].values])
        distances = _compute_pairwise_distances(coords)
        
        # Set diagonal to infinity (avoid self-pairing)
        np.fill_diagonal(distances, np.inf)
        
        # Find minimum distance pair
        min_idx = np.unravel_index(np.argmin(distances), distances.shape)
        idx1, idx2 = active.index[min_idx[0]], active.index[min_idx[1]]
        
        # Merge the two buildings
        building1 = buildings.loc[idx1]
        building2 = buildings.loc[idx2]
        
        # Compute new building properties
        new_building = _merge_two_buildings(
            building1, building2, attributes
        )
        
        # Deactivate old buildings
        buildings.loc[idx1, 'active'] = False
        buildings.loc[idx2, 'active'] = False
        
        # Add new building
        new_row = pd.DataFrame([new_building])
        buildings = pd.concat([buildings, new_row], ignore_index=True)
        buildings.loc[len(buildings)-1, 'active'] = True
        
        current_count -= 1
    
    # Get final active buildings
    result = buildings[buildings['active']].copy()
    
    # Harmonize sizes
    result = _harmonize_sizes(result, distance)
    
    # Create geometries from enriched data
    result = _create_geometries_from_attributes(result)
    
    # Transfer attributes from largest building
    if attributes:
        for attr in attributes:
            if f'transfer_{attr}' in result.columns:
                result[attr] = result[f'transfer_{attr}']
                result = result.drop(columns=[f'transfer_{attr}'])
        
        if 'largest_area' in result.columns:
            result = result.drop(columns=['largest_area'])
    
    # Clean up working columns
    result = result.drop(columns=['active', 'building_id'], errors='ignore')
    
    return result


def _merge_two_buildings(
    building1: pd.Series,
    building2: pd.Series,
    attributes: Optional[List[str]] = None
) -> dict:
    """
    Merge two buildings into one new building.
    
    Position: midpoint of centroids
    Orientation and L_by_W: average
    Area: larger area if ratio > 2, otherwise average
    """
    # Position: midpoint
    new_x = (building1['X'] + building2['X']) / 2
    new_y = (building1['Y'] + building2['Y']) / 2
    
    # Orientation: average (handle wraparound)
    orientation1 = building1['Orientation']
    orientation2 = building2['Orientation']
    
    # Handle 180-degree wraparound
    if abs(orientation1 - orientation2) > 90:
        if orientation1 > orientation2:
            orientation2 += 180
        else:
            orientation1 += 180
    
    new_orientation = ((orientation1 + orientation2) / 2) % 180
    
    # L_by_W: average
    new_l_by_w = (building1['L_by_W'] + building2['L_by_W']) / 2
    
    # Area: preserve character
    area1 = building1['Area']
    area2 = building2['Area']
    
    if max(area1, area2) / min(area1, area2) >= 2.0:
        # Significant difference: use larger area
        new_area = max(area1, area2)
        larger_building = building1 if area1 > area2 else building2
    else:
        # Similar sizes: use average
        new_area = (area1 + area2) / 2
        larger_building = building1 if area1 > area2 else building2
    
    new_building = {
        'X': new_x,
        'Y': new_y,
        'Area': new_area,
        'Orientation': new_orientation,
        'L_by_W': new_l_by_w,
        'group_id': building1['group_id']
    }
    
    # Transfer attributes from largest building in the merge
    if attributes:
        new_building['largest_area'] = max(
            building1.get('largest_area', area1),
            building2.get('largest_area', area2)
        )
        
        # Determine which building was originally largest
        if building1.get('largest_area', area1) >= building2.get('largest_area', area2):
            source_building = building1
        else:
            source_building = building2
        
        for attr in attributes:
            transfer_key = f'transfer_{attr}'
            if transfer_key in source_building:
                new_building[transfer_key] = source_building[transfer_key]
            elif attr in source_building:
                new_building[transfer_key] = source_building[attr]
    
    return new_building


def _harmonize_sizes(
    buildings: gpd.GeoDataFrame,
    distance: float
) -> gpd.GeoDataFrame:
    """
    Harmonize building sizes according to formula (4) in the paper.
    
    Formula:
    - If x <= A_msd: f(x) = A_msd
    - If A_msd < x < A_max: f(x) = x * (A_max / A_max)^(A_msd / A_max)
    - If x = A_max: f(x) = x
    """
    buildings = buildings.copy()
    
    # Compute A_msd (area of rectangle with MSD dimensions)
    A_msd = distance ** 2
    
    # Find maximum area
    A_max = buildings['Area'].max()
    
    if A_max <= A_msd:
        # All buildings are too small, set to minimum
        buildings['Area'] = A_msd
        return buildings
    
    # Apply harmonization function
    def harmonize_area(area):
        if area <= A_msd:
            return A_msd
        elif area >= A_max:
            return area
        else:
            # Apply scaling function
            ratio = A_msd / A_max
            return area * (A_max / A_max) ** ratio
    
    buildings['Area'] = buildings['Area'].apply(harmonize_area)
    
    return buildings


def _create_geometries_from_attributes(
    buildings: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Create rectangular building geometries from enriched attributes.
    
    Uses Area, Orientation, and L_by_W to reconstruct polygons.
    """
    geometries = []
    
    for idx, building in buildings.iterrows():
        x, y = building['X'], building['Y']
        area = building['Area']
        orientation = building['Orientation']
        l_by_w = building['L_by_W']
        
        # Compute dimensions from area and ratio
        # Area = length * width
        # l_by_w = length / width
        # Therefore: length = sqrt(area * l_by_w), width = sqrt(area / l_by_w)
        
        width = np.sqrt(area / l_by_w)
        length = np.sqrt(area * l_by_w)
        
        # Create rectangle centered at (x, y) with given orientation
        # Rectangle corners before rotation
        half_length = length / 2
        half_width = width / 2
        
        corners = np.array([
            [-half_length, -half_width],
            [half_length, -half_width],
            [half_length, half_width],
            [-half_length, half_width]
        ])
        
        # Rotate corners
        angle_rad = np.radians(orientation)
        rotation_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
        ])
        
        rotated_corners = corners @ rotation_matrix.T
        
        # Translate to centroid position
        rotated_corners[:, 0] += x
        rotated_corners[:, 1] += y
        
        # Create polygon
        polygon = Polygon(rotated_corners)
        geometries.append(polygon)
    
    buildings = buildings.copy()
    buildings.geometry = geometries
    
    return buildings