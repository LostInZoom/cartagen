import geopandas as gpd
from shapely.geometry import Point, LineString
import numpy as np
import pandas as pd

# --- Constants & Parameters ---
# The paper uses a cushioning parameter β (beta) in [0;1].
# SizePZ = max_dep / beta + max_dep (Size of Propagation Zone) [cite: 132]
# SizeIZ = max_dep / beta (Radius of Influence Zone) [cite: 150]
# Let's use SizePZ as the direct input and calculate beta from it.
# SizeIZ is derived from SizePZ for a given max_dep: SizeIZ = SizePZ - max_dep.

def calculate_size_iz(max_dep, propagation_distance):
    """Calculate SizeIZ (Radius of Influence Zone) based on max_dep and SizePZ (propagation_distance)."""
    # propagation_distance is SizePZ. SizeIZ = SizePZ - max_dep
    return propagation_distance - max_dep

def get_vertices(gdf):
    """Extract all unique vertices from a GeoDataFrame of (multi)lines/polygons."""
    vertices = []
    
    # Placeholder for actual complex geometry iteration
    # In a real implementation, you would iterate over all geometries,
    # then over their components (exteriors, interiors, segments), and extract all unique points.
    
    for geom in gdf.geometry:
        if geom.geom_type in ['LineString', 'Point']:
            if geom.geom_type == 'Point':
                vertices.append(geom)
            else:
                vertices.extend(list(geom.coords))
        elif geom.geom_type in ['MultiLineString', 'MultiPolygon', 'Polygon']:
            # This is where a lot of complexity lives for multi-geometries
            pass # Skipping for conceptual implementation
            
    # Convert to a unique list of Shapely Points with an ID/index
    # For this conceptual example, we'll return a simple structure.
    
    # Returning a GeoDataFrame of all vertices is necessary for join/spatial indexing
    # We'll use a simplified list of coordinates for this conceptual part.
    all_coords = []
    # Simplified extraction for single LineStrings (common for roads)
    for i, row in gdf.iterrows():
        if row.geometry.geom_type == 'LineString':
             all_coords.extend([(row.geometry.coords[j], i, j) for j in range(len(row.geometry.coords))]) # (coord, obj_idx, vertex_idx)
    
    # Create a DataFrame for easier processing
    vertices_df = pd.DataFrame(all_coords, columns=['coords', 'obj_idx', 'vertex_idx'])
    vertices_df[['x', 'y']] = vertices_df['coords'].apply(lambda c: pd.Series([c[0], c[1]]))
    vertices_gdf = gpd.GeoDataFrame(
        vertices_df, 
        geometry=gpd.points_from_xy(vertices_df.x, vertices_df.y),
        crs=gdf.crs
    )
    return vertices_gdf

# --- Core Geometric Functions (Highly Simplified Placeholders) ---

def interpolate_displacement_vectors(initial_geom, final_geom, interval_distance, crs):
    """
    Computes initial displacement vectors using interpolation at regular intervals.
    Origin point on initial_geom, extremity on final_geom (same curvilinear abscissa ratio).
    
    Returns a GeoDataFrame of LineStrings (vectors) and a Series of vector lengths.
    """
    # This function is non-trivial and requires geometry-specific interpolation logic.
    # Placeholder: Assuming simple LineString initiators for demonstration.
    
    if initial_geom.geom_type != 'LineString' or final_geom.geom_type != 'LineString':
        # Must handle complex types (e.g., MultiLineString) in a real implementation
        return gpd.GeoDataFrame(), 0.0

    # Simplified mock for a few vectors
    vectors_data = []
    num_steps = int(initial_geom.length / interval_distance)
    
    for i in range(num_steps + 1):
        ratio = i / num_steps
        
        # Homolog points using same ratio of curvilinear abscissa [cite: 111]
        p_initial = initial_geom.interpolate(initial_geom.length * ratio)
        p_final = final_geom.interpolate(final_geom.length * ratio)
        
        vector = LineString([p_initial, p_final])
        dep_ini = vector.length
        
        vectors_data.append({
            'geometry': vector, 
            'origin': p_initial, 
            'dep_ini': dep_ini,
            'direction': np.array([p_final.x - p_initial.x, p_final.y - p_initial.y]) / dep_ini if dep_ini > 0 else np.array([0.0, 0.0])
        })
    
    vectors_gdf = gpd.GeoDataFrame(vectors_data, crs=crs)
    max_dep = vectors_gdf['dep_ini'].max() if not vectors_gdf.empty else 0.0
    return vectors_gdf, max_dep

def visibility_condition_check(vertex, vector_origin, initiator_geom):
    """
    Checks if the vector_origin is visible from the vertex, i.e., not hidden by any initiator object. 
    This is a ray-tracing-like operation, highly complex for general geometries.
    Placeholder: always True for conceptual demo.
    """
    return True

# --- Main Propagation Logic ---

def compute_propagation_crow_flies(
    movable_objects_gdf: gpd.GeoDataFrame, 
    initiator_initial_gdf: gpd.GeoDataFrame, 
    initiator_final_gdf: gpd.GeoDataFrame,
    frozen_objects_gdf: gpd.GeoDataFrame,
    propagation_distance: float,
    interval_distance=2.0) -> gpd.GeoDataFrame:
    """
    Implementation of the displacement propagation algorithm "as the crow flies", by Legrand *et al.* :footcite:p:`legrand:2005`.
    This function propagates the displacement defined by the movement of initiator objects to nearby movable objects,
    while respecting frozen objects that should not be moved.

    Be careful, the initiators should be simple LineStrings.

    Parameters
    ----------
        movable_objects_gdf : The geopandas Geodataframe containing the movable objects.
        initiator_initial_gdf : The geopandas Geodataframe containing the initial geometry of the initiators.
        initiator_final_gdf : The geopandas Geodataframe containing the final geometry of the initiators.
        frozen_objects_gdf : The geopandas Geodataframe containing the frozen objects, on which no propagation is carried out.
        propagation_distance : SizePZ in the paper by Legrand et al., i.e. the distance around the initiator to propagate displacement.
        interval_distance : For initial vector calculation.
    
    Returns
    -------
    A geopandas GeoDataFrame with the same structure as `movable_objects_gdf`, but with updated geometries after propagation.

    References
    ----------
    .. footbibliography::
    """
    # 1. Compute Initial Displacement Vectors and Max Displacement
    # Assuming one initiator for simplicity; in reality, iterate over all.
    initiator_initial = initiator_initial_gdf.iloc[0].geometry
    initiator_final = initiator_final_gdf.iloc[0].geometry
    
    initial_vectors_gdf, max_dep = interpolate_displacement_vectors(
        initiator_initial, initiator_final, interval_distance, crs=movable_objects_gdf.crs
    )
    
    if max_dep == 0.0:
        return movable_objects_gdf.copy() # No displacement to propagate

    # Calculate Influence Zone Radius (SizeIZ) [cite: 150]
    # Note: SizePZ = propagation_distance [cite: 132]
    size_iz = calculate_size_iz(max_dep, propagation_distance)
    
    # 2. Compute Propagation Zone (Buffer)
    propagation_zone = initiator_initial.buffer(propagation_distance) # [cite: 127]
    
    # 3. Identify Movable Vertices within the Propagation Zone
    movable_vertices_gdf = get_vertices(movable_objects_gdf)
    
    # Filter vertices *outside* the propagation zone (no displacement)
    vertices_in_zone = movable_vertices_gdf.clip(propagation_zone)
    
    # Exclude vertices of frozen objects (null displacement) 
    # This simplified version just uses the *geometry* of frozen objects to filter
    if not frozen_objects_gdf.empty:
         # Use spatial join or overlay difference to exclude
         # Note: A proper check should see if the vertex *belongs* to a frozen object edge/node.
         vertices_in_zone = vertices_in_zone[~vertices_in_zone.geometry.intersects(frozen_objects_gdf.geometry.unary_union)]

    # 4. Compute Propagation for each Vertex
    propagated_vertices_data = []

    for idx, vertex_row in vertices_in_zone.iterrows():
        vertex = vertex_row.geometry
        
        # 4.1. Selection of Influencing Initial Displacement Vectors
        
        # Distance Condition: Origin of vector inside circle of radius SizeIZ around vertex [cite: 153]
        # Equivalent: Vertex inside Influence Zone (SizeIZ radius) of the vector's origin [cite: 150]
        influencing_vectors = initial_vectors_gdf[initial_vectors_gdf.geometry.apply(lambda v: vertex.distance(Point(v.coords[0])) <= size_iz)]

        final_displacement_vector = np.array([0.0, 0.0])
        total_weight = 0.0
        
        for v_idx, vector_row in influencing_vectors.iterrows():
            # 4.1. Visibility Condition (Skipping for conceptual demo, see placeholder)
            if not visibility_condition_check(vertex, vector_row.origin, initiator_initial):
                continue

            dist_i_vertex = vertex.distance(vector_row.origin)
            dep_ini = vector_row['dep_ini']
            
            # Distance must be > 0 to calculate weight; if dist_i_vertex=0, the displacement is dep_ini.
            if dist_i_vertex == 0.0:
                # Vertex is on the initiator's initial position (an initial vector origin)
                # It should receive the full initial displacement.
                # Since these points should theoretically be handled as part of the initiator's *new* geometry,
                # we'll use a very small epsilon for weight calculation to avoid ZeroDivisionError, or simply
                # use full displacement and continue, or handle it as a frozen point.
                # Here, we'll treat it as a point that receives the full displacement.
                vector_propagation = dep_ini * vector_row['direction']
                weight = 1.0 / (1e-6)**2 # Max weight

            else:
                # 4.2. Compute Displacement Vector due to one Initial Vector (Quadratic Cushioning)
                # cushioning(i, vertex) = dist²(i, vertex) / SizeIZ² [cite: 175]
                cushioning_coeff = (dist_i_vertex / size_iz)**2 
                
                # dep_i(vertex) = dep_ini * (1 - cushioning(i, vertex)) [cite: 168, 177]
                dep_vertex_length = dep_ini * (1 - cushioning_coeff)
                
                vector_propagation = dep_vertex_length * vector_row['direction']
                
                # 4.3. Weight for Aggregation (Weighted Average)
                # weight(i, vertex) = 1/dist²(i, vertex) 
                weight = 1.0 / (dist_i_vertex**2) 
            
            final_displacement_vector += vector_propagation * weight
            total_weight += weight

        # 4.4. Aggregation (Weighted Average)
        if total_weight > 0:
            aggregated_displacement = final_displacement_vector / total_weight
        else:
            # Vertex not influenced by any vector (e.g., hidden or too far after filtering)
            aggregated_displacement = np.array([0.0, 0.0])
        
        # New Position
        new_x = vertex.x + aggregated_displacement[0]
        new_y = vertex.y + aggregated_displacement[1]
        
        propagated_vertices_data.append({
            'obj_idx': vertex_row['obj_idx'],
            'vertex_idx': vertex_row['vertex_idx'],
            'original_geom': vertex_row.geometry,
            'new_geom': Point(new_x, new_y)
        })

    # 5. Reconstruct New Geometries
    propagated_vertices_df = pd.DataFrame(propagated_vertices_data)
    
    # Merge the new positions back into the original movable_objects structure
    new_geometries = movable_objects_gdf.copy()
    
    # Placeholder for geometry reconstruction
    # This involves:
    # 1. Iterating through the original geometries of 'movable_objects_gdf'.
    # 2. For each geometry, collecting its vertices' new positions from 'propagated_vertices_df'.
    # 3. Creating a new Shapely geometry (LineString, Polygon, etc.) from these new points.
    
    for obj_idx, obj_group in propagated_vertices_df.groupby('obj_idx'):
        original_line = movable_objects_gdf.loc[obj_idx].geometry
        new_coords = []
        
        # Reconstruct coordinates in original order
        sorted_group = obj_group.sort_values('vertex_idx')
        for _, row in sorted_group.iterrows():
            new_coords.append((row['new_geom'].x, row['new_geom'].y))
            
        # Handle vertices outside the propagation zone (they keep original position)
        if original_line.geom_type == 'LineString':
            original_coords = list(original_line.coords)
    
            
            # Assuming all vertices of the line are in the propagated set for the demo.
            if len(new_coords) == len(original_coords):
                new_geom = LineString(new_coords)
                new_geometries.at[obj_idx, 'geometry'] = new_geom
            else:
                final_coords = []
                for i, orig_coord in enumerate(original_coords):
                    if i in sorted_group['vertex_idx'].values:
                        new_point = sorted_group[sorted_group['vertex_idx'] == i]['new_geom'].values[0]
                        final_coords.append((new_point.x, new_point.y))
                    else:
                        final_coords.append(orig_coord)
                new_geom = LineString(final_coords)
                new_geometries.at[obj_idx, 'geometry'] = new_geom
    
    return new_geometries