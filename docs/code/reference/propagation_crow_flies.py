import geopandas as gpd
from shapely.geometry import LineString
import cartagen as c4
from matplotlib import pyplot as plt

crs_wgs84 = "EPSG:4326"

# A simple road segment (LineString) as initiator
init_initial = LineString([(0, 0), (10, 0), (20, 1)])
init_final = LineString([(0, 1), (10, 1.5), (20, 2)])

initiator_initial_gdf = gpd.GeoDataFrame(gpd.DataFrame([{'name': 'init', 'geometry': init_initial}]), crs=crs_wgs84)
initiator_final_gdf = gpd.GeoDataFrame(gpd.DataFrame([{'name': 'init', 'geometry': init_final}]), crs=crs_wgs84)

# Movable objects (e.g., a nearby boundary line)
movable_geom = LineString([(0, 3), (10, 3), (20, 4)])
movable_objects_gdf = gpd.GeoDataFrame(gpd.DataFrame([{'name': 'movable_line', 'geometry': movable_geom}]), crs=crs_wgs84)

# Frozen objects (None for this simple case)
frozen_objects_gdf = gpd.GeoDataFrame(geometry=[])

# --- Algorithm Inputs ---
# Distance of propagation (SizePZ in the paper) [cite: 132]
distance_of_propagation = 5.0 
# Interval for initial vector computation [cite: 111]
vector_interval_distance = 2.0 
vectors_gdf, max_dep = c4.__interpolate_displacement_vectors(init_initial, init_final, vector_interval_distance, crs=crs_wgs84)

# propagate the displacement to the movable objects
propagated_gdf = c4.compute_propagation_crow_flies(
        movable_objects_gdf, 
        initiator_initial_gdf, 
        initiator_final_gdf,
        frozen_objects_gdf,
        distance_of_propagation,
        vector_interval_distance)

# --- Visualization ---
base = movable_objects_gdf.plot(facecolor='none', edgecolor='blue')
propagated_gdf.plot(ax=base, facecolor='none', edgecolor='red')
initiator_final_gdf.plot(ax=base, facecolor='none', edgecolor='red')
initiator_initial_gdf.plot(ax=base, facecolor='none', edgecolor='blue')
for v_idx, vector_row in vectors_gdf.iterrows():
    vect_geom = vector_row['geometry']
    if vect_geom.length > 0:
        base.annotate('', xy=(vect_geom.coords[1][0], vect_geom.coords[1][1]), 
                        xytext=(vect_geom.coords[0][0], vect_geom.coords[0][1]),
                        arrowprops=dict(arrowstyle='->', color='orange'))
base.axes.get_xaxis().set_visible(False)
base.axes.get_yaxis().set_visible(False)        
plt.show()