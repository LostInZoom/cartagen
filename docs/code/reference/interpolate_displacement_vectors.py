from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
import shapely
from shapely import LineString
from shapely.wkt import loads
import cartagen as c4

# A simple road segment (LineString) as initiator
initial = LineString([(0, 0), (10, 0), (20, 1)])
final = LineString([(0, 1), (10, 1.5), (20, 2)])

crs = 3857
initiator_initial_gdf = gpd.GeoDataFrame([{'name': 'init', 'geometry': initial}], crs=crs)
initiator_final_gdf = gpd.GeoDataFrame([{'name': 'init', 'geometry': final}], crs=crs)

# Interval for initial vector computation [cite: 111]
interval = 2.0 
vectors_gdf, max_dep = c4.interpolate_displacement_vectors(initial, final, interval, crs=crs)

fig = plt.figure(1, (10, 4))

#############################################################

sub1 = fig.add_subplot(211)
sub1.set_aspect('equal')
sub1.set_title('interval=2.0', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

for line in initiator_initial_gdf.geometry:
    path1 = Path(numpy.asarray(line.coords)[:, :2])
    sub1.add_patch(PathPatch(path1, facecolor="none", edgecolor='blue', linewidth=1))

for line in initiator_final_gdf.geometry:
    path1 = Path(numpy.asarray(line.coords)[:, :2])
    sub1.add_patch(PathPatch(path1, facecolor="none", edgecolor='red', linewidth=1))

for v_idx, vector_row in vectors_gdf.iterrows():
    vect_geom = vector_row['geometry']
    if vect_geom.length > 0:
        sub1.annotate('', xy=(vect_geom.coords[1][0], vect_geom.coords[1][1]), 
                        xytext=(vect_geom.coords[0][0], vect_geom.coords[0][1]),
                        arrowprops=dict(arrowstyle='->', color='red'))

sub1.autoscale_view()

sub2 = fig.add_subplot(212)
sub2.set_aspect('equal')
sub2.set_title('interval=5.0', pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

# Interval for initial vector computation [cite: 111]
interval = 5 
vectors_gdf, max_dep = c4.interpolate_displacement_vectors(initial, final, interval, crs=crs)

for line in initiator_initial_gdf.geometry:
    path1 = Path(numpy.asarray(line.coords)[:, :2])
    sub2.add_patch(PathPatch(path1, facecolor="none", edgecolor='blue', linewidth=1))

for line in initiator_final_gdf.geometry:
    path1 = Path(numpy.asarray(line.coords)[:, :2])
    sub2.add_patch(PathPatch(path1, facecolor="none", edgecolor='red', linewidth=1))

for v_idx, vector_row in vectors_gdf.iterrows():
    vect_geom = vector_row['geometry']
    if vect_geom.length > 0:
        sub2.annotate('', xy=(vect_geom.coords[1][0], vect_geom.coords[1][1]), 
                        xytext=(vect_geom.coords[0][0], vect_geom.coords[0][1]),
                        arrowprops=dict(arrowstyle='->', color='red'))

sub2.autoscale_view()

plt.show()