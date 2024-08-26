from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
import shapely
from shapely.wkt import loads
import cartagen as c4
from shapely import LineString, Polygon, Point
from cartagen.utils.debug import plot_debug, geojson_debug, geojson_to_variable

line = loads('LineString (-187714.88334235321963206 5374960.62342192325741053, -187618.28458588907960802 5375091.72173426765948534, -187494.08618472088710405 5374864.02466545905917883, -187354.01798784791026264 5375088.96176979690790176, -187229.12959556211717427 5374852.98480757791548967, -187120.80099009876721539 5374972.35327092278748751)')

fig = plt.figure(1, (12, 16))

sub1 = fig.add_subplot(211)
sub1.set_title('keep_vertices=False, original vertices in blue, resampled line vertices in red.', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(212)
sub2.set_title('keep_vertices=True, original vertices in blue, resampled line vertices in red.', pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

resampled1 = c4.resample_line(line, 50, keep_vertices=False)
resampled2 = c4.resample_line(line, 50, keep_vertices=True)

coords1 = list(line.coords)
coords2 = list(resampled1.coords)
coords3 = list(resampled2.coords)

path1 = Path(numpy.asarray(line.coords)[:, :2])
sub1.add_patch(PathPatch(path1, facecolor="none", edgecolor='gray', linewidth=1))

path2 = Path(numpy.asarray(resampled1.coords)[:, :2])
sub1.add_patch(PathPatch(path2, facecolor="none", edgecolor='red', linewidth=1))

path3 = Path(numpy.asarray(resampled2.coords)[:, :2])
sub2.add_patch(PathPatch(path3, facecolor="none", edgecolor='red', linewidth=1))

for c in coords1:
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="blue", markersize=10)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="blue", markersize=10)

for c in coords2:
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for c in coords3:
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

sub1.autoscale_view()
plt.show()