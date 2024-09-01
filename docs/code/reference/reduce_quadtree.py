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

depth1 = 1
depth2 = 3
depth3 = 5
width, height = 40, 40
N = 60
coords = numpy.random.randn(N, 2) * height/4 + (width, height)
values = numpy.random.randint(1, 10, N)
points = [ Point(*coord) for coord in coords ]

gdf = gpd.GeoDataFrame([ {'value': values[i], 'geometry': points[i]} for i in range(0, len(points)) ])

o1, g1 = c4.reduce_quadtree(gdf, depth1, 'selection', 'value', True)
o2, g2 = c4.reduce_quadtree(gdf, depth1, 'aggregation', 'value', True)

o3, g3 = c4.reduce_quadtree(gdf, depth2, 'selection', 'value', True)
o4, g4 = c4.reduce_quadtree(gdf, depth2, 'aggregation', 'value', True)

o5, g5 = c4.reduce_quadtree(gdf, depth3, 'selection', 'value', True)
o6, g6 = c4.reduce_quadtree(gdf, depth3, 'aggregation', 'value', True)

fig = plt.figure(1, (12, 16))

sub1 = fig.add_subplot(321)
sub1.set_xticks([])
sub1.set_yticks([])
sub1.set_ylabel("depth={0}".format(depth1))
sub1.set_xlabel("mode='selection'")
sub1.xaxis.set_label_position('top') 

sub2 = fig.add_subplot(322)
sub2.axes.get_yaxis().set_visible(False)
sub2.set_xticks([])
sub2.set_xlabel("mode='aggregation'")
sub2.xaxis.set_label_position('top') 

sub3 = fig.add_subplot(323)
sub3.axes.get_xaxis().set_visible(False)
sub3.set_yticks([])
sub3.set_ylabel("depth={0}".format(depth2))

sub4 = fig.add_subplot(324)
sub4.axes.get_xaxis().set_visible(False)
sub4.axes.get_yaxis().set_visible(False)

sub5 = fig.add_subplot(325)
sub5.axes.get_xaxis().set_visible(False)
sub5.set_yticks([])
sub5.set_ylabel("depth={0}".format(depth3))

sub6 = fig.add_subplot(326)
sub6.axes.get_xaxis().set_visible(False)
sub6.axes.get_yaxis().set_visible(False)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub4.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub5.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub6.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)

for l in g1.geometry(depth1):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub1.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for l in g2.geometry(depth1):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub2.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for l in g3.geometry(depth2):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub3.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for l in g4.geometry(depth2):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub4.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for l in g5.geometry(depth3):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub5.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for l in g6.geometry(depth3):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub6.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for point in o1.to_dict('records'):
    c = point['geometry'].coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=3)

for point in o2.to_dict('records'):
    c = point['geometry'].coords[0]
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['count'])

for point in o3.to_dict('records'):
    c = point['geometry'].coords[0]
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=3)

for point in o4.to_dict('records'):
    c = point['geometry'].coords[0]
    sub4.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['count'])

for point in o5.to_dict('records'):
    c = point['geometry'].coords[0]
    sub5.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=3)

for point in o6.to_dict('records'):
    c = point['geometry'].coords[0]
    sub6.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['count'])

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()
sub4.autoscale_view()
sub5.autoscale_view()
sub6.autoscale_view()

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()