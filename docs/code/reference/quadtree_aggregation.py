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

o1, g1 = c4.quadtree_aggregation(gdf, depth1, 'value', True)
o2, g2 = c4.quadtree_aggregation(gdf, depth2, 'value', True)
o3, g3 = c4.quadtree_aggregation(gdf, depth3, 'value', True)

fig = plt.figure(1, (8, 12))

sub1 = fig.add_subplot(311)
sub1.set_aspect('equal')
sub1.set_title("depth={0}".format(depth1), pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(312)
sub2.set_aspect('equal')
sub2.set_title("depth={0}".format(depth2), pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

sub3 = fig.add_subplot(313)
sub3.set_aspect('equal')
sub3.set_title("depth={0}".format(depth3), pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
   
for l in g1.geometry(depth1):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub1.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for l in g2.geometry(depth2):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub2.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for l in g3.geometry(depth3):
    path = Path(numpy.asarray(l.coords)[:, :2])
    sub3.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=.5))

for point in o1.to_dict('records'):
    c = point['geometry'].coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['cell_count']*3)

for point in o2.to_dict('records'):
    c = point['geometry'].coords[0]
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['cell_count']*3)

for point in o3.to_dict('records'):
    c = point['geometry'].coords[0]
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['cell_count']*3)

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()

plt.show()