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

ratio1, ratio2, ratio3 = 0.7, 0.4, 0.1
width, height = 40, 40
N = 60
coords = numpy.random.randn(N, 2) * height/4 + (width, height)
values = numpy.random.randint(1, 10, N)
points = [ Point(*coord) for coord in coords ]

gdf = gpd.GeoDataFrame([ {'value': values[i], 'geometry': points[i]} for i in range(0, len(points)) ])

o1 = c4.reduce_kmeans(gdf, ratio1, 'selection', 'value')
o2 = c4.reduce_kmeans(gdf, ratio1, 'aggregation', 'value')

o3 = c4.reduce_kmeans(gdf, ratio2, 'selection', 'value')
o4 = c4.reduce_kmeans(gdf, ratio2, 'aggregation', 'value')

o5 = c4.reduce_kmeans(gdf, ratio3, 'selection', 'value')
o6 = c4.reduce_kmeans(gdf, ratio3, 'aggregation', 'value')

plt.rcParams.update({
    'font.size': 12,
    'axes.labelpad': 10
})

fig = plt.figure(1, (12, 16))

sub1 = fig.add_subplot(321)
sub1.set_xticks([])
sub1.set_yticks([])
sub1.set_ylabel("ratio={0}".format(ratio1))
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
sub3.set_ylabel("ratio={0}".format(ratio2))

sub4 = fig.add_subplot(324)
sub4.axes.get_xaxis().set_visible(False)
sub4.axes.get_yaxis().set_visible(False)

sub5 = fig.add_subplot(325)
sub5.axes.get_xaxis().set_visible(False)
sub5.set_yticks([])
sub5.set_ylabel("ratio={0}".format(ratio3))

sub6 = fig.add_subplot(326)
sub6.axes.get_xaxis().set_visible(False)
sub6.axes.get_yaxis().set_visible(False)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)
    sub4.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)
    sub5.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)
    sub6.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)

for point in o1.to_dict('records'):
    c = point['geometry'].coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for point in o2.to_dict('records'):
    c = point['geometry'].coords[0]
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['count']*3)

for point in o3.to_dict('records'):
    c = point['geometry'].coords[0]
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for point in o4.to_dict('records'):
    c = point['geometry'].coords[0]
    sub4.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['count']*3)

for point in o5.to_dict('records'):
    c = point['geometry'].coords[0]
    sub5.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for point in o6.to_dict('records'):
    c = point['geometry'].coords[0]
    sub6.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['count']*3)

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()
sub4.autoscale_view()
sub5.autoscale_view()
sub6.autoscale_view()

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()