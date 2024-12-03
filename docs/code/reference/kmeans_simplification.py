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

o1 = c4.kmeans_simplification(gdf, ratio1)
o2 = c4.kmeans_simplification(gdf, ratio2)
o3 = c4.kmeans_simplification(gdf, ratio3)

fig = plt.figure(1, (8, 12))

sub1 = fig.add_subplot(311)
sub1.set_aspect('equal')
sub1.set_title("ratio={0}".format(ratio1), pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(312)
sub2.set_aspect('equal')
sub2.set_title("ratio={0}".format(ratio2), pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

sub3 = fig.add_subplot(313)
sub3.set_aspect('equal')
sub3.set_title("ratio={0}".format(ratio3), pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)

for point in o1.to_dict('records'):
    if point['selected_kmeans']:
        c = point['geometry'].coords[0]
        sub1.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for point in o2.to_dict('records'):
    if point['selected_kmeans']:
        c = point['geometry'].coords[0]
        sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for point in o3.to_dict('records'):
    if point['selected_kmeans']:
        c = point['geometry'].coords[0]
        sub3.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()

plt.show()