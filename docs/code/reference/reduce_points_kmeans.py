from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
import shapely
from shapely.wkt import loads
import cartagen4py as c4
from shapely import LineString, Polygon, Point
from cartagen4py.utils.debug import plot_debug, geojson_debug, geojson_to_variable

width, height = 40, 40
N = 30
coords = numpy.random.randn(N, 2) * height/4 + (width, height)
values = numpy.random.randint(1, 10, N)

points = [ Point(*coord) for coord in coords ]

simplified1 = c4.reduce_points_kmeans(points, 0.1, True)
simplified2 = c4.reduce_points_kmeans(points, 0.4, True)

fig = plt.figure(1, (12, 24))

sub1 = fig.add_subplot(211)
sub1.set_title('shrink_ratio=0.1', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(212)
sub2.set_title('shrink_ratio=0.4', pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="blue", markersize=5)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="blue", markersize=5)

for s in simplified1:
    c = s.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for s in simplified2:
    c = s.coords[0]
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

sub1.autoscale_view()
sub2.autoscale_view()
plt.show()