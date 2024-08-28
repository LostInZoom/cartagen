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

depth1 = 2
depth2 = 3
width, height = 40, 40
N = 60
coords = numpy.random.randn(N, 2) * height/4 + (width, height)
values = numpy.random.randint(1, 10, N)
points = [ Point(*coord) for coord in coords ]

gdf = gpd.GeoDataFrame([ {'geometry': x} for x in points ])
gdf['value'] = values

output1, qtree1 = c4.reduce_quadtree(gdf, depth1, 'selection', attribute='value')
output2, qtree2 = c4.reduce_quadtree(gdf, depth2, 'selection', attribute='value')

fig = plt.figure(1, (12, 36))

sub1 = fig.add_subplot(311)
sub1.set_title('a) The set of points', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(312)
sub2.set_title('b) Reduction with a depth of {0}'.format(depth1), pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

sub3 = fig.add_subplot(313)
sub3.set_title('b) Reduction with a depth of {0}'.format(depth2), pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

qtree1.draw(sub2, depth1)
qtree2.draw(sub3, depth2)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=5)

for point in output1:
    c = point[0].coords[0]
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

for point in output2:
    c = point[0].coords[0]
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=5)

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()

plt.show()