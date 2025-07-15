from matplotlib import pyplot as plt
from matplotlib import colormaps
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from geopandas import GeoDataFrame
import shapely
from shapely.wkt import loads
import cartagen as c4
from shapely import LineString, Polygon, Point
from cartagen.utils.debug import plot_debug, geojson_debug, geojson_to_variable

polygon = Polygon([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
line = polygon.boundary

lsmoothed = c4.gaussian_smoothing(line, 5, 1)
psmoothed = c4.gaussian_smoothing(polygon, 5, 1)

fig = plt.figure(1, (10, 10))

sub1 = fig.add_subplot(221)
sub1.set_aspect('equal')
sub1.set_title('a) Original line', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

path = Path(numpy.asarray(line.coords)[:, :2])
sub1.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
for point in line.coords:
    sub1.plot(point[0], point[1], linestyle="", marker='o', color="gray")
sub1.autoscale_view()

#############################################################################

sub2 = fig.add_subplot(222)
sub2.set_aspect('equal')
sub2.set_title('b) Original polygon', pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

poly = Path.make_compound_path(Path(numpy.asarray(polygon.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in polygon.interiors])
sub2.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))
for point in line.coords:
    sub2.plot(point[0], point[1], linestyle="", marker='o', color="gray")
sub2.autoscale_view()

#############################################################################

sub3 = fig.add_subplot(223)
sub3.set_aspect('equal')
sub3.set_title('c) Smoothed line', pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

path = Path(numpy.asarray(lsmoothed.coords)[:, :2])
sub3.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
for point in lsmoothed.coords:
    sub3.plot(point[0], point[1], linestyle="", marker='o', color="gray")
sub3.autoscale_view()

#############################################################################

sub4 = fig.add_subplot(224)
sub4.set_aspect('equal')
sub4.set_title('d) Smoothed polygon', pad=10, family='sans-serif')
sub4.axes.get_xaxis().set_visible(False)
sub4.axes.get_yaxis().set_visible(False)

poly = Path.make_compound_path(Path(numpy.asarray(psmoothed.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in psmoothed.interiors])
sub4.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))
for point in psmoothed.exterior.coords[:-1]:
    sub4.plot(point[0], point[1], linestyle="", marker='o', color="gray")
sub4.autoscale_view()

plt.show()