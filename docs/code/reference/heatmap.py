from matplotlib import pyplot as plt
from matplotlib import colormaps
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
import shapely
from shapely.wkt import loads
import cartagen as c4
from shapely import LineString, Polygon, Point
from cartagen.utils.debug import plot_debug, geojson_debug, geojson_to_variable

cell_size = 1
radius1, radius2, radius3 = 2, 5, 10
width, height = 40, 40
N = 60
coords = numpy.random.randn(N, 2) * height/4 + (width, height)
values = numpy.random.randint(1, 10, N)
points = [ Point(*coord) for coord in coords ]

gdf = gpd.GeoDataFrame([ {'value': values[i], 'geometry': points[i]} for i in range(0, len(points)) ])

o1 = c4.heatmap(gdf, cell_size, radius1, 'value')
o2 = c4.heatmap(gdf, cell_size, radius2, 'value')
o3 = c4.heatmap(gdf, cell_size, radius3, 'value')

r1 = o1.to_dict('records')
r2 = o2.to_dict('records')
r3 = o3.to_dict('records')

fig = plt.figure(1, (12, 12))

sub1 = fig.add_subplot(221)
sub1.set_title("a) The set of points", pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(222)
sub2.set_title("b) radius={0}".format(radius1), pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

sub3 = fig.add_subplot(223)
sub3.set_title("c) radius={0}".format(radius2), pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

sub4 = fig.add_subplot(224)
sub4.set_title("d) radius={0}".format(radius3), pad=10, family='sans-serif')
sub4.axes.get_xaxis().set_visible(False)
sub4.axes.get_yaxis().set_visible(False)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=2)

cmap = colormaps['coolwarm']

d1 = [ x['density'] for x in r1 ]
norm1 = ((d1 - numpy.min(d1)) / (numpy.max(d1) - numpy.min(d1))).tolist()
for i1, p in enumerate(r1):
    geom = p['geometry']
    poly = Path.make_compound_path(Path(numpy.asarray(geom.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in geom.interiors])
    sub2.add_patch(PathPatch(poly, facecolor=cmap(norm1[i1]), edgecolor='none'))

d2 = [ x['density'] for x in r2 ]
norm2 = ((d2 - numpy.min(d2)) / (numpy.max(d2) - numpy.min(d2))).tolist()
for i2, p in enumerate(r2):
    geom = p['geometry']
    poly = Path.make_compound_path(Path(numpy.asarray(geom.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in geom.interiors])
    sub3.add_patch(PathPatch(poly, facecolor=cmap(norm2[i2]), edgecolor='none'))

d3 = [ x['density'] for x in r3 ]
norm3 = ((d3 - numpy.min(d3)) / (numpy.max(d3) - numpy.min(d3))).tolist()
for i3, p in enumerate(r3):
    geom = p['geometry']
    poly = Path.make_compound_path(Path(numpy.asarray(geom.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in geom.interiors])
    sub4.add_patch(PathPatch(poly, facecolor=cmap(norm3[i3]), edgecolor='none'))

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()
sub4.autoscale_view()

plt.show()