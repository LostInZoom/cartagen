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

width, height = 40, 40
N = 60
coords = numpy.random.randn(N, 2) * height/4 + (width, height)
values = numpy.random.randint(1, 10, N)
points = [ Point(*coord) for coord in coords ]

gdf = gpd.GeoDataFrame([ {'value': values[i], 'geometry': points[i]} for i in range(0, len(points)) ])

o1, g1 = c4.labelgrid_selection(gdf, 5, 5, 'value', 'square', grid=True)
o2, g2 = c4.labelgrid_selection(gdf, 5, 5, 'value', 'diamond', grid=True)
o3, g3 = c4.labelgrid_selection(gdf, 5, 5, 'value', 'hexagonal', grid=True)

fig = plt.figure(1, (8, 26))

sub1 = fig.add_subplot(311)
sub1.set_aspect('equal')
sub1.set_title("grid='square'", pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(312)
sub2.set_aspect('equal')
sub2.set_title("grid='diamond'", pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

sub3 = fig.add_subplot(313)
sub3.set_aspect('equal')
sub3.set_title("grid='hexagonal'", pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

for point in points:
    c = point.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="gray", markersize=1)

for cell in g1.to_dict('records'):
    c = cell['geometry']
    poly = Path.make_compound_path(Path(numpy.asarray(c.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in c.interiors])
    sub1.add_patch(PathPatch(poly, facecolor="none", edgecolor='lightgray'))

for cell in g2.to_dict('records'):
    c = cell['geometry']
    poly = Path.make_compound_path(Path(numpy.asarray(c.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in c.interiors])
    sub2.add_patch(PathPatch(poly, facecolor="none", edgecolor='lightgray'))

for cell in g3.to_dict('records'):
    c = cell['geometry']
    poly = Path.make_compound_path(Path(numpy.asarray(c.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in c.interiors])
    sub3.add_patch(PathPatch(poly, facecolor="none", edgecolor='lightgray'))

for point in o1.to_dict('records'):
    if point['selected_labelgrid']:
        c = point['geometry'].coords[0]
        sub1.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['cell_count']*3)

for point in o2.to_dict('records'):
    if point['selected_labelgrid']:
        c = point['geometry'].coords[0]
        sub2.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['cell_count']*3)

for point in o3.to_dict('records'):
    if point['selected_labelgrid']:
        c = point['geometry'].coords[0]
        sub3.plot(c[0], c[1], linestyle="", marker='o', color="red", markersize=point['cell_count']*3)

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()

plt.show()