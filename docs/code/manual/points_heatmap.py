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

land = gpd.GeoDataFrame.from_file('../data/land_australia.geojson')
cities = gpd.GeoDataFrame.from_file('../data/cities_australia.geojson')

heatmap = c4.heatmap(cities, 50000, 500000, 'pop_max', clip=land)
r = heatmap.to_dict('records')

fig = plt.figure(1, (10, 9))

sub1 = fig.add_subplot(111)
sub1.set_aspect('equal')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

cmap = colormaps['coolwarm']
d1 = [ x['density'] for x in r ]
norm1 = ((d1 - numpy.min(d1)) / (numpy.max(d1) - numpy.min(d1))).tolist()
for i1, p in enumerate(r):
    geom = p['geometry']
    if geom.geom_type == 'MultiPolygon':
        for g in geom.geoms:
            poly = Path.make_compound_path(Path(numpy.asarray(g.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in g.interiors])
            sub1.add_patch(PathPatch(poly, facecolor=cmap(norm1[i1]), edgecolor='none'))
    else:
        poly = Path.make_compound_path(Path(numpy.asarray(geom.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in geom.interiors])
        sub1.add_patch(PathPatch(poly, facecolor=cmap(norm1[i1]), edgecolor='none'))

sub1.autoscale_view()
plt.show()