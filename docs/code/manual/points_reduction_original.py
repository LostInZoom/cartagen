from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from shapely.wkt import loads
import cartagen as c4

land = gpd.GeoDataFrame.from_file('../data/land_australia.geojson')
cities = gpd.GeoDataFrame.from_file('../data/cities_australia.geojson')

fig = plt.figure(1, (10, 9))

sub1 = fig.add_subplot(111)
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

for l in land.geometry:
    poly = Path.make_compound_path(Path(numpy.asarray(l.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in l.interiors])
    sub1.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))

for c in cities.geometry:
    sub1.plot(c.coords[0][0], c.coords[0][1], linestyle="", marker='o', color="black", markersize=2)

sub1.autoscale_view()
plt.show()