from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from shapely.wkt import loads
from shapely import Polygon
import cartagen as c4

p1 = Polygon([(0,0), (2,0), (2,2), (0,2)])
p2 = Polygon([(3,0.5), (5,0.8), (4.5,2.5), (2.8,2)]) # Un peu incliné

agg = c4.building_amalgamation(p1, p2)

fig = plt.figure(1, (12, 12))

sub1 = fig.add_subplot(111)
sub1.set_aspect('equal')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

poly = Path.make_compound_path(Path(numpy.asarray(p1.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in p1.interiors])
sub1.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))

poly = Path.make_compound_path(Path(numpy.asarray(p2.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in p2.interiors])
sub1.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))

poly = Path.make_compound_path(Path(numpy.asarray(agg.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in agg.interiors])
sub1.add_patch(PathPatch(poly, facecolor="none", edgecolor='red', linewidth=1.5))

sub1.autoscale_view()
plt.show()