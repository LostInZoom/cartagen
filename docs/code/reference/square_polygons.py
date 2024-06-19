from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from shapely.wkt import loads
import cartagen4py as c4

buildings = [loads('Polygon ((394999.5 6272975, 395006.70000000001164153 6272962.40000000037252903, 395010.59999999997671694 6272957.5, 394996.59999999997671694 6272944.40000000037252903, 394991 6272949, 394999.20000000001164153 6272956.29999999981373549, 394996.09999999997671694 6272959.70000000018626451, 394998.29999999998835847 6272961.29999999981373549, 394992 6272969.40000000037252903, 394999.5 6272975))'), loads('Polygon ((395007.29999999998835847 6272975.79999999981373549, 395013.20000000001164153 6272981, 395021.20000000001164153 6272969.59999999962747097, 395024.20000000001164153 6272971.90000000037252903, 395031 6272963.79999999981373549, 395020.79999999998835847 6272957.40000000037252903, 395007.29999999998835847 6272975.79999999981373549))')]

fig = plt.figure(1)
sub1 = fig.add_subplot(111)
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

squared = c4.square_polygons(buildings)

for building in buildings:
    poly1 = Path.make_compound_path(Path(numpy.asarray(building.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in building.interiors])
    sub1.add_patch(PathPatch(poly1, facecolor="gray", edgecolor='none'))

for sq in squared:
    poly2 = Path.make_compound_path(Path(numpy.asarray(sq.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in sq.interiors])
    sub1.add_patch(PathPatch(poly2, facecolor="none", edgecolor='red', linewidth=1))

sub1.autoscale_view()
plt.show()