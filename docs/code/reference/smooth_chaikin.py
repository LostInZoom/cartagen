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

line = loads('LINESTRING (290682.99807003745809197 5645165.06180777028203011, 290884.54200811864575371 5645085.63542265351861715, 290805.12668338674120605 5645198.35101594962179661, 290734.80616105266381055 5645320.61769961565732956, 290764.8178957705385983 5645355.40702610928565264, 290960.98510244640056044 5645213.5829611960798502, 290974.72192761034239084 5645135.65565379150211811, 291066.68295895465416834 5645116.24101543426513672, 291153.86838414391968399 5645149.98778768535703421, 291315.10353460890473798 5645031.44846810027956963, 291311.96432496851775795 5645279.56204771902412176, 291359.73151846788823605 5645338.67530056741088629, 291281.36259694944601506 5645388.55487615056335926, 291223.85494800563901663 5645337.88606899976730347, 291160.13567147561116144 5645481.55892462935298681, 291323.76419099263148382 5645488.86732831597328186)')
a = [1, 2, 3, 4]

fig = plt.figure(1, (10, 10))

sub1 = fig.add_subplot(221)
sub1.set_title(f'iteration={a[0]}', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

path = Path(numpy.asarray(line.coords)[:, :2])
sub1.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
path = Path(numpy.asarray(c4.smooth_chaikin(line, a[0]).coords)[:, :2])
sub1.add_patch(PathPatch(path, facecolor="none", edgecolor='red', linewidth=1))

sub1.autoscale_view()

#############################################################################

sub2 = fig.add_subplot(222)
sub2.set_title(f'iteration={a[1]}', pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

path = Path(numpy.asarray(line.coords)[:, :2])
sub2.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
path = Path(numpy.asarray(c4.smooth_chaikin(line, a[1]).coords)[:, :2])
sub2.add_patch(PathPatch(path, facecolor="none", edgecolor='red', linewidth=1))

sub2.autoscale_view()

#############################################################################

sub3 = fig.add_subplot(223)
sub3.set_title(f'iteration={a[2]}', pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

path = Path(numpy.asarray(line.coords)[:, :2])
sub3.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
path = Path(numpy.asarray(c4.smooth_chaikin(line, a[2]).coords)[:, :2])
sub3.add_patch(PathPatch(path, facecolor="none", edgecolor='red', linewidth=1))

sub3.autoscale_view()

#############################################################################

sub4 = fig.add_subplot(224)
sub4.set_title(f'iteration={a[3]}', pad=10, family='sans-serif')
sub4.axes.get_xaxis().set_visible(False)
sub4.axes.get_yaxis().set_visible(False)

path = Path(numpy.asarray(line.coords)[:, :2])
sub4.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
path = Path(numpy.asarray(c4.smooth_chaikin(line, a[3]).coords)[:, :2])
sub4.add_patch(PathPatch(path, facecolor="none", edgecolor='red', linewidth=1))

sub4.autoscale_view()

plt.show()