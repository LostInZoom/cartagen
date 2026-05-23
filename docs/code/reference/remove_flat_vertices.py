from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
import shapely
from shapely.wkt import loads
import cartagen as c4

building = loads('Polygon ((483162.60420556488679722 6044749.90797633398324251, 483167.29305910837138072 6044741.50208449736237526, 483158.54631337692262605 6044736.39747273270040751, 483164.02214775263564661 6044727.26683164667338133, 483160.10401825612643734 6044724.26727752014994621, 483158.94143492548028007 6044723.54596658423542976, 483156.44691042724298313 6044721.9793161153793335, 483152.05093409406254068 6044730.33950771857053041, 483148.83314704877557233 6044727.98296463023871183, 483143.6222503071767278 6044734.96569846104830503, 483162.60420556488679722 6044749.90797633398324251))')
removed = c4.remove_flat_vertices(building, 15)

fig = plt.figure(1, (12, 10))

#############################################################

sub1 = fig.add_subplot(121)
sub1.set_aspect('equal')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(122)
sub2.set_aspect('equal')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

poly = Path.make_compound_path(Path(numpy.asarray(building.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in building.interiors])
sub1.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray'))

poly = Path.make_compound_path(Path(numpy.asarray(removed.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in removed.interiors])
sub2.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray'))

for p in building.boundary.coords:
    sub1.plot(p[0], p[1], linestyle="", marker='o', color="gray")

for p in removed.boundary.coords:
    sub2.plot(p[0], p[1], linestyle="", marker='o', color="gray")

sub1.autoscale_view()
plt.show()