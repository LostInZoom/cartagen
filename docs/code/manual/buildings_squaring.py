from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
import shapely
from shapely.wkt import loads
import cartagen as c4

buildings = [
    loads('Polygon ((286245.81888399075251073 6249040.75842231791466475, 286250.11808830907102674 6249031.19752174336463213, 286260.70718818274326622 6249036.57995702791959047, 286255.94798312475904822 6249047.05115610081702471, 286249.44365562056191266 6249043.6687324745580554, 286247.1403015754185617 6249048.82868163380771875, 286233.22071577416500077 6249042.21097574569284916, 286236.29076149506727234 6249035.53389179147779942, 286245.81888399075251073 6249040.75842231791466475))'),
    loads('Polygon ((286254.76777102699270472 6249023.44433117751032114, 286261.51791547652101144 6249009.48524082824587822, 286269.38724610098870471 6249012.87513171602040529, 286273.68480439885752276 6249003.6186762535944581, 286284.72294759127544239 6249010.06847814749926329, 286281.19538637605728582 6249017.19932754337787628, 286284.97258918144507334 6249019.95831919647753239, 286279.91083130962215364 6249030.27553514018654823, 286274.46097477461444214 6249028.2680484177544713, 286272.77178870484931394 6249032.06208168156445026, 286254.76777102699270472 6249023.44433117751032114))'),
    loads('Polygon ((286228.65751166746485978 6249023.67769361659884453, 286236.62014532077591866 6249009.87721555121243, 286247.05744614545255899 6249015.25872234255075455, 286240.61083843186497688 6249029.2196182981133461, 286228.65751166746485978 6249023.67769361659884453))'),
]

fig = plt.figure(1, (12, 5))

#############################################################

sub1 = fig.add_subplot(121)
sub1.set_aspect('equal')
sub1.set_title('a) Least square squaring', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(122)
sub2.set_aspect('equal')
sub2.set_title('b) Naive squaring', pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

for building in buildings:
    poly = Path.make_compound_path(Path(numpy.asarray(building.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in building.interiors])
    sub1.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))
    sub2.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))
    
for building in buildings:
    generalized = c4.square_polygon_ls(building)
    poly = Path.make_compound_path(Path(numpy.asarray(generalized.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in generalized.interiors])
    sub1.add_patch(PathPatch(poly, facecolor="none", edgecolor='red', linewidth=1.5))

    generalized = c4.square_polygon_naive(building, orient='swo', remove_flat=False)
    poly = Path.make_compound_path(Path(numpy.asarray(generalized.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in generalized.interiors])
    sub2.add_patch(PathPatch(poly, facecolor="none", edgecolor='red', linewidth=1.5))

sub1.autoscale_view()
sub2.autoscale_view()
plt.show()