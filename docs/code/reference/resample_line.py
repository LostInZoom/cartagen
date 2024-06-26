from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from shapely.wkt import loads
import cartagen4py as c4

line = loads('LineString (-187714.88334235321963206 5374960.62342192325741053, -187618.28458588907960802 5375091.72173426765948534, -187494.08618472088710405 5374864.02466545905917883, -187354.01798784791026264 5375088.96176979690790176, -187229.12959556211717427 5374852.98480757791548967, -187120.80099009876721539 5374972.35327092278748751)')
fig = plt.figure(1, (12, 8))

sub1 = fig.add_subplot(111)
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

resampled = c4.resample_line(line, 50)
coords1 = list(line.coords)
coords2 = list(resampled.coords)

path1 = Path(numpy.asarray(line.coords)[:, :2])
sub1.add_patch(PathPatch(path1, facecolor="none", edgecolor='gray', linewidth=1))

for c in coords1:
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="black")

for c in coords2:
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="red")

sub1.autoscale_view()
plt.show()