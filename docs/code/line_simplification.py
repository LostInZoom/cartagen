from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import geopandas as gpd
import numpy

import cartagen4py as c4

lines = gpd.read_file('data/test_lines.geojson')
lines = lines.to_dict('records')

l = []
for line in lines:
    l.append(line['geometry'])

fig = plt.figure(1, (12, 4))

sub1 = fig.add_subplot(131)
sub1.set_title('a) Douglas-Peucker', family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

for line in l:
    path = Path(numpy.asarray(line.coords)[:, :2])
    pathdouglas = Path(numpy.asarray(c4.douglas_peucker(line, 5).coords)[:, :2])

    sub1.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
    sub1.add_patch(PathPatch(pathdouglas, facecolor="none", edgecolor='red', linewidth=1))
    sub1.autoscale_view()

sub2 = fig.add_subplot(132)
sub2.set_title('b) Visvalingam-Whyatt')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

for line in l:
    path = Path(numpy.asarray(line.coords)[:, :2])
    pathdouglas = Path(numpy.asarray(c4.visvalingam_whyatt(line, 200).coords)[:, :2])

    sub2.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
    sub2.add_patch(PathPatch(pathdouglas, facecolor="none", edgecolor='red', linewidth=1))
    sub2.autoscale_view()

sub3 = fig.add_subplot(133)
sub3.set_title('c) Raposo')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

for line in l:
    path = Path(numpy.asarray(line.coords)[:, :2])
    pathdouglas = Path(numpy.asarray(c4.raposo(line, 15000.0, 100000.0).coords)[:, :2])

    sub3.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))
    sub3.add_patch(PathPatch(pathdouglas, facecolor="none", edgecolor='red', linewidth=1))
    sub3.autoscale_view()


plt.show()