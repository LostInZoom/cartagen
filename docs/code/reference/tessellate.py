from matplotlib import pyplot as plt
from matplotlib import colormaps
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy as np
import geopandas as gpd
from shapely.wkt import loads
import cartagen as c4

extent = [ 0, 0, 10, 10 ]

g1 = c4.tessellate(extent, 1, shape='square')
g2 = c4.tessellate(extent, 1, shape='diamond')
g3 = c4.tessellate(extent, 1, shape='hexagonal')

fig = plt.figure(1, (12, 4))

sub1 = fig.add_subplot(131)
sub1.set_title("a) shape='square'", pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(132)
sub2.set_title("b) shape='diamond'", pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

sub3 = fig.add_subplot(133)
sub3.set_title("c) shape='hexagonal'", pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

for c in g1:
    poly = Path.make_compound_path(Path(np.asarray(c.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in c.interiors])
    sub1.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray', linewidth=0.5))

for c in g2:
    poly = Path.make_compound_path(Path(np.asarray(c.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in c.interiors])
    sub2.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray', linewidth=0.5))

for c in g3:
    poly = Path.make_compound_path(Path(np.asarray(c.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in c.interiors])
    sub3.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray', linewidth=0.5))

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()
plt.show()