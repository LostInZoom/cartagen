from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from shapely.wkt import loads
import cartagen as c4

a = (1, 1)
b = (1, 2)
c = (2, 1)

interpolation = c4.circle_interpolation(a, b, c, quad_segs=8)

fig = plt.figure(1, (8, 8))

sub1 = fig.add_subplot(111)
sub1.set_aspect('equal')
sub1.set_title('quad_segs=8', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub1.plot(a[0], a[1], linestyle="", marker='o', color="black")
sub1.annotate('a', (1.01, 1.01))
sub1.plot(b[0], b[1], linestyle="", marker='o', color="black")
sub1.annotate('b', b, (1.01, 1.97))
sub1.plot(c[0], c[1], linestyle="", marker='o', color="black")
sub1.annotate('c', (1.97, 1.01))

for i in range(1, len(interpolation) - 1):
    p = interpolation[i]
    sub1.plot(p[0], p[1], linestyle="", marker='o', color="red")

sub1.autoscale_view()
plt.show()