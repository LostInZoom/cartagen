from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
import shapely
from shapely.wkt import loads
import cartagen as c4

points = [
    loads('POINT (-187538.24561624724 5374649.667424928)'), 
    loads('POINT (-187514.09592713122 5374668.987176221)'), 
    loads('POINT (-187493.16619656398 5374646.447466379)'), 
]

fig = plt.figure(1, (12, 4))

#############################################################

sub1 = fig.add_subplot(121)
sub1.set_title("a) offset=10 cap_style='round'", pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

offset = c4.offset_line(shapely.LineString(points), 10)

for p in points:
    sub1.plot(p.coords[0][0], p.coords[0][1], linestyle="", marker='o', color="black")

for o in offset:
    for p in o['projected']:
        sub1.plot(p[0], p[1], linestyle="", marker='o', color="red")

sub1.autoscale_view()

# #############################################################

sub2 = fig.add_subplot(122)
sub2.set_title("b) offset=10 cap_style='flat'", pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

offset = c4.offset_line(shapely.LineString(points), 10, cap_style='flat')

for p in points:
    sub2.plot(p.coords[0][0], p.coords[0][1], linestyle="", marker='o', color="black")

for o in offset:
    for p in o['projected']:
        sub2.plot(p[0], p[1], linestyle="", marker='o', color="red")

sub2.autoscale_view()

plt.show()