from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from shapely.wkt import loads
import cartagen4py as c4

line = loads('LineString (-187545.95508891576901078 5374353.07019249629229307, -187519.07677357134525664 5374362.02963094506412745, -187501.1578966750530526 5374375.46878861729055643, -187484.3589495847991202 5374393.38766551297158003, -187478.75930055469507352 5374413.54640202131122351, -187476.51944094264763407 5374435.94499814230948687, -187477.63937074868590571 5374460.58345387410372496, -187483.23901977876084857 5374488.58169902488589287, -187482.11908997275168076 5374506.50057592149823904, -187469.79986210656352341 5374528.89917204156517982, -187454.12084482228965499 5374546.81804893724620342, -187435.08203811998828314 5374555.7774873860180378, -187409.32365258157369681 5374552.41769796796143055, -187394.76456510333810002 5374533.37889126501977444, -187393.64463529732893221 5374512.10022495128214359, -187404.84393335750792176 5374490.82155863661319017, -187419.40302083574351855 5374470.6628221282735467, -187429.48238908991334029 5374443.78450678382068872, -187433.96210831397911534 5374413.54640202131122351, -187437.32189773203572258 5374385.54815687146037817, -187435.08203811998828314 5374362.02963094506412745, -187429.48238908991334029 5374342.99082424212247133, -187413.80337180563947186 5374318.35236851032823324, -187394.76456510333810002 5374298.19363200198858976)')

fig = plt.figure(1, (12, 4))

offset = 10
sub1 = fig.add_subplot(131)
sub1.set_title('a) Offset = {0}'.format(offset), pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

path1 = Path(numpy.asarray(line.coords)[:, :2])
path2 = Path(numpy.asarray(c4.min_break(line, offset).coords)[:, :2])
sub1.add_patch(PathPatch(path1, facecolor="none", edgecolor='gray', linewidth=1))
sub1.add_patch(PathPatch(path2, facecolor="none", edgecolor='red', linewidth=1))
sub1.autoscale_view()

offset = 30
sub2 = fig.add_subplot(132)
sub2.set_title('a) Offset = {0}'.format(offset), pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

path1 = Path(numpy.asarray(line.coords)[:, :2])
path2 = Path(numpy.asarray(c4.min_break(line, offset).coords)[:, :2])
sub2.add_patch(PathPatch(path1, facecolor="none", edgecolor='gray', linewidth=1))
sub2.add_patch(PathPatch(path2, facecolor="none", edgecolor='red', linewidth=1))
sub2.autoscale_view()

offset = 50
sub3 = fig.add_subplot(133)
sub3.set_title('a) Offset = {0}'.format(offset), pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

path1 = Path(numpy.asarray(line.coords)[:, :2])
path2 = Path(numpy.asarray(c4.min_break(line, offset).coords)[:, :2])
sub3.add_patch(PathPatch(path1, facecolor="none", edgecolor='gray', linewidth=1))
sub3.add_patch(PathPatch(path2, facecolor="none", edgecolor='red', linewidth=1))
sub3.autoscale_view()

plt.show()