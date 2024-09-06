from matplotlib import pyplot as plt
from matplotlib import colormaps
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy as np
import geopandas as gpd
from shapely.wkt import loads
from shapely import Point
import cartagen as c4

width, height = 40, 40
N = 60
coords = np.random.randn(N, 2) * height/4 + (width, height)
values = np.random.randint(1, 10, N)
points = [ Point(*coord) for coord in coords ]

points_gdf = gpd.GeoDataFrame(geometry=points)
p1, g1 = c4.partition_grid(points_gdf, 10, shape='square')
p2, g2 = c4.partition_grid(points_gdf, 20, shape='diamond')
p3, g3 = c4.partition_grid(points_gdf, 10, shape='hexagonal')

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

cmap = colormaps['tab20']

c1 = cmap(np.linspace(0, 1, len(p1)))
for i, color in enumerate(c1):
    for p in p1[i]:
        coords = points[p].coords[0]
        sub1.plot(coords[0], coords[1], linestyle="", marker='o', color=color, markersize=3)
    for c in g1:
        poly = Path.make_compound_path(Path(np.asarray(c.exterior.coords)[:, :2]),
            *[Path(np.asarray(ring.coords)[:, :2]) for ring in c.interiors])
        sub1.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray', linewidth=0.5))

c2 = cmap(np.linspace(0, 1, len(p2)))
for i, color in enumerate(c2):
    for p in p2[i]:
        coords = points[p].coords[0]
        sub2.plot(coords[0], coords[1], linestyle="", marker='o', color=color, markersize=3)
    for c in g2:
        poly = Path.make_compound_path(Path(np.asarray(c.exterior.coords)[:, :2]),
            *[Path(np.asarray(ring.coords)[:, :2]) for ring in c.interiors])
        sub2.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray', linewidth=0.5))

c3 = cmap(np.linspace(0, 1, len(p3)))
for i, color in enumerate(c3):
    for p in p3[i]:
        coords = points[p].coords[0]
        sub3.plot(coords[0], coords[1], linestyle="", marker='o', color=color, markersize=3)
    for c in g3:
        poly = Path.make_compound_path(Path(np.asarray(c.exterior.coords)[:, :2]),
            *[Path(np.asarray(ring.coords)[:, :2]) for ring in c.interiors])
        sub3.add_patch(PathPatch(poly, facecolor="none", edgecolor='gray', linewidth=0.5))

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()
plt.show()