from matplotlib import pyplot as plt
import cartagen4py as c4
import geopandas as gp
from shapely.geometry import Point
import numpy as np

DPI = 72
width, height = 40, 40

N = 60
coords = np.random.randn(N, 2) * height/4 + (width, height)
values = np.random.randint(1, 10, N)
points = [Point(*coord) for coord in coords]
p1 = gp.GeoSeries(points)
gdf = gp.GeoDataFrame(geometry=gp.GeoSeries(p1))
gdf['value'] = values

output, qtree = c4.quadtree_point_set_reduction(gdf, 2, 'selection', attribute='value')

xmin, ymin, xmax, ymax = qtree.envelope.bounds
length = max(xmax - xmin, ymax - ymin)
ax = plt.subplot()
ax.set_xlim(xmin, xmin + length)
ax.set_ylim(ymin, ymin + length)
qtree.draw(ax)
ax.scatter([p[0].x for p in output], [p[0].y for p in output], s=8, c='red')
ax.scatter([p.x for p in points], [p.y for p in points], s=4)
ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
plt.show()