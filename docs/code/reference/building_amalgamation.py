from matplotlib import pyplot as plt

from shapely.geometry import Polygon
import cartagen as c4

p1 = Polygon([(0,0), (2,0), (2,2), (0,2)])
p2 = Polygon([(3,0.5), (5,0.8), (4.5,2.5), (2.8,2)]) # Un peu incliné

agg = c4.building_amalgamation(p1, p2)

fig, ax = plt.subplots()
for p, c in [(p1, 'blue'), (p2, 'green'), (agg, 'red')]:
    x, y = p.exterior.xy
    ax.fill(x, y, alpha=0.3, fc=c, ec='black')
plt.show()