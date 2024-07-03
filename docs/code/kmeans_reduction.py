from matplotlib import pyplot as plt
import cartagen4py as c4
import geopandas as gp
from shapely.geometry import Point

points = [Point(1,1), Point(1,2), Point(0,1), Point(2,1), Point(2,2), Point(5,5), Point(8,10), Point(10,10), Point(10,8), 
              Point(16,10), Point(16,9), Point(14,11)]
simplified = c4.reduce_points_kmeans(points, 0.25, True)
p1 = gp.GeoSeries(points)
p2 = gp.GeoSeries(simplified)
base = p1.plot()
p2.plot(ax=base, color='red')
plt.show()