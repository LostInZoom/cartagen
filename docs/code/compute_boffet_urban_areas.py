from matplotlib import pyplot as plt
import cartagen4py as c4
import geopandas as gp

zipfile = 'data/buildings.zip'
df = gp.read_file(zipfile)
polygons = df.geometry

urbanareas = c4.boffet_areas(polygons,25.0, 10.0)
p1 = gp.GeoSeries(polygons)
p2 = gp.GeoSeries(urbanareas)
base = p2.plot(facecolor='none', edgecolor='red')
p1.plot(ax=base, facecolor='gray', edgecolor='none')
base.set_title('Urban areas (Boffet, 2000)', pad=10, family='sans-serif')
base.axes.get_xaxis().set_visible(False)
base.axes.get_yaxis().set_visible(False)
plt.show()