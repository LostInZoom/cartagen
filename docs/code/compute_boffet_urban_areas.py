from matplotlib import pyplot as plt
import cartagen4py as c4
import geopandas as gp

zipfile = 'data/buildings.zip'
df = gp.read_file(zipfile)
polygons = df.geometry

urbanareas = c4.compute_boffet_urban_areas(polygons,25.0, 10.0)
p1 = gp.GeoSeries(polygons)
p2 = gp.GeoSeries(urbanareas)
base = p2.plot()
p1.plot(ax=base, color='red')
plt.show()