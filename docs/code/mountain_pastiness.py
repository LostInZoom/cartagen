from matplotlib import pyplot as plt
import cartagen as c4
import geopandas as gpd
import shapely

network = gpd.read_file("data/pastiness.geojson")
crs = network.crs
net = network.to_dict('records')

width = { 0: 2, 1: 6, 2: 12 }
tolerance = 60

lines = []
left, right = None, None
for n in net:
    geom = n['geometry']
    lines = c4.coalescence_splitting(geom, tolerance)
    right = { 'geometry': shapely.dilate_line(geom, tolerance) }
    left = { 'geometry': shapely.dilate_line(geom, -tolerance) }

lgpd = gpd.GeoDataFrame([left], crs=crs)
rgpd = gpd.GeoDataFrame([right], crs=crs)
original = lgpd.plot(color='gray', linewidth=1)
rgpd.plot(ax=original, color='gray', linewidth=1)
for line in lines:
    paste = gpd.GeoDataFrame([line], crs=network.crs)
    paste.plot(ax=original, color='red', linewidth=width[line['coalescence']])

original.axes.get_xaxis().set_visible(False)
original.axes.get_yaxis().set_visible(False)

plt.show()