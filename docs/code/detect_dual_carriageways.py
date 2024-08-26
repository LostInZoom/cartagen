from matplotlib import pyplot as plt
import cartagen as c4
import geopandas as gpd

network = gpd.read_file("data/detect_dual_carriageways.geojson")
dual_carriageways = c4.detect_dual_carriageways(network)

original = network.plot(edgecolor='gray',linewidth=1)
dual_carriageways.plot(ax=original, color='red', linewidth=1)
original.axes.get_xaxis().set_visible(False)
original.axes.get_yaxis().set_visible(False)

plt.show()