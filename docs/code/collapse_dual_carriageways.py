from matplotlib import pyplot as plt
import cartagen4py as c4
import geopandas as gpd

network = gpd.read_file("data/detect_dual_carriageways.geojson")
dual_carriageways = c4.detect_dual_carriageways(network)
collapsed = c4.collapse_dual_carriageways(network, dual_carriageways, sigma=2)

original = network.plot(edgecolor='gray',linewidth=1)
collapsed.plot(ax=original, color='red', linewidth=1)
original.axes.get_xaxis().set_visible(False)
original.axes.get_yaxis().set_visible(False)

plt.show()