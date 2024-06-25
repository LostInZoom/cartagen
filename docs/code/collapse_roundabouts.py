from matplotlib import pyplot as plt
import cartagen4py as c4
import geopandas as gpd

network = gpd.read_file("data/roundabouts.geojson")
roundabouts = c4.detect_roundabouts(network)
branching = c4.detect_branching_crossroads(network)

collapsedr = c4.collapse_roundabouts(network, roundabouts, branching)
collapsedb = c4.collapse_branching_crossroads(network, branching, roundabouts)

original = network.plot(edgecolor='gray',linewidth=1)
collapsedb.plot(ax=original, edgecolor='red', linewidth=1)
original.axes.get_xaxis().set_visible(False)
original.axes.get_yaxis().set_visible(False)

plt.show()