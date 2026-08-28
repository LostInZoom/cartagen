import geopandas as gpd
import matplotlib.pyplot as plt
import cartagen as c4

roads = gpd.GeoDataFrame.from_file('../data/enriched_network.geojson')
selected = c4.mesh_density_selection(roads, 0.1, col_importance='traffic_vo')
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

roads.plot(ax=axes[0], color="steelblue", linewidth=1, aspect="equal")
axes[0].set_title("Initial road network")
axes[0].set_axis_off()

selected.plot(ax=axes[1], color="darkorange", linewidth=1, aspect="equal")
axes[1].set_title("Selected roads")
axes[1].set_axis_off()

plt.tight_layout()
plt.show()