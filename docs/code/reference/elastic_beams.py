import geopandas as gpd
from shapely.wkt import loads
import cartagen as c4
from matplotlib import pyplot as plt

roads = [
        loads('LINESTRING (435657.1 6227901.0, 435663.8 6227869.8, 435669.7 6227843.7, ' \
        '435671.2 6227827.5, 435670.6 6227806.6, 435657.9 6227774.5, 435654.6 6227750.9, ' \
        '435648.7 6227690.1, 435646.1 6227602.4, 435659.6 6227552.6, 435666.4 6227529.8, ' \
        '435665.9 6227506.1, 435643.6 6227462.3, 435623.3 6227432.8, 435581.2 6227408.4)'),
        loads('LINESTRING (435677.3 6227890.9, 435688.3 6227816.7, 435690.8 6227802.3, ' \
        '435693.4 6227791.4, 435688.3 6227772.8, 435679.0 6227730.6, 435671.4 6227701.9, ' \
        '435664.7 6227644.6, 435662.0 6227601.6, 435671.4 6227566.0, 435681.6 6227546.7, ' \
        '435700.0 6227526.5, 435716.1 6227503.7, 435719.5 6227492.7, 435722.9 6227472.5, ' \
        '435727.1 6227450.5, 435727.0 6227432.8, 435727.1 6227417.6, 435723.7 6227404.1, ' \
        '435722.9 6227404.0)')
    ]
lines = gpd.GeoDataFrame([{'geometry': g, 'importance': 1} for g in roads])
p2 = c4.beams_displacement(lines, 50, verbose = True, iterations=30)

base = lines.plot(facecolor='none', edgecolor='blue')
base.axes.get_xaxis().set_visible(False)
base.axes.get_yaxis().set_visible(False)
p2.plot(ax=base, facecolor='none', edgecolor='red')
base.autoscale_view()
plt.show()