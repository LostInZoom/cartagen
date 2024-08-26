from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from geopandas import GeoDataFrame
import shapely
from shapely.wkt import loads
import cartagen as c4
from shapely import LineString, Polygon, Point
from cartagen.utils.debug import plot_debug, geojson_debug, geojson_to_variable

polygons = [
   loads('Polygon ((-181999.78383096933248453 5373723.77217177487909794, -181933.01503285119542852 5373634.16299416404217482, -181613.30242908780928701 5373857.4101169491186738, -181726.15867044063634239 5373993.88870265055447817, -181819.22560953520587645 5373919.64375638216733932, -181743.60545600406476296 5373822.13367635291069746, -181868.50148552606697194 5373734.31990322005003691, -181911.28392042327322997 5373788.39354557823389769, -181999.78383096933248453 5373723.77217177487909794))')
]

fig = plt.figure(1, (12, 10))

#############################################################

sub1 = fig.add_subplot(111)
# sub1.set_title('buffer=1.0 edge_length=1.0', pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

for polygon in polygons:
    if polygon.geom_type == 'MultiPolygon':
        for polyg in polygon.geoms:
            poly1 = Path.make_compound_path(Path(numpy.asarray(polyg.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in polyg.interiors])
            sub1.add_patch(PathPatch(poly1, facecolor="lightgray", edgecolor='none'))
    else:
        poly2 = Path.make_compound_path(Path(numpy.asarray(polygon.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in polygon.interiors])
        sub1.add_patch(PathPatch(poly2, facecolor="lightgray", edgecolor='none'))

    d = c4.open_polygon(polygon, 30, 1)

    if d.geom_type == 'MultiPolygon':
        for polyg in d.geoms:
            poly3 = Path.make_compound_path(Path(numpy.asarray(polyg.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in polyg.interiors])
            sub1.add_patch(PathPatch(poly3, facecolor="none", edgecolor='red', linewidth=1.5))
    else:
        poly3 = Path.make_compound_path(Path(numpy.asarray(d.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in d.interiors])
        sub1.add_patch(PathPatch(poly3, facecolor="none", edgecolor='red', linewidth=1.5))
    

sub1.autoscale_view()
plt.show()