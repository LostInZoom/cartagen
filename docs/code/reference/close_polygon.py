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
   loads('MultiPolygon (((-181971.10639871368766762 5373438.9800818981602788, -181871.04461579365306534 5373515.03858671616762877, -181830.39398684300249442 5373471.48148276377469301, -181900.21001316665206105 5373416.14983465988188982, -181883.24522059736773372 5373395.05042703542858362, -181921.19564371986780316 5373368.60838820412755013, -181971.10639871368766762 5373438.9800818981602788)),((-181862.91726404300425202 5373412.03444804064929485, -181814.67939064223901369 5373447.85257737524807453, -181781.85609685847884975 5373412.60387528780847788, -181833.1266669845499564 5373374.85333522595465183, -181862.91726404300425202 5373412.03444804064929485)))')
]

fig = plt.figure(1, (12, 10))

#############################################################

sub1 = fig.add_subplot(111)
sub1.set_aspect('equal')
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

    d = c4.close_polygon(polygon, 10)

    if d.geom_type == 'MultiPolygon':
        for polyg in d.geoms:
            poly3 = Path.make_compound_path(Path(numpy.asarray(polyg.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in polyg.interiors])
            sub1.add_patch(PathPatch(poly3, facecolor="none", edgecolor='red', linewidth=1.5))
    else:
        poly3 = Path.make_compound_path(Path(numpy.asarray(d.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in d.interiors])
        sub1.add_patch(PathPatch(poly3, facecolor="none", edgecolor='red', linewidth=1.5))
    

sub1.autoscale_view()
plt.show()