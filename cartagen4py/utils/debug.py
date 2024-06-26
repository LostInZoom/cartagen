import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import geopandas as gpd
import numpy
import random

def plot_debug(geom, *geoms):
    """
    Plot lines for debugging purposes.
    """
    def add(sub, geom, color, linewidth=1):
        if geom.geom_type == 'Polygon':
            poly1 = Path.make_compound_path(Path(numpy.asarray(geom.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in geom.interiors])
            sub.add_patch(PathPatch(poly1, facecolor=color, edgecolor=color, linewidth=linewidth))
        elif geom.geom_type == 'LineString':
            path = Path(numpy.asarray(geom.coords)[:, :2])
            sub.add_patch(PathPatch(path, facecolor="none", edgecolor=color, linewidth=linewidth))
        elif geom.geom_type == 'Point':
            sub.plot(geom.coords[0][0], geom.coords[0][1], linestyle="", marker='o', color=color)

    fig = plt.figure(1)
    sub = fig.add_subplot(111)
    sub.axes.get_xaxis().set_visible(False)
    sub.axes.get_yaxis().set_visible(False)

    color = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
    if isinstance(geom, list):
        for g in geom:
            if isinstance(g, dict):
                add(sub, g['geometry'], color)
            else:
                add(sub, g, color)
    elif isinstance(geom, gpd.GeoDataFrame):
        for g in geom.geometry:
            add(sub, g, color)
    else:
        add(sub, geom, color)
    
    for geom in geoms:
        color = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
        if isinstance(geom, list):
            for g in geom:
                if isinstance(g, dict):
                    add(sub, g['geometry'], color)
                else:
                    add(sub, g, color)
        elif isinstance(geom, gpd.GeoDataFrame):
            for g in geom.geometry:
                add(sub, g, color)
        else:
            add(sub, geom, color)

    sub.autoscale_view()
    plt.show()

def geojson_to_variable(geojson):
    """
    Generate a variable string to directly load the shapely geometry inside the Python file.
    """
    v = gpd.read_file(geojson)
    v = v.to_dict('records')

    s = "variable = ["
    for o in v:
        geom = "loads('{0}'), \n".format(o['geometry'].wkt)
        s += geom
    
    s += "]"
    print(s)