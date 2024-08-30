import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import geopandas as gpd
import numpy
import random

def plot_debug(geom, *geoms):
    """
    Plot geometries for debugging purposes.
    """
    def add(sub, geom, color1, color2, linewidth=1):
        if geom is not None:
            if geom.geom_type[:5] == 'Multi':
                for g in geom.geoms:
                    add_to_plot(sub, g, color1, color2, linewidth=1)
            else:
                add_to_plot(sub, geom, color1, color2, linewidth=1)

    def add_to_plot(sub, geom, color1, color2, linewidth=1):
        if geom.geom_type == 'Polygon':
            poly1 = Path.make_compound_path(Path(numpy.asarray(geom.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in geom.interiors])
            sub.add_patch(PathPatch(poly1, facecolor=color1, edgecolor=color2, linewidth=linewidth))
        elif geom.geom_type == 'LineString':
            path = Path(numpy.asarray(geom.coords)[:, :2])
            sub.add_patch(PathPatch(path, facecolor="none", edgecolor=color1, linewidth=linewidth))
        elif geom.geom_type == 'Point':
            sub.plot(geom.coords[0][0], geom.coords[0][1], linestyle="", marker='o', color=color1)

    fig = plt.figure(1)
    sub = fig.add_subplot(111)
    sub.axes.get_xaxis().set_visible(False)
    sub.axes.get_yaxis().set_visible(False)

    color1 = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
    color2 = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
    if isinstance(geom, list):
        for g in geom:
            if isinstance(g, dict):
                add(sub, g['geometry'], color1, color2)
            else:
                add(sub, g, color1, color2)
    elif isinstance(geom, gpd.GeoDataFrame):
        for g in geom.geometry:
            add(sub, g, color1, color2)
    else:
        add(sub, geom, color1, color2)
    
    for geom in geoms:
        color1 = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
        color2 = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
        if isinstance(geom, list):
            for g in geom:
                if isinstance(g, dict):
                    add(sub, g['geometry'], color1, color2)
                else:
                    add(sub, g, color1, color2)
        elif isinstance(geom, gpd.GeoDataFrame):
            for g in geom.geometry:
                add(sub, g, color1, color2)
        else:
            add(sub, geom, color1, color2)

    sub.autoscale_view(tight=True)
    fig.tight_layout()
    plt.show()

def geojson_debug(geom, *geoms):
    """
    Plot geometries for debugging purposes.
    """
    name = 0
    skip = False

    gdf = []
    if isinstance(geom, list):
        for g in geom:
            if isinstance(g, dict):
                gdf.append(g)
            else:
                gdf.append({'geometry' : g})
    elif isinstance(geom, gpd.GeoDataFrame):
        geom.to_file('{0}.geojson'.format(name), driver='GeoJSON')
        name += 1
        skip = True
    else:
        gdf.append({'geometry' : geom})

    if not skip:
        gdf = gpd.GeoDataFrame(gdf, crs='3857')
        gdf.to_file('{0}.geojson'.format(name), driver='GeoJSON')
        name += 1

    for geom in geoms:
        skip = False
        gdf = []
        if isinstance(geom, list):
            for g in geom:
                if isinstance(g, dict):
                    gdf.append(g)
                else:
                    gdf.append({'geometry' : g})
        elif isinstance(geom, gpd.GeoDataFrame):
            geom.to_file('{0}.geojson'.format(name), driver='GeoJSON')
            name += 1
            skip = True
        else:
            gdf.append({'geometry' : geom})

        if not skip:
            gdf = gpd.GeoDataFrame(gdf, crs='3857')
            gdf.to_file('{0}.geojson'.format(name), driver='GeoJSON')
            name += 1

def geojson_to_variable(geojson, write=False, attributes=None):
    """
    Generate a variable string to directly load the shapely geometry inside the Python file.
    """
    v = gpd.read_file(geojson)
    v = v.to_dict('records')

    f = None
    if write:
        f = open('myfile.txt', 'w')

    s = "variable = ["
    for o in v:
        if attributes is not None:
            entry = "{"
            for k, v in o.items():
                if k in attributes:
                    entry += "'{0}': {1},".format(k, v)
            entry += "'geometry': loads('{0}')".format(o['geometry'].wkt)
            entry += "}, \n"
            if write:
                f.write(entry)
        else:
            entry = "loads('{0}'), \n".format(o['geometry'].wkt)
            if write:
                f.write(entry)
            s += entry
    
    s += "]"
    if not write:
        print(s)