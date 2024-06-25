import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import numpy

def plot_debug(geom, *geoms):
    """
    Plot lines for debugging purposes.
    """    

    fig = plt.figure(1, (12, 12))
    sub1 = fig.add_subplot(111)
    sub1.axes.get_xaxis().set_visible(False)
    sub1.axes.get_yaxis().set_visible(False)
    path = Path(numpy.asarray(geom.coords)[:, :2])
    sub1.add_patch(PathPatch(path, facecolor="none", edgecolor='gray', linewidth=1))

    for g in geoms:
        pathn = Path(numpy.asarray(g.coords)[:, :2])
        sub1.add_patch(PathPatch(pathn, facecolor="none", edgecolor='red', linewidth=1))

    sub1.autoscale_view()
    plt.show()