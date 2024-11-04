import shapely
import geopandas as gpd

from cartagen.algorithms.lines.coalescence import coalescence_splitting

from cartagen.utils.geometry.dilation import dilate_line 
from cartagen.utils.debug import plot_debug, geojson_debug

def galbe(line, width):
    """
    Apply GALBE to the provided line.
    """

    split = coalescence_splitting(line, width)
    left, right = dilate_line(line, width)

    splitted = [ s['geometry'] for s in split ]

    plot_debug(left, right, *splitted)

    # return split