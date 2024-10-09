import shapely
import geopandas as gpd

from cartagen.algorithms.lines.coalescence import coalescence_splitting

def galbe(line, width):
    """
    Apply GALBE to the provided line.
    """

    split = coalescence_splitting(line, width)
    print(split)
    return split