import geopandas as gpd
import shapely

from cartagen4py.algorithms.lines.line_smoothing import *

def detect_bends(line, sigma, threshold):
    """
    Detect bends from a given LineString.
    """

    smoothed = gaussian_smoothing(line, sigma, threshold)

    return smoothed