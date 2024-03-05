import geopandas as gpd
import shapely

from cartagen4py.algorithms.lines.line_smoothing import *
from cartagen4py.utils.geometry.angle import *
from cartagen4py.utils.geometry.line import *

# from test_functions import *

def detect_bends(line, sigma, threshold):
    """
    Detect bends from a given LineString.
    """

    # Create the gaussian smoothed line
    smoothed = gaussian_smoothing(line, sigma, threshold)

    # Get inflexion points indexes
    indexes = inflexion_points(smoothed)
    
    # make_gdf([{'geometry': shapely.Point(list(smoothed.coords)[x])} for x in indexes], 'inflexion')

    return smoothed