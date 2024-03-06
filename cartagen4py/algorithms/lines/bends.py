import shapely

from cartagen4py.algorithms.lines.line_smoothing import *
from cartagen4py.utils.geometry.bends import *

def accordion(line, sigma, sample):
    """
    Apply the accordion algorithm to a linestring.
    """

    # Smooth the line to avoid unnecessary micro inflexion points
    smoothed = gaussian_smoothing(line, sigma, sample, densify=False)

    # Detect individual bends inside the smoothed line
    bs = BendSerie(smoothed)