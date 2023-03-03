# This file contains line smoothing algorithms

from shapely import LineString, Point

# Compute the gaussian smoothing of a set of a LineString. Sigma is the gaussian filter parameter, and threshold is the subsampling parameter
# This code is a port from the GaussianFilter class in the GeOxygene Java library.
def gaussian_smoothing(line, sigma, threshold):
    # TODO
    return LineString()