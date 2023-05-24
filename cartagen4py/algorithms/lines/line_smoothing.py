# This file contains line smoothing algorithms

from shapely.geometry import LineString, Point
from math import pi, sqrt, exp
from cartagen4py.utils.geometry.line import densify_geometry, to_2d, get_index_of_vertex

# Compute the gaussian smoothing of a set of a LineString. Sigma is the gaussian filter parameter, and threshold is the subsampling parameter
# This code is a port from the GaussianFilter class in the GeOxygene Java library.
# TODO The code needs to fixed
def gaussian_smoothing(line, sigma, threshold):
    # first resample the line, making sure there is a maximum distance between two consecutive vertices
    resampled = densify_geometry(line, threshold)

    interval = round(4*sigma/threshold)
    # if the interval is longer than the input line, we change the interval
    if(interval >= len(resampled.coords)):
        interval = len(resampled.coords) - 1
        sigma = interval * threshold / 4
    
    # compute gaussian coefficients
    c2 = -1.0 / (2.0 * sigma * sigma)
    c1 = 1.0 / (sigma * sqrt(2.0 * pi))
    # compute gassian weights and their sum
    weights = []
    sum = 0
    for k in range (0,interval+1):
        weight = c1 * exp(c2*k*k)
        weights.append(weight)
        sum += weight
        if k>0:
            sum += weight
    
    # extend the line at its first and last points with central inversion
    extended = __extend(resampled,interval)
    smoothed_coords = []
    for i in range(0,len(resampled.coords)):
        x = 0
        y = 0
        for k in range(-interval,interval+1):
            p1 = extended.coords[i-k+interval]
            x += weights[abs(k)]*p1[0] / sum
            y += weights[abs(k)]*p1[1] / sum
        smoothed_coords.append((x,y))
    
    # only return the points matching the input points in the resulting filtered line
    final_coords = []
    for point in line.coords:
        # get the index of point in resampled
        index = get_index_of_vertex(resampled, Point(to_2d(point[0],point[1], point[2])))
        final_coords.append(smoothed_coords[index])
    
    return LineString(final_coords)

# Extend the given set of points at its fist and last points of k points using central inversion.
def __extend(line, interval):
    nb_vert = len(line.coords)
    first = line.coords[0]
    last = line.coords[nb_vert-1]
    new_coords = []
    for i in range(0,nb_vert+2*interval):
        position = i - interval
        if(i < interval):
            p = line.coords[position]
            new_coords.append(__central_inversion(first,p))
        else:
            if(position >= nb_vert):
                beyond = position - nb_vert +1
                p = line.coords[nb_vert-1-beyond]
                new_coords.append(__central_inversion(last,p))
            else:
                p = line.coords[position]
                new_coords.append(p)

    return LineString(new_coords)

# Compute the central inversion of a position. origin is the center of symmetry and p is the point to inverse.
def __central_inversion(origin, p):
    x = 2 * origin[0] - p[0]
    y = 2 * origin[1] - p[1]
    return (x,y)
