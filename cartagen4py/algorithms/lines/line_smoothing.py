# This file contains line smoothing algorithms

from shapely.geometry import LineString, Point
from math import pi, sqrt, exp
from cartagen4py.utils.geometry.line import densify_geometry, to_2d, get_index_of_nearest_vertex

def gaussian_smoothing(line, sigma, sample, densify=True):
    """
    Compute the gaussian smoothing of a set of a LineString.
    Parameters
    ----------
    line : shapely LineString
        The line to smooth.
    sigma : float
        Gaussian filter strength.
    sample : int
        The length in meter between each nodes after resampling the line.
    densify : boolean optional
        Whether the resulting line should keep the new node density. Default to True.
    """
    # First resample the line, making sure there is a maximum distance between two consecutive vertices
    resampled = densify_geometry(line, sample)

    # Calculate the interval (number of vertex to take into consideration when smoothing)
    interval = round(4 * sigma / sample)
    # If the interval is longer than the input line, we change the interval and recalculate the sigma
    if interval >= len(resampled.coords):
        interval = len(resampled.coords) - 1
        sigma = interval * sample / 4
    
    # Compute gaussian coefficients
    c2 = -1.0 / (2.0 * sigma * sigma)
    c1 = 1.0 / (sigma * sqrt(2.0 * pi))

    # Compute the gaussian weights and their sum
    weights = []
    total = 0
    for k in range (0, interval + 1):
        weight = c1 * exp(c2 * k * k)
        weights.append(weight)
        total += weight
        if k > 0:
            total += weight
    
    # Extend the line at its first and last points with central inversion
    extended = __extend(resampled, interval)

    smoothed_coords = []
    for i in range(0, len(resampled.coords)):
        x, y = 0, 0
        for k in range(-interval , interval + 1):
            p1 = extended.coords[i - k + interval]
            x += weights[abs(k)] * p1[0] / total
            y += weights[abs(k)] * p1[1] / total
        smoothed_coords.append((x,y))

    if densify:
        return LineString(smoothed_coords)
    else:
        # Only return the points matching the input points in the resulting filtered line
        final_coords = []
        # Stores for index of already treated vertices
        done = []
        # Loop through initial vertices
        for point in list(line.coords):
            # Set the distance to infinite
            distance = float("inf")
            nearest = None
            # Loop through smoothed coordinates
            for i in range(len(smoothed_coords)):
                # Check that the index has not been already added
                if i not in done:
                    # Calculate distance from the point
                    d = Point(smoothed_coords[i]).distance(Point(point))
                    if d < distance:
                        # Update distance and nearest index if below existing
                        distance, nearest = d, i

            # If a nearest point has been found, add it to the new line
            if nearest is not None:
                final_coords.append(smoothed_coords[nearest])
                # Add the index as treated already
                done.append(nearest)
            else:
                final_coords.append(point)

        return LineString(final_coords)

# Extend the given set of points at its first and last points of k points using central inversion.
def __extend(line, interval):
    # Get the coordinates of the vertices
    coords = list(line.coords)
    # Get the first and last vertex
    first, last = coords[0], coords[-1]

    # Get the index of the penultimate vertex
    # -2 is to avoid taking the last vertex
    pen = len(coords) - 2

    # Set the start of the line as the central inversion of n first vertices (n = interval)
    result = [__central_inversion(first, coords[i]) for i in range(interval, 0, -1)]

    # Add the full line as the middle part of the line
    result.extend(coords)

    # Add the end of the line as the central inversion of n last vertices (n = interval)
    result.extend([__central_inversion(last, coords[i]) for i in range(pen, pen - interval, -1)])

    return LineString(result)

# Compute the central inversion of a position. origin is the center of symmetry and p is the point to inverse.
def __central_inversion(origin, p):
    x = 2 * origin[0] - p[0]
    y = 2 * origin[1] - p[1]
    return (x,y)
