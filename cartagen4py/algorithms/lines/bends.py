import numpy as np
import shapely
from shapely import LineString, Point

from cartagen4py.utils.geometry.bends import *
from cartagen4py.utils.geometry.line import *


def accordion(line, width, sigma=None, sample=None):
    """
    Apply the accordion algorithm to a linestring.
    Return a displaced linestring.
    Parameters
    ----------
    line : shapely LineString
        The line to apply the accordion algorithm.
    width : float
        The width in meters of the casing of the symbol.
    sigma : float optional
        Gaussian smoothing arguments -> Gaussian filter strength to derive inflexion points.
        Default to None which will let the gaussian smoothing algorithm derived it from the line.
    sample : int optional
        Gaussian smoothing arguments -> The length in meter between each nodes after resampling the line.
        Default to None which will let the gaussian smoothing algorithm derived it from the line.
    """
    # Detect individual bends inside the smoothed line
    bs = BendSerie(line, sigma, sample)

    # Storage for the future distorded individual lines
    distorded = []

    test = []
    
    for bend in bs.bends:
        # Get the translation vector
        v = __get_vector(bend, width)

        # Retrieve start and end of bend coordinates
        coords = list(bend.bend.coords)
        x1, y1, x2, y2 = coords[0][0], coords[0][1], coords[-1][0], coords[-1][1]

        # Calculate the vector formed by the start and end point of the bend
        vab = np.array([x2 - x1, y2 - y1])

        # Calculate the scalar product between the translation vector
        # and the vector formed by the start and end point of the bend
        scalar = v[0] * vab[0] + v[1] * vab[1]

        # If the scalar product is negative, inverse the translation vector
        if scalar < 0:
            v = np.array([-v[0], -v[1]])

        # Calculate the vector formed by the start and the summit of the bend
        vas = np.array([bend.summit.x - x1, bend.summit.y - y1])

        # Calculate the vector product of vas and vab
        product = vas[0] * vab[1] - vas[1] * vab[0]

        # Check that the closed bend is clockwise
        clockwise = not shapely.is_ccw(LineString(coords + [(x1, y1)]))

        # Inverse the vector by comparating product and clockwise flag
        if clockwise and product > 0:
            v = np.array([-v[0], -v[1]])
        elif clockwise is False and product < 0:
            v = np.array([-v[0], -v[1]])

        # Storage for the translated points
        translated = []

        total = 0
        length = bend.bend.length
        # Loop through coordinates of the bend vertices starting from the end
        for i in range(len(coords) - 1, -1, -1):
            # Get the coordinates of the current point
            p = coords[i]
            # If it's the end point of the bend
            if i == len(coords) - 1:
                # Displace the point by the full vector
                translated.append([p[0] + v[0], p[1] + v[1]])
            else:
                # Here, it's not the end point
                # Increment the current processed length of the line
                total += Point(p).distance(Point(coords[i + 1]))
                # Calculate the portion of the processed length
                factor = (length - total) / length    
                # Calculate a new reduced vector    
                vn = np.array([v[0] * factor, v[1] * factor])
                # Displace the point by the reduced vector and add it at the start of the list
                translated.insert(0, [p[0] + vn[0], p[1] + vn[1]])
        
        # Add the list of translated points to the full list of the whole bend serie
        distorded.append(LineString(translated))

    # Set the displaced line to be the first distorded line
    displaced = distorded[0]
    # Set the next start to be the end of the first line
    lstart = displaced.coords[-1]

    # Loop through remaining distorded lines
    for i in range(1, len(distorded)):
        # Get the line
        dline = distorded[i]
        # Get the start of the line
        dstart = dline.coords[0]
        # Translate the line with dx and dy which are the difference between line start and end of previous line
        tline = shapely.affinity.translate(dline, lstart[0] - dstart[0], lstart[1] - dstart[1])
        # Merge the displaced lines with the current displaced line
        displaced = merge_linestrings(displaced, tline)
        # Set the last start as the current displaced line end point 
        lstart = displaced.coords[-1]

    # Get the middle of the provided bend serie
    bm = get_line_middle_point(line)
    # Get the middle of the distorded line
    dm = get_line_middle_point(displaced)

    # Translate the distorded line by the difference between both middle points
    return shapely.affinity.translate(displaced, bm.x - dm.x, bm.y - dm.y)

def __get_vector(bend, width):
    # Check on both sides if the translation crosses the bend
    def check_translation(base, bend, width, angle, dx, dy):
        # Translate the base point using the vector
        x = base[0] + dx
        y = base[1] + dy
        # Calculate the line between the base point and the translated point
        tline = LineString([Point(base), Point([x, y])])
        # Extend the line by twice the width of the bend on both sides
        eline = extend_line_by_length(tline, bend.width * 2)

        # Create the bend without the base vertex to avoid intersection there
        crossline = LineString([ x for x in list(bend.bend.coords) if x != base ])
        # Calculate the intersection
        intersection = shapely.intersection(eline, crossline)

        # Check that intersection is a point.
        if intersection.geom_type == 'Point':
            # Get the distance between the base point and the intersection
            dist = Point(base).distance(intersection)
            n = np.sqrt(pow(dx, 2) + pow(dy, 2))
            return np.array([(dx / n) * (width - dist), (dy / n) * (width - dist)])
        else:
            # If the interseection is not a point
            # Check that the end point has not been tried already
            if base == list(bend.bend.coords)[-1]:
                # Here, no intersection where found
                return None
            else:
                # Restart the function using the end point as the base
                return check_translation(list(bend.bend.coords)[-1], bend, width, angle, dx, dy)

    # Calculate the angle perpendicular from the orientation of the bend
    o = bend.orientation
    if o < np.pi:
        angle = bend.orientation + np.pi / 2
    else:
        angle = bend.orientation - np.pi / 2
   
    # Create x, y of the vector from the angle and the norm
    dx = bend.width * 2 * np.cos(angle)
    dy = bend.width * 2 * np.sin(angle)

    return check_translation(list(bend.bend.coords)[0], bend, width, angle, dx, dy)