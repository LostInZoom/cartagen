import numpy as np
import shapely
from shapely.ops import split, nearest_points, snap, transform
from shapely.geometry import LineString, Point, Polygon, MultiPoint

from cartagen4py.utils.geometry.angle import *
from cartagen4py.utils.tessellation.hexagonal import *

def visvalingam_whyatt(line, area_tolerance):
    """
    Visvalingam-Whyatt algorithm (1993)
    """
    def contains_another_point(line, pt, index):
        first = line.coords[index-1]
        last = line.coords[index+1]
        triangle_coords = [first, pt, last, first]
        triangle = Polygon(triangle_coords)
        for vertex in line.coords:
            if (vertex == first):
                continue
            if (vertex == pt):
                continue
            if (vertex == last):
                continue
            if (triangle.contains(Point(vertex))):
                return True
        return False
    def compute_area_point(line, pt, index):
        first = line.coords[index-1]
        last = line.coords[index+1]
        triangle = [first, pt, last, first]
        polygon = Polygon(triangle)
        return polygon.area

    final_coords = list(line.coords)
    while (True):
        pt_area_min = None
        area_min = line.envelope.area
        current_line = LineString(tuple(final_coords))
        for i in range(1,len(final_coords)-2):
            pt = final_coords[i]
            if(contains_another_point(current_line, pt, i)):
                continue
            area = compute_area_point(current_line, pt, i)
            if (area < area_min):
                area_min = area
                if(area < area_tolerance):
                    pt_area_min = pt
        if(pt_area_min is None):
            break
        # remove pt from the line
        final_coords.remove(pt_area_min)

    return LineString(tuple(final_coords))

def raposo_simplification(line, initial_scale, final_scale, centroid=True, tobler=False):
    """
    Raposo simplification algorithm (2010): uses an hexagonal tessallation, with a size related to the final scale, 
    and the algorithm only retains one vertex per hexagonal cell. Be careful, it uses the scale as parameter. If the centroid parameter 
    is True, the vertices inside an hexagon cell are replaced by the centroid of the cell; if it is False, they are replaced by the nearest
    vertex to the centroid of the cell.
    """
    width = 0
    if tobler:
        width = final_scale * 5 / 4000 
    else:
        firstFactor = line.length / (len(line.coords)-1)
        secondFactor = final_scale / initial_scale
        width = firstFactor * secondFactor

    # compute hexagon tessellation
    tessellation = HexagonalTessellation(line.envelope, width)

    current_index = 0
    final_coords = []
    # append the first point of the line, but without the z coordinate
    twod_point = transform(lambda *args: args[:2], Point(line.coords[0]))
    final_coords.append(twod_point.coords[0])
    previous_cell = None
    while (current_index < len(line.coords)-1):
        current_cell = None
        # now loop on the vertices from current index
        # builds a point cloud as a multi-point geometry with all line vertices
        # contained in the current cell.
        point_cloud = []
        for i in range(current_index,len(line.coords)):
            # get the cells containing the point
            point = line.coords[i]
            containing_cells = tessellation.get_containing_cells(point)
            if(current_cell is None):
                if(previous_cell in containing_cells):
                    containing_cells.remove(previous_cell)
                current_cell = containing_cells[0]
                point_cloud.append(point)
                continue

            if (current_cell in containing_cells):
                point_cloud.append(point)
                current_index = i+1
            else:
                current_index = i
                previous_cell = current_cell
                break
        multipoint = MultiPoint(point_cloud)
        if (centroid):
            # replace the points by the centroid of the vertices in the cell
            final_coords.append(multipoint.centroid.coords[0])
        else:
            # find the nearest vertex to the centroid
            nearest = nearest_points(multipoint,multipoint.centroid)
            final_coords.append(nearest[0].coords[0])
        if (current_index == len(line.coords) - 1):
            twod_final_point = transform(lambda *args: args[:2], Point(line.coords[current_index]))
            final_coords.append(twod_final_point.coords[0])
    # add the final point if it is not in the line
    twod_final_point = transform(lambda *args: args[:2], Point(line.coords[len(line.coords) - 1]))
    if(twod_final_point.coords[0] not in final_coords):
        final_coords.append(twod_final_point.coords[0])
    
    # print(final_coords)
    return LineString(final_coords)

# Li-Openshaw simplification algorithm (1993). The simplification factor is the size of the regular grid applied on the line
def __li_openshaw_simplification(line, cell_size):
    # TODO: li openshaw line simplification
    return line
    

def gaussian_smoothing(line, sigma=None, sample=None, densify=True, preserve_extremities=False):
    """
    Compute the gaussian smoothing of a set of a LineString.
    Parameters
    ----------
    line : shapely LineString
        The line to smooth.
    sigma : float optional
        Gaussian filter strength.
        Default value to 30, which is a high value.
    sample : float optional
        The length in meter between each nodes after resampling the line.
        If not provided, the sample is derived from the line and is the average distance between
        each consecutive vertex.
    densify : boolean optional
        Whether the resulting line should keep the new node density. Default to True.
    """
    # Extend the given set of points at its first and last points of k points using central inversion.
    def extend(line, interval):
        # Compute the central inversion of a position. origin is the center of symmetry and p is the point to inverse.
        def central_inversion(origin, p):
            x = 2 * origin[0] - p[0]
            y = 2 * origin[1] - p[1]
            return (x,y)
        
        # Get the coordinates of the vertices
        coords = list(line.coords)
        # Get the first and last vertex
        first, last = coords[0], coords[-1]

        # Get the index of the penultimate vertex
        # -2 is to avoid taking the last vertex
        pen = len(coords) - 2

        # Set the start of the line as the central inversion of n first vertices (n = interval)
        result = [central_inversion(first, coords[i]) for i in range(interval, 0, -1)]

        # Add the full line as the middle part of the line
        result.extend(coords)

        # Add the end of the line as the central inversion of n last vertices (n = interval)
        result.extend([central_inversion(last, coords[i]) for i in range(pen, pen - interval, -1)])

        return LineString(result)
   
    if sigma is None:
        sigma = 30

    if sample is None:
        coords = list(line.coords)
        distances = []
        for i in range(0, len(coords) - 1):
            v1, v2 = coords[i], coords[i + 1]
            distances.append(shapely.Point(v1).distance(shapely.Point(v2)))
        avg = (sum(distances) / len(distances))
        sample = avg

    # First resample the line, making sure there is a maximum distance between two consecutive vertices
    resampled = resample_line(line, sample)

    # Calculate the interval (number of vertex to take into consideration when smoothing)
    interval = round(4 * sigma / sample)
    # If the interval is longer than the input line, we change the interval and recalculate the sigma
    if interval >= len(resampled.coords):
        interval = len(resampled.coords) - 1
        sigma = interval * sample / 4
    
    # Compute gaussian coefficients
    c2 = -1.0 / (2.0 * sigma * sigma)
    c1 = 1.0 / (sigma * np.sqrt(2.0 * pi))

    # Compute the gaussian weights and their sum
    weights = []
    total = 0
    for k in range (0, interval + 1):
        weight = c1 * np.exp(c2 * k * k)
        weights.append(weight)
        total += weight
        if k > 0:
            total += weight
    
    # Extend the line at its first and last points with central inversion
    extended = extend(resampled, interval)

    smoothed_coords = []
    for i in range(0, len(resampled.coords)):
        x, y = 0, 0
        for k in range(-interval , interval + 1):
            p1 = extended.coords[i - k + interval]
            x += weights[abs(k)] * p1[0] / total
            y += weights[abs(k)] * p1[1] / total
        smoothed_coords.append((x,y))

    if densify:
        final_coords = smoothed_coords
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
        
    if preserve_extremities:
        final_coords[0] = Point(coords[0])
        final_coords[-1] = Point(coords[-1])
    
    return LineString(final_coords)

def get_bend_side(line):
    """
    Return the side of the interior of the bend regarding the line direction. (left or right)
    Parameters
    ----------
    line : shapely LineString
        The line to get the bend side from.
    """
    # Get the list of nodes
    coords = list(line.coords)

    # The total angle of the bend
    total = 0

    # Get the start point of the bend
    start = shapely.Point(coords[0])

    # Loop through the nodes of the linestring, starting on the seconde node
    for i in range(1, len(coords) - 1):
        # Get the current and the next node coordinates
        c1 = coords[i]
        c2 = coords[i + 1]

        # Add to the total the angle between the starting point and those two nodes
        angle = angle_3_pts(start, shapely.Point(c2), shapely.Point(c1))

        # Set the angle to be between 0 and 2*pi
        if angle % (2 * np.pi) >= 0:
            angle = abs(angle % (2 * np.pi))
        else:
            angle = angle % (2 * np.pi) + 2 * np.pi

        # Set the angle to be between -pi and pi
        if angle > np.pi:
            angle = angle - 2 * np.pi

        total += angle

    # If the total is above 0, the bend is left sided, otherwise it right sided
    if total < 0:
        return 'right'
    else:
        return 'left'

def get_shortest_edge_length(geom: LineString):
    min_length = float('inf')
    segments = get_linestring_segments(geom)
    for segment in segments:
        print(segment)
        length = Point(segment[0]).distance(Point(segment[1]))
        if(length < min_length):
            min_length = length

    return min_length

def get_linestring_segments(line: LineString):
    segments = []
    prev_coord = line.coords[0]
    for coord in line.coords[1:]:
        segments.append((prev_coord,coord))
        prev_coord = coord

    return segments


# for a 3d sequence of coordinates, returns the 2d sequence. Useful to convert 3d points to 2d points
def to_2d(x, y, z):
    return tuple(filter(None, [x, y]))


# for an array of 3d sequences of coordinates, returns an array containing the 2d sequences. Useful to convert 3d geometries to 2d geometries
def array_to_2d(array):
    array_2d = []
    for tuple in array:
       array_2d.append(to_2d(tuple[0],tuple[1],tuple[2]))
    return array_2d

def polygons_3d_to_2d(polygons):
    polygons_2d = []
    for polygon in polygons:
       # outer ring
       coords_2d = array_to_2d(polygon.exterior.coords)
       # inner rings
       interiors = []
       for interior in polygon.interiors:
          interiors.append(array_to_2d(interior.coords))
       polygons_2d.append(Polygon(coords_2d,interiors))
    return polygons_2d

def resample_line(line, step):
    """
    Densifies a line with vertices every 'step' along the line. Keep the start and end vertex.
    Parameters
    ----------
    line : shapely LineString
        The line to densify.
    step : float
        The step (in meters) to resample the geometry
    """

    # Get the length of the line
    length = line.length 
    
    # Storage for final vertices, starting with the start of the line
    xy = [(line.coords[0])]
    
    for distance in np.arange(step, int(length), step):
        # Interpolate a point every step along the old line
        point = line.interpolate(distance)
        # Add the tuple of coordinates
        xy.append((point.x, point.y))
    
    # Add the last point of the line if it doesn't already exist
    if xy[-1] != line.coords[-1]:
        xy.append(line.coords[-1])
    
    # Here, we return a new line with densified points.
    return LineString(xy)

# returns the index of a vertex in a line. Returns -1 if the point given is not a vertex of the line
def get_index_of_vertex(line, vertex, tolerance = 0.01):
    for i in range(0, len(line.coords)):
        point = Point(line.coords[i])
        if point.equals_exact(vertex, tolerance):
            return i
    return -1

# returns the index of the nearest vertex in a line. Returns -1 if the point given is not a vertex of the line
def get_index_of_nearest_vertex(line, vertex):
    mindist = float("inf")
    nearest_index = -1
    for i in range(0, len(line.coords)):
        point = Point(line.coords[i])
        dist = point.distance(vertex)
        if dist < mindist:
            nearest_index = i
            mindist = dist
    return nearest_index

def geometry_flatten(geom):
  if hasattr(geom, 'geoms'):  # Multi<Type> / GeometryCollection
    for g in geom.geoms:
      yield from geometry_flatten(g)
  elif hasattr(geom, 'interiors'):  # Polygon
    yield geom.exterior
    yield from geom.interiors
  else:  # Point / LineString
    yield geom

def geometry_len(geom):
  return sum(len(g.coords) for g in geometry_flatten(geom))

# Gets the nearest point of the geometry to a point. If the point is a vertex, it is not chosen.
def get_nearest_vertex(geometry, point):
    min_dist = float('inf')
    nearest = None
   
    for vertex in geometry.coords:
        vertex_pt = Point(vertex)
        distance = vertex_pt.distance(point)
        if distance == 0:
            continue
        if distance < min_dist:
            min_dist = distance
            nearest = vertex_pt
    
    return nearest

def extend_line_with_point(line, point, position='start'):
    """
    Extend the line with a given point, depending on the position, adds it at the start or the end
    """
    new_line = []
    p = [point.x, point.y]

    if position == 'start':
        new_line.append(p)

    for n in line.coords:
        new_line.append(n)

    if position == 'end':
        new_line.append(p)

    return LineString(new_line)

def extend_line_by_length(line, length, position='both'):
    """
    Extend the line by a given length.
    This function extend the line by adding a vertice at the start and/or the end of the line and extend it by a given length
    following the direction of the first and/or last segment of the line.
    Position is by default set to 'both' which will add the length to both
    the start and the end of the line. Can be set to 'start' or 'end'.
    """

    def __calculate_coords(a, b, length):
        xa, ya, xb, yb = a[0], a[1], b[0], b[1]
        ab = np.sqrt((yb - ya)**2 + (xb -xa)**2)
        xc = xb - (((xb - xa) / ab) * (ab + length))
        yc = yb - (((yb - ya) / ab) * (ab + length))
        return (xc, yc)

    # Get the list of vertex
    vertex = list(line.coords)

    new_line = []

    if position == 'both' or position == 'start':
        new_line.append(__calculate_coords(vertex[0], vertex[1], length))

    new_line.extend(vertex)

    if position == 'both' or position == 'end':
        new_line.append(__calculate_coords(vertex[-1], vertex[-2], length))

    return LineString(new_line)

def get_line_middle_point(line):
    """
    Return the point located along the line at half its length.
    """
    return line.interpolate(line.length / 2)

        
def split_line_at_point(line, point):
    """
    Split a line at a given point along this line.
    Return the two new linestrings.
    Return None if the line and the point doesn't intersect.
    """
    projected = nearest_points(point, line)[1]
    splitted = split(snap(line, projected, 0.0001), projected)

    lines = []
    for s in splitted.geoms:
        lines.append(s)

    return lines[0], lines[1]

def get_segment_center(segment):
    """
    Return the center of the given shapely segment (a linestring with only two vertices) as a shapely Point.
    """
    c = segment.coords

    if len(c) > 2:
        raise Exception("The provided segment has more than two vertices.")

    x1, y1, x2, y2 = c[0][0], c[0][1], c[1][0], c[1][1]
    return shapely.Point([(x1 + x2) / 2, (y1 + y2) / 2])

def project_point_on_line(point, line):
    """
    Project a point on a line. Return the projected point.
    """
    dist = line.project(point)
    # Create the point projected on the line.
    return shapely.Point(list(line.interpolate(dist).coords))

def merge_linestrings(line1, line2):
    """
    Merge two linestring that are connected by their start or end points.
    Return the new linestring. The first line defines the direction of the merged linestring.
    """
    coords1, coords2 = line1.coords, line2.coords

    start1, end1 = coords1[0], coords1[-1]
    start2, end2 = coords2[0], coords2[-1]
    
    # Lines connected by...
    # ...start and start...
    if start1 == start2:
        # Reverse the second line, and concatenate the first (drop the start)
        return shapely.LineString(list(coords2[::-1]) + list(coords1[1:]))
    # ...start and end...
    elif start1 == end2:
        # Concatenate the second and the first line (drop the start)
        return shapely.LineString(list(coords2) + list(coords1[1:]))
    # ...end and start...
    elif end1 == start2:
        # Concatenate the first and the second line (drop the start)
        return shapely.LineString(list(coords1) + list(coords2[1:]))
    # ...end and end
    elif end1 == end2:
        # Reverse the second line (drop the start) and add it to the first line
        return shapely.LineString(list(coords1) + list(coords2[::-1][1:]))
    # Here, linstrings are not connected
    else:
        raise Exception("Provided LineStrings are not connected by their start or end vertex.")

def get_inflexion_points(line, min_dir=120.0):
    """
    This algorithm extract inflexion points from a LineString object and return a list of index of inflexion points.
    Parameters
    ----------
    line : shapely LineString
        The line to extract the inflexion points from.
    min_dir : float
        The minimum direction change (in degrees) between two consecutive inflexion points.
        This parameter allows to remove micro inflexions from the results.
    """

    coords = list(line.coords)

    # Storage for inflexion points
    inflexion = []

    # Storage for parts to check when micro inflexion points are detected
    part = []

    # Stores for the previous angle and the previous direction
    prevangle = None
    prevdir = None

    # Loop through vertices
    for i in range(1, len(coords) - 1):
        # Get previous, current and next point
        p1, p2, p3 = coords[i - 1], coords[i], coords[i + 1]

        # Calculate angles formed by p1 and p2, and p1 and p3
        a1 = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])
        a2 = np.arctan2(p3[1] - p1[1], p3[0] - p1[0])

        # Calculate the direction of p1 p2
        direction = np.degrees(a1) % 360
        # Calculate the angle formed by the three points
        alpha = a2 - a1

        # Store if angle is positive or negative
        angle = 1 if alpha > 0 else -1

        # Check that previous angle has been calculated
        if prevangle is not None:
            # Check if the angle is not the same as the previous -> inflexion point
            if angle != prevangle:
                # Check that the previous direction has been calculated
                if prevdir is not None:
                    # If the absolute difference between the previous direction and the current one is above the direction threshold
                    # It means that this is not a micro inflexion
                    if abs(prevdir - direction) > min_dir:
                        # Adding the middle of the part list, i.e. the middle of the micro inflexions
                        inflexion.append(part[len(part) // 2])
                        # Restart the part list with the current inflexion point
                        part = [i]
                    else:
                        # Add the current point to the part list
                        part.append(i)
                # If not, append the point to the part
                else:
                    part.append(i)

                # Set the previous angle as the direction of the current 
                prevdir = direction
        
        # Set previous angle as current
        prevangle = angle

    # Append the last part's middle point
    inflexion.append(part[len(part) // 2])

    return inflexion
