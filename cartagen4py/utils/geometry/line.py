import numpy as np
import shapely
from shapely.ops import split, nearest_points, snap
from shapely.geometry import Point, Polygon, LineString

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

# densifies a line with vertices every 'step' along the line
def densify_geometry (line_geometry, step):

    # step: add a vertice every step in whatever unit your coordinate reference system use.
    
    length_m=line_geometry.length # get the length of the line
    
    xy=[] # to store new tuples of coordinates
    
    for distance_along_old_line in np.arange(0,int(length_m),step): 
    
        point = line_geometry.interpolate(distance_along_old_line) # interpolate a point every step along the old line
        xp,yp = point.x, point.y # extract the coordinates
    
        xy.append((xp,yp)) # and store them in xy list
    
    return LineString(xy) # Here, we finally create a new line with densified points.

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
        
def split_line_at_point(line, point):
    """
    Split a line at a given point along this line.
    Return the two new linestrings.
    Return None if the line and the point doesn't intersect.
    """
    if line.distance(point) < 1e-8:
        return None

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