from shapely.geometry import Point, Polygon, LineString
import numpy as np

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
def get_index_of_vertex(line, vertex):
    for i in range(0, len(line.coords)):
        point = Point(line.coords[i])
        if point.equals(vertex):
            return i
    return -1

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