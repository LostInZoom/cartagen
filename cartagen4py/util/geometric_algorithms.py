from shapely.geometry import Point, Polygon, LineString
import math

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

# resample a LineString with at least one vertex every max_distance
def resampling_line(line, max_distance):
    resampled = []
    prev_point = line.coords[0]
    resampled.append(prev_point)
    for i in range (1,len(line.coords)-1):
        next_point = line.coords[i]
        dist = LineString([prev_point,next_point]).length
        fseg = dist / max_distance
        # TODO

    return LineString(resampled)

# for a 3d sequence of coordinates, returns the 2d sequence. Useful to convert 3d points to 2d points
def to_2d(x, y, z):
    return tuple(filter(None, [x, y]))