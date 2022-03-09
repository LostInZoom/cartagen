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

def main():
    line = LineString([Point(0.0, 1.0), (2.0, 3.0), Point(5.0, 5.0)])
    length = get_shortest_edge_length(line)
    print(length)

if __name__ == "__main__":
    main()