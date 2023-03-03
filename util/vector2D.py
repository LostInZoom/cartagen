import math
from shapely.geometry import Point

class Vector2D:
    x=0.0
    y=0.0

    # constructor with shapely points
    def __init__(self, point1, point2):
        self.x = point2.x - point1.x
        self.y = point2.y - point1.y

    def __init__(self, segment):
        self.x = segment.point2[0] - segment.point1[0]
        self.y = segment.point2[1] - segment.point1[1]
    
    def translate(self, point):
        return Point(point.x + self.x, point.y + self.y)