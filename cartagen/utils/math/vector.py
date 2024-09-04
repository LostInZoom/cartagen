import math
import numpy as np
from shapely.geometry import Point

class Vector2D:
    """
    Create a 2D vector.
    """
    def __init__(self, point1, point2):
        self.x = point2.x - point1.x
        self.y = point2.y - point1.y

    @classmethod
    def from_point(cls, point):
        """
        Create a vector between 0 and the given point.
        """
        return cls(Point(0, 0), point)

    @classmethod
    def from_points(cls, point1, point2):
        """
        Create a vector between the two given points.
        """
        return cls(point1, point2)

    @classmethod
    def from_segment(cls, segment):
        """
        Create a vector from a segment.
        """
        point1 = Point(segment.point1[0], segment.point1[1])
        point2 = Point(segment.point2[0], segment.point2[1])
        return cls(point1, point2)

    @classmethod
    def from_angle(cls, angle, norm):
        """
        Create a vector from an angle and a norm.
        """
        return cls(Point(0, 0), Point(norm * np.cos(angle), norm * np.sin(angle)))

    def copy(self):
        """
        Copy the vector.
        """
        return Vector2D.from_point(Point(self.x, self.y))

    
    def translate(self, point):
        """
        Translate the given point with the vector.
        """
        return Point(point.x + self.x, point.y + self.y)

    def get_norm(self):
        """
        Return the norm of the vector.
        """
        return np.sqrt(pow(self.x, 2) + pow(self.y, 2))

    def normalize(self):
        """
        Normalize the vector.
        """
        n = self.get_norm()
        self.x = self.x / n
        self.y = self.y / n

    def opposite(self):
        """
        Return the opposite of the vector.
        """
        return Vector2D.from_point(Point(-self.x, -self.y))

    def angle(self, vector):
        """
        Calculate the angle with the provided vector.
        """
        value = self.scalar_product(vector) / (self.get_norm() * vector.get_norm())
        if value < -1 or value > 1:
            value = round(value)

        if not (self.get_norm() * vector.get_norm()) == 0:
            angle = np.arccos(value)
            if not np.isnan(angle):
                return angle
        return 0

    def change_norm(self, norm):
        """
        Return a new vector with a new norm.
        """
        vector = self.copy()
        vector.normalize()
        vector.scalar_multiplication(norm)
        return vector

    def const_product(self, constant):
        """
        Return a new vector as the product of the current one with a constant.
        """
        return Vector2D.from_point(Point(self.x * constant, self.y * constant))

    def scalar_product(self, v):
        """
        Return the scalar product with a vector.
        """
        return self.x * v.x + self.y * v.y

    def scalar_multiplication(self, value):
        """
        Multiply the vector by a given value.
        """
        self.x = self.x * value
        self.y = self.y * value

    def add(self, vector):
        """
        Add the vector and return a new one.
        """
        return Vector2D.from_point(Point(self.x + vector.x, self.y + vector.y))

    def product(self, vector):
        """
        Calculate the vector product between this vector and a given one.
        """
        return self.x * vector.y - self.y * vector.x

    def project(self, direction):
        """
        Project the vector in a given direction.
        """
        if direction >= (np.pi * 2):
            direction = direction - np.pi * 2
        unit = Vector2D.from_angle(direction, 1)
        norm = self.scalar_product(unit)
        return unit.change_norm(norm)

    def rotate(self, angle):
        """
        Returns the rotated vector by the angle value.
        """
        cos = np.cos(angle)
        sin = np.sin(angle)
        x = self.x * cos - self.y * sin
        y = self.x * sin - self.y * cos
        return Vector2D.from_point(Point(x, y))