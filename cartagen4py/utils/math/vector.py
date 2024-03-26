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
        return cls(segment.point1, segment.point2)

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

    def change_norm(self, norm):
        """
        Return a new vector with a new norm.
        """
        vector = self.copy()
        vector.normalize()
        vector.scalar_multiplication(norm)
        return vector

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
        # if direction % (2 * np.pi) >= 0:
        #     direction = abs(direction % (2 * np.pi))
        # else:
        #     direction = direction % (2 * np.pi) + 2 * np.pi
        unit = Vector2D.from_angle(direction, 1)
        norm = self.scalar_product(unit)
        return unit.change_norm(norm)