import math
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, LineString

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

    def get_angle(self):
        """
        Calcule l'angle du vecteur en radians par rapport à l'axe X positif.
        L'angle retourné est compris dans l'intervalle $(-\pi, \pi]$.
        """
        return np.arctan2(self.y, self.x)

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

def interpolate_displacement_vectors(initial, final, interval, crs=3857):
    """
    Calculate initial displacement vectors using interpolation at regular intervals.

    Origin point on initial, extremity on final (same curvilinear abscissa ratio).
    Returns a GeoDataFrame of LineStrings (vectors) and a Series of vector lengths.

    Parameters
    ----------
    initial : LineString
        Geometry before displacement.
    final : LineString
        Geometry after displacement.
    interval : float
        Interval used for the interpolation.
    crs : int, optional
        The CRS that will be used to generate the GeoDataFrame.

    Returns
    -------
    tuple of (GeoDataFrame of LineString, Series of vector length)
    """
    # This function is non-trivial and requires geometry-specific interpolation logic.
    # Placeholder: Assuming simple LineString initiators for demonstration.
    
    if initial.geom_type != 'LineString' or final.geom_type != 'LineString':
        # Must handle complex types (e.g., MultiLineString) in a real implementation
        return gpd.GeoDataFrame(), 0.0

    # Simplified mock for a few vectors
    vectors_data = []
    num_steps = int(initial.length / interval)
    
    for i in range(num_steps + 1):
        ratio = i / num_steps
        
        # Homolog points using same ratio of curvilinear abscissa [cite: 111]
        p_initial = initial.interpolate(initial.length * ratio)
        p_final = final.interpolate(final.length * ratio)
        
        vector = LineString([p_initial, p_final])
        dep_ini = vector.length
        
        vectors_data.append({
            'geometry': vector, 
            'origin': p_initial, 
            'dep_ini': dep_ini,
            'direction': np.array([p_final.x - p_initial.x, p_final.y - p_initial.y]) / dep_ini if dep_ini > 0 else np.array([0.0, 0.0])
        })
    
    vectors_gdf = gpd.GeoDataFrame(vectors_data, crs=crs)
    max_dep = vectors_gdf['dep_ini'].max() if not vectors_gdf.empty else 0.0
    return vectors_gdf, max_dep