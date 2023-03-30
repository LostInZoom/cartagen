# This is an implementation of the least squares method proposed by Lars E. Harrie (1999)

import numpy, shapely
from cartagen4py.util.partitioning.network import network_partition

class ConstraintMethod:
    """
    Initialize constraint method object, with default constraints
    """
    def __init__(self, default_distance=20, max_iteration=1000, network_partitioning=True, verbose=False):
        self.MAX_ITER = max_iteration
        self.DEFAULT_DISTANCE = default_distance
        self.NETWORK_PARTITIONING = network_partitioning
        self.VERBOSE = verbose

        self.__ALLOWED_WEIGHTS = {
            'Point': ['movement', 'spatial'],
            'LineString': ['movement', 'stiffness', 'curvature', 'spatial'],
            'Polygon': ['movement', 'stiffness', 'curvature', 'spatial']
        }

        self._DEFAULT_WEIGHTS = {

        }
        
        self.__OBJECTS_NUMBER = 0

        self.__POINTS = []
        self.__IMMOVABLE_POINTS = []

        self.__LINES = []
        self.__STIFF_LINES = []
        self.__IMMOVABLE_LINES = []

        self.__POLYGONS = []
        self.__STIFF_POLYGONS = []
        self.__IMMOVABLE_POLYGONS = []

        self.__OBJECTS = {}

        

        self.__current_coordinates_points = {}
        self.__current_points_index = []

    def add(self, *objects, distance=None, movable=None, flexible=None, **weights):
        """
        Add one or multiple geographic object for preparation before the generalization.
        Parameters
        ----------
        object : Geopandas GeoSerie
            A GeoSerie of geographic objects, can be points, lines or polygons (if multigeometry are provided, they will be exploded.)
        movable : boolean, optional
            Indicates whether the geometry of the object is allowed to move. By default, points and polygons are movable, lines are not.
        flexible : boolean, optional
            Indicates whether the geometry of the object can be altered. Only applies to lines and polygons. By default, polygons and lines are not flexible.
        weights : Specify a weight value for a specific constraint:
            movement : int
                If an object is able to move, specify the weight of the movement constraint.
            stiffness : int
                If a line or polygon is able to move but is not flexible, specify the weight of the stiffness constraint.
            curvature : int
                If a line or polygon is flexible, specify the weight of the curvature constraint.
            conflict : int
                If a point, line or polygon spatial conflicts with other objects (i.e. is within the distance parameter), specify the weight of the spatial conflict constraint.
        """

        for o in objects:
            # Explode geometries to avoid multigeometries
            exploded = o.explode(ignore_index=True)
            # Retrieve the geometry type from the first object
            geometry = exploded.geom_type[0]
            # Check if specified weights apply
            self.__check_weights(geometry, weights)

            # Preparing the objects
            self.__prepare_objects(exploded)

            # Checks if provided arguments are correct and sets parameters
            if geometry == 'Point':
                if flexible == True:
                    raise Exception('A point object cannot be flexible, its internal geometry being only one point.')
            if movable is None:
                if geometry == 'LineString':
                    movable = False
                else:
                    movable = True
            if flexible is None:
                flexible = False

            # Increment the objects number
            self.__OBJECTS_NUMBER += 1
        

    def generalize(self):
        """
        Launch the constraint generalization.
        Can take different types of entry
        """
        if self.__OBJECTS_NUMBER < 1:
            raise Exception('No objects provided, cannot generalize.')

        for i in range(self.MAX_ITER):
            pass

    def __check_weights(self, geomtype, weights):
        """
        Check whether a weight is allowed for a specific geometry type
        """
        allowed = self.__ALLOWED_WEIGHTS[geomtype]
        for w in weights:
            if w not in allowed:
                raise Exception('{0} weight does not exists or cannot be applied to {1}.\n'.format(w, geomtype) +
                    'Available weights for {0}: {1}'.format(geomtype, ', '.join(allowed)))

    def __prepare_objects(self, objects):
        """
        Prepare the geometries of the provided object
        """
        # Creating a list of a list of the object index and its geometry
        shapes = []
        for i, o in objects.iterrows():
            shapes.append(o.geometry)

        self.__build_points_from_shape(shapes)
        self.__build_shapes_from_points(shapes)
        self.__build_observation_matrice(list(self.__current_coordinates_points))

    def __build_points_from_shape(self, shapes):
        """
        Populate the current dictionnary where :
        - The key is the tuple of each point coordinates of the objects (x, y)
        - The value is a list of each shape index it belongs to.
        """
        for i, s in enumerate(shapes):
            coordinates = self.__get_coordinates(s)
            for p in coordinates:
                if p in self.__current_coordinates_points:
                    if i not in self.__current_coordinates_points[p]:
                        self.__current_coordinates_points[p].append(i)
                else:
                    self.__current_coordinates_points[p] = [i]
    
    def __build_shapes_from_points(self, shapes):
        """
        Populate the current list where each geometry is made of a list of index
        of the points it is made of.
        """
        for s in shapes:
            points = []
            coordinates = self.__get_coordinates(s)
            for p in coordinates:
                points.append(self.__get_rank(self.__current_coordinates_points, p))
            self.__current_points_index.append(points)

    def __get_rank(self, points, point):
        """
        Return the index 
        """
        for i, p in enumerate(points):
            if p == point:
                return i
        return -1

    def __get_coordinates(self, shape):
        """
        Returns a list of coordinates from a shapely geometry.
        """
        if shape.geom_type == 'Polygon':
            return shape.exterior.coords
        else:
            return shape.coords

    def __build_observation_matrice(self, points):
        """
        Build the observation matrice with all unique points.
        """
        nb = len(points)
        self.Y = numpy.zeros(2 * nb)

        for i, p in enumerate(points):
            self.Y[2*i] = p[0]
            self.Y[2*i+1] = p[1]

    def __build_weight_matrice(self, points, weights):
        """
        Build the weight matrice.
        """