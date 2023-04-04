# This is an implementation of the least squares method proposed by Lars E. Harrie (1999)

import numpy, shapely
from cartagen4py.util.partitioning.network import network_partition

class ConstraintMethod:
    """
    Initialize constraint method object, with default constraints
    Parameters
        ----------
        default_distance : int optional
            This is the default distance for the buffers created to detect spatial conflicts.
            This value will be used with all objects added to the generalization if no distance is provided.
        
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

        self.__DEFAULT_WEIGHTS = {
            # By default, points can move and are influenced by spatial conflicts
            'Point': {
                'movement': 50,
                'spatial': 50
            },
            # By default, lines are immovable (-1 correspond to positive infinite)
            'LineString': {
                'movement': -1,
                'spatial': 50
            },
            # By default, polygons can move and are stiff, i.e. their internal geometry must be invariant
            'Polygon': {
                'movement': 50,
                'stiffness': 200,
                'spatial': 50
            }
        }
        
        self.__OBJECTS = []
        self.__WEIGHTS = []
        self.__DISTANCES = []
        self.__SHAPE_COUNT = 0
        
        self.__current_coordinates_points = {}
        self.__current_points_index = []

    def generalize(self):
        """
        Launch the constraint generalization on the added objects.
        """

        # Checks if objects are present
        if len(self.__OBJECTS) < 1:
            raise Exception('No objects provided, cannot generalize.')

        shapes = []
        for obj in self.__OBJECTS:
            for i, o in obj.iterrows():
                shapes.append(o.geometry)

        print(self.__WEIGHTS)
        print(self.__DISTANCES)
        
        self.__build_points_from_shape(shapes)
        self.__build_shapes_from_points(shapes)
        unique = list(self.__current_coordinates_points)
        self.__build_observation_matrice(unique)

        for i in range(self.MAX_ITER):
            pass

    def add(self, *objects, distance=None, **weights):
        """
        Add one or multiple geographic object for preparation before the generalization. If multiple objects are provided, they must have the same geometry type.
        Parameters
        ----------
        object : **Geopandas**, *GeoSerie*.
            One or multiple GeoSerie of geographic objects, can be points, lines or polygons (if multigeometry are provided, they will be exploded).
            If multiple objects are provided, they must be the same geometry type because the same constraints will be applied.
        distance : **int**, *optional*.
            The distance around the object to detect spatial conflicts. If no distance is provided, the default one from the class is used.
        weights : **int**, *optional*.
            Specify a weight value for a specific constraint. A value of -1 is infinite. Possible weights:
            movement : **int**.
                If an object is able to move, specify the weight of the movement constraint.
            stiffness : **int**.
                If a line or polygon is able to move but is not flexible, specify the weight of the stiffness constraint. A point cannot have a stiffness constraint.
            curvature : **int**.
                If a line or polygon is flexible, specify the weight of the curvature constraint. A point cannot have a curvature constraint.
            conflict : **int**.
                If a point, line or polygon spatially conflicts with other objects (i.e. is within the distance parameter), specify the weight of the spatial conflict constraint.
        """

        geometry = None
        for i, obj in enumerate(objects):
            # Explode geometries to avoid multigeometries
            o = obj.explode(ignore_index=True)
            
            # If multiple objects, checking if they have the same geometry type
            if i > 0 and o.geom_type[0] != geometry:
                raise Exception("Provided objects have different geometry type.")
            
            # Retrieve the geometry type from the first object
            geometry = o.geom_type[0]
            
            # Check if weights are correctly specified
            w = self.__reconstruct_weights(geometry, **weights)
            self.__WEIGHTS.append(w)

            # Create a new field and populate it with a unique id based on the total number of shapes
            count = self.__SHAPE_COUNT
            o.insert(0, 'cid', range(count, count + len(o)))

            if distance is None:
                self.__DISTANCES.append(self.DEFAULT_DISTANCE)
            else:
                self.__DISTANCES.append(distance)
            
            # Increment the number of shapes
            self.__SHAPE_COUNT += len(o)

            # Add the object for generalization
            self.__OBJECTS.append(o)

    def __reconstruct_weights(self, geomtype, **weights):
        """
        Reconstruct weights as a dict. Check whether weights are correcty set and/or apply default weights when needed
        """
        w = {}
        # If no weights are provided, set weights to default
        if len(weights) < 1:
            w = self.__DEFAULT_WEIGHTS[geomtype]
        else:
            # Retrieve allowed weights for the given geometry type
            allowed = self.__ALLOWED_WEIGHTS[geomtype]
            for name, value in weights.items():
                # If weight is not allowed, raise an error
                if name not in allowed:
                    raise Exception('{0} weight does not exists or cannot be applied to {1}.\n'.format(name, geomtype) +
                        'Available weights for {0}: {1}'.format(geomtype, ', '.join(allowed)))
                # Check if the weight is an integer
                if isinstance(value, int):
                    # If it's a negative integer other than -1, raise an error
                    if value < 0 and value != -1:
                        raise Exception('Provided weight ({0}) is negative.'.format(value))
                # Raise an error if the weight is not an integer
                else:
                    raise Exception('Provided weight ({0}) is not an integer.'.format(value))

                # Add the weight to the dict
                w[name] = value

                # Check if both stiffness and curvature are provided to raise an error
                if geomtype in ['LineString', 'Polygon']:
                    if (w.keys() >= {'curvature', 'stiffness'}):
                        raise Exception('{0} cannot have both stiffness and curvature constraints.'.format(geomtype))
        return w

            

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
        Build the observation matrice.
        """
        nb = len(points)
        self.Y = numpy.zeros(2 * nb)

        for i, p in enumerate(points):
            self.Y[2*i] = p[0]
            self.Y[2*i+1] = p[1]

    def __build_weighting_matrice(self, points, weights):
        """
        Build the weighting matrice.
        """