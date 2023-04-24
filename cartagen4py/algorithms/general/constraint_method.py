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
    def __init__(self, default_distance=5, max_iteration=1000, verbose=False):
        self.MAX_ITER = max_iteration
        self.DEFAULT_DISTANCE = default_distance
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
        
        self.__points = []
        self.__shapes = []

        # Nodes and links will be respectively populated with node to node and node to link spatial conflicts
        self.__nodes = []
        self.__links = []

        # Nodes that will accept constraints
        self.__move_points = []
        self.__stiff_points = []
        self.__curve_points = []

    def generalize(self, network_partitioning=True):
        """
        Launch the constraint generalization on the added objects.
        """
        # Checks if objects are present
        if len(self.__OBJECTS) < 1:
            raise Exception('No objects provided, cannot generalize.')

        objects = []
        for obj in self.__OBJECTS:
            single = []
            for i, o in obj.iterrows():
                single.append(o.geometry)
            objects.append(single)
        
        # Prepare points and shapes using their index
        self.__generate_point_list(objects)

        # Flag spatially conflicting points depending on the distance parameter
        self.__calculate_spatial_conflicts()

        # Build the observation matrix
        self.__build_Y()

        # Build the weighing matrice
        # self.__build_P(self.__WEIGHTS)

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

            # Check if movement weight is provided
            if 'movement' not in w.keys():
                raise Exception('{0} must have a movement constraint. A weight of -1 can be used to avoid this object from moving.'.format(geomtype))
            
            # Check if spatial weight is provided
            if 'spatial' not in w.keys():
                raise Exception('{0} must have a spatial constraint.'.format(geomtype))

            # Linestring and polygon specific verifications
            if geomtype in ['LineString', 'Polygon']:
                # Check if the movement weight is -1 (infinite)
                if w['movement'] == -1:
                    # Raise an error if a stiffness or a curvature is also provided
                    if 'stiffness' in w.keys():
                        raise Exception('Immovable {0} (movement weight equal -1) cannot have a stiffness constraint.'.format(geomtype))
                    if 'curvature' in w.keys():
                        raise Exception('Immovable {0} (movement weight equal -1) cannot have a curvature constraint.'.format(geomtype))
                else:
                    # Check if neither stiffness nor curvature is provided
                    if 'stiffness' not in w.keys() and 'curvature' not in w.keys():
                        raise Exception('Movable {0} must have a stiffness or a curvature constraint.'.format(geomtype))
                    # Check if both stiffness and curvature are provided
                    if (w.keys() >= {'curvature', 'stiffness'}):
                        raise Exception('Movable {0} cannot have both stiffness and curvature constraints.'.format(geomtype))
        return w

    def __calculate_spatial_conflicts(self):
        """
        Retrieve conflicting pairs of nodes and nodes or nodes and links.
        """

        # Check if the considered point is spatially conflicting with other points
        def check_conflicts(p, oid, sid, pid):
            point = shapely.Point(p)
            d = self.__DISTANCES[oid]
            # Loop through all objects
            for oid1, o in enumerate(self.__shapes):
                # Loop through all shapes
                for sid1, s in enumerate(o):
                    # If the shape is not the same as the concerned point
                    if oid != oid1 and sid != sid1:
                        # Retrieve the distance value associated with this object
                        d1 = self.__DISTANCES[oid1]
                        # Select the highest distance of the two
                        min_distance = max(d, d1)
                        # Retrieve the geometry of the shape
                        shape = self.__OBJECTS[oid1].geometry[sid1]
                        # Checks if that shape is within the minimum distance before retrieving pairs of nodes and links
                        if shapely.dwithin(point, shape, min_distance):
                            # Checks if the shape is a point
                            if len(s) == 1:
                                self.__nodes.append((p, (oid, sid, pid), s[0], (oid1, sid1, 0), min_distance, shapely.distance(point, shape)))
                            else:
                                # Loop through all the points of the LineString or Polygon shape
                                for pid1, p1 in enumerate(s):
                                    # Check if it's not the last point of the shape
                                    if pid1 < (len(s) - 1):
                                        # Retrieve the current and next point, calculate the line between them
                                        point1 = shapely.Point(p1)
                                        p2 = self.__shapes[oid1][sid1][pid1 + 1]
                                        point2 = shapely.Point(p2)
                                        line = shapely.LineString([point1, point2])
                                        # Calculate the distance between the point and point 1 and 2, and with the line
                                        pdist1 = shapely.distance(point, point1)
                                        pdist2 = shapely.distance(point, point2)
                                        ldist = shapely.distance(point, line)
                                        # First, checks if one of the three distance is below the minimum distance threshold
                                        if pdist1 < min_distance or pdist2 < min_distance or ldist < min_distance:
                                            # If the line distance is the smallest, insert a node to link spatial conflict
                                            if ldist < pdist1 and ldist < pdist2:
                                                self.__links.append((p, (oid, sid, pid), p1, (oid1, sid1, pid1), p2, (oid1, sid1, pid1 + 1), min_distance, ldist))
                                            else:
                                                # Determine which node is the closest
                                                if pdist1 < pdist2:
                                                    self.__nodes.append((p, (oid, sid, pid), p1, (oid1, sid1, pid1), min_distance, pdist1))
                                                else:
                                                    self.__nodes.append((p, (oid, sid, pid), p2, (oid1, sid1, pid1 + 1), min_distance, pdist2))

        # Loop through all objects
        for oid, o in enumerate(self.__shapes):
            # Loop through all shapes
            for sid, s in enumerate(o):
                # Loop through all points
                for pid, p in enumerate(s):
                    # Check conflicts with other points
                    check_conflicts(p, oid, sid, pid)
            

    def __build_Y(self):
        """
        Build the observation vector.
        """
        points = self.__points
        nodes = self.__nodes
        links = self.__links

        for oid, o in enumerate(self.__shapes):
            weights = self.__WEIGHTS[oid]

        
        size = 2 * len(points) + 2 * len(nodes) + 3 * len(links)
        
        self.__Y = numpy.zeros(size)

        for i, o in enumerate(points):
            self.__Y[(2*i)] = o[0][0]
            self.__Y[(2*i+1)] = o[0][1] 
        
        offset = 2*len(o)

        print(self.__Y)

    def __build_P(self, weights):
        """
        Build the weighting matrice.
        """
        points = list(self.__points)
        self.__P = numpy.empty()
        for i, o in enumerate(points):
            weight = weights[i]
            for w in weight:
                if weight[w] >= 0:
                    m = numpy.full()

    def __generate_point_list(self, objects):
        """
        Generate a list of points composed by a tuple of coordinates (x, y), the id of the object, the id of the shape of the object, and the id of the point inside the shape.
        Generate a nested list of objects composed of all points it is made of, sorted by shapes.
        """
        for oid, shapes in enumerate(objects):
            o = []
            for sid, shape in enumerate(shapes):
                s = []
                points = self.__get_coordinates(shape)
                for pid, p in enumerate(points):
                    self.__points.append((p, oid, sid, pid))
                    s.append(p)
                o.append(s)
            self.__shapes.append(o)


    # def __generate_points_by_shape(self, objects):
    #     """
    #     Create a list of dictionnaries where :
    #     - The key is the tuple of each point coordinates of the objects (x, y)
    #     - The value is a list of each shape index it belongs to.
    #     """
    #     for shapes in objects:
    #         points = {}
    #         for i, s in enumerate(shapes):
    #             coordinates = self.__get_coordinates(s)
    #             for p in coordinates:
    #                 if p in points:
    #                     if i not in points[p]:
    #                         points[p].append(i)
    #                 else:
    #                     points[p] = [i]
    #         self.__points.append(points)
    
    # def __generate_shapes_by_points(self):
    #     """
    #     Populate a list where each geometry is made of a list of index
    #     of the points it is made of.
    #     """
    #     points = self.__points

    #     for i, shapes in enumerate(objects):
    #         shape = []
    #         for s in shapes:
    #             points = []
    #             coordinates = self.__get_coordinates(s)
    #             for p in coordinates:
    #                 points.append(self.__get_rank(self.__points[i], p))
    #             shape.append(points)
    #         self.__shapes.append(shape)

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