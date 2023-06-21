# This is an implementation of the least squares method proposed by Lars E. Harrie (1999)

import shapely, pprint, geopandas
import numpy as np
from cartagen4py.utils.partitioning.network import network_partition

class ConstraintMethod:
    """
    TODO: Remove infinite movement constraint
    Initialize constraint method object, with default constraints
    Parameters
    ----------
    default_distance : int optional
        This is the default distance for the buffers created to detect spatial conflicts.
        This value will be used with all objects added to the generalization if no distance is provided.
        
    """
    def __init__(self, default_distance=5, max_iteration=1000, norm_tolerance=0.05, verbose=False):
        self.MAX_ITER = max_iteration
        self.DEFAULT_DISTANCE = default_distance
        self.NORM_TOLERANCE = norm_tolerance
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
        self.__RESULTS = []
        self.__WEIGHTS = []
        self.__DISTANCES = []
        self.__SHAPE_COUNT = 0
        
        self.__points = []
        self.__shapes = []

        # Nodes and links will be respectively populated with node to node and node to link spatial conflicts
        self.__nodes = []
        self.__links = []

        # Nodes that will accept constraints
        self.__constraints = {
            'movement': [],
            'stiffness': [],
            'curvature': []
        }

    def generalize(self, network_partitioning=True):
        """
        Launch the constraint generalization on the added objects.
        """
        # Checks if objects are present
        if len(self.__OBJECTS) < 1:
            raise Exception('No objects provided, cannot generalize.')
        
        # Prepare points and shapes using their index
        points = self.__generate_point_list()

        # Flag spatially conflicting points depending on the distance parameter
        self.__calculate_spatial_conflicts(points)

        # Debugging tool to visualize spatial conflicts
        self.__export_spatial_conflicts(points)

        # Build the observation matrix
        self.__build_Y(points)

        # Build the weighing matrice
        self.__build_W(points)

        nb_points = len(points)
        # Calculate dx as long as the maximum iteration is not reached
        for i in range(self.MAX_ITER):
            # Apply the least square method
            dx = self.__compute_dx(points)
            # Reshape the matrix
            points += dx.reshape(nb_points, 2)
            # Calculate the norm
            norm = np.linalg.norm(dx, ord=np.inf)
            if self.VERBOSE:
                print(norm)
            # Break the loop if the norm is below the tolerance threshold
            if norm < self.NORM_TOLERANCE:
                break
            else:
                # Recalculate distances between conflicting nodes after the first iteration
                self.__update_distances(points)

        # Recreate objects geometries
        result = self.__reconstruct_geometries(points)

        return result

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
                If a polygon is able to move but is not flexible, specify the weight of the stiffness constraint. A point or line cannot have a stiffness constraint.
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

            # Add the object for safe keeping
            self.__OBJECTS.append(o)
            # Add the object for generalization
            self.__RESULTS.append(o)

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

            # # Linestring and polygon specific verifications
            # if geomtype in ['LineString', 'Polygon']:
            #     # Check if the movement weight is -1 (infinite)
            #     if w['movement'] == -1:
            #         # Raise an error if a stiffness or a curvature is also provided
            #         if 'stiffness' in w.keys():
            #             raise Exception('Immovable {0} (movement weight equal -1) cannot have a stiffness constraint.'.format(geomtype))
            #         if 'curvature' in w.keys():
            #             raise Exception('Immovable {0} (movement weight equal -1) cannot have a curvature constraint.'.format(geomtype))
            #     else:
            #         # Check if neither stiffness nor curvature is provided
            #         if 'stiffness' not in w.keys() and 'curvature' not in w.keys():
            #             raise Exception('Movable {0} must have a stiffness or a curvature constraint.'.format(geomtype))
            #         # Check if both stiffness and curvature are provided
            #         if (w.keys() >= {'curvature', 'stiffness'}):
            #             raise Exception('Movable {0} cannot have both stiffness and curvature constraints.'.format(geomtype))
        return w

    def __build_Y(self, points):
        """
        Build the observation vector.
        """
        stiffness = self.__constraints['stiffness']
        nodes = self.__nodes
        links = self.__links

        # Define the size of the vector
        size = 2 * len(points) + 2 * len(stiffness) + len(nodes) + len(links)
        # Create a zero filled np array of the wanted shape
        Y = np.zeros(size)

        # Add x and y value to the begining of the vector representing the movement constraint
        for i, p in enumerate(points):
            Y[2 * i] = p[0]
            Y[2 * i + 1] = p[1]

        offset = 2 * len(points)
        for i, p in enumerate(stiffness):
            current, dx, dy = self.__calculate_stiffness(p, points)
            Y[offset + 2 * i] = dx
            Y[offset + 2 * i + 1] = dy

        offset += 2 * len(stiffness)
        # Add minimum distance between close nodes
        for i, n in enumerate(nodes):
            Y[offset + i] = n[2]

        # Add minimum distance between close nodes and links
        offset += len(nodes)
        for i, l in enumerate(links):
            Y[offset + i] = l[3]

        self.__Y = Y

    def __build_B(self, points):
        """
        Build the matrix B -> Y - S(X).
        """
        stiffness = self.__constraints['stiffness']
        nodes = self.__nodes
        links = self.__links

        S = np.zeros(2 * len(points) + 2 * len(stiffness) + len(nodes) + len(links))

        for i, p in enumerate(points):
            S[2 * i] = p[0]
            S[2 * i + 1] = p[1]

        offset = 2 * len(points)
        for i, p in enumerate(stiffness):
            current, dx, dy = self.__calculate_stiffness(p, points)
            S[offset + 2 * i] = dx
            S[offset + 2 * i + 1] = dy

        offset += 2 * len(stiffness)
        for i, n in enumerate(nodes):
            S[offset + i] = n[3]

        offset += len(nodes)
        for i, l in enumerate(links):
            S[offset + i] = l[4]

        return self.__Y - S

    def __build_W(self, points):
        """
        Build the weighting matrix.
        """
        movement = []
        # Loop through moving points
        for i, p in enumerate(points):
            w = self.__constraints['movement'][i][1]
            # Add the weight value to the movement list two times for x and y
            movement.extend([w, w])
        # Create a full matrix 
        m = np.full(len(movement), movement)

        constraints = []
        # Loop through the non movement constraint
        for cname, clist in self.__constraints.items():
            # Check if the constraint is not movement and that some points are affected
            if cname != 'movement' and len(clist) > 0:
                p = []
                # Loop through each point concerned by the constraint
                for c in clist:
                    p.extend([c[1], c[1]])
                # Create the matrix with the weights
                constraints.append(np.full(len(p), p))

        conflicts = []
        # Loop through node to node conflicts
        for n in self.__nodes:
            # Add the weight to the list
            conflicts.append(n[4])
        # Loop through node to link conflicts
        for l in self.__links:
            # Add the weight to the list
            conflicts.append(l[5])
        # Create the spatial conflict weight matrix
        spatial = np.full(len(conflicts), conflicts)

        # Create the diagonal weighting matrix by concatenating movement matrix
        # with others constraints and spatial conflicts matrices
        self.__W = np.diag(np.concatenate((m, *constraints, spatial)))

    def __compute_dx(self, points):
        """
        Apply the least square method to the model.
        """
        A = self.__build_A(points)
        B = self.__build_B(points)

        atp = A.T @ self.__W
        atpa = atp @ A
        atpb = atp @ B

        dx = np.linalg.solve(atpa, atpb)
        return dx

    def __build_A(self, points):
        """
        Build the jacobian matrix of the model.
        """
        identity = np.identity(2 * len(points))
        stiffness = self.__build_stiffness(points)
        A = np.vstack((identity, stiffness))
        spatial = self.__build_spatial(points)
        A = np.vstack((A, spatial))
        return A

    def __build_stiffness(self, points):
        """
        Create a matrix for the stiffness constraint.
        """
        stiff_points = self.__constraints['stiffness']
        shapes = self.__shapes
        m = np.zeros((2 * len(stiff_points), 2 * len(points)))

        for i, p in enumerate(stiff_points):
            current, dx, dy = self.__calculate_stiffness(p, points)
            m[2 * i][2 * current] = dx
            m[2 * i + 1][2 * current + 1] = dy

        return m


    def __build_curvature(self):
        """
        Create a matrix for the curvature constraint.
        """

    def __build_spatial(self, points):
        """
        Create a matrix for the spatial conflicts constraint.
        """
        nodes = self.__nodes
        links = self.__links

        # Create the matrix
        m = np.zeros((len(nodes) + len(links), 2 * len(points)))

        # Loop through all node to node conflicts
        for i, n in enumerate(nodes):
            # Retrieve both points coordinates
            n1, n2 = points[n[0]], points[n[1]]
            x1, x2, y1, y2 = n1[0], n2[0], n1[1], n2[1]
            # Calculate the vector formed by those two points
            v = np.array([x2 - x1, y2 - y1])
            # Calculate the vector norm
            norm = np.linalg.norm(v)
            # Filling the matrix with the influence of x and y on the vector norm
            xi = np.abs(v[0]) / norm
            yi = np.abs(v[1]) / norm
            m[i][n[0]] = xi
            m[i][n[0] + 1] = yi
            m[i][n[1]] = xi
            m[i][n[1] + 1] = yi

        offset = len(nodes)
        # Loop through all node to link conflicts
        for i, n in enumerate(links):
            # Retrieve the three points coordinates
            n0, n1, n2 = points[n[0]], points[n[1]], points[n[2]]
            x0, y0 = n0[0], n0[1]
            x1, y1 = n1[0], n1[1]
            x2, y2 = n2[0], n2[1]

            # Calculate the line formed by n1 and n2
            a = (y2 - y1) / (x1 - x2)
            b = 1.0
            c = x1 * (y1 - y2) / (x1 - x2) - y1

            # Calculate the equation factors
            h = 0.001
            u = (np.abs(a * (x0 + h) + b * y0 + c) - np.abs(a * (x0 - h) + b * y0 + c)) \
                / np.sqrt(np.square(a) + np.square(b)) * 2 * h
            v = (np.abs(a * x0 + b * (y0 + h) + c) - np.abs(a * x0 + b * (y0 - h) + c)) \
                / np.sqrt(np.square(a) + np.square(b)) * 2 * h
            w = 1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 - y1) / (x1 + h - x2) - y1 + (x1 + h) * (y1 - y2) / (x1 + h - x2)) \
                / np.sqrt(np.square(y2 - y1) / np.square(x1 + h - x2) + 1) \
                - np.abs(y0 + x0 * (y2 - y1) / (x1 - h - x2) - y1 + (x1 - h) * (y1 - y2) / (x1 - h - x2)) \
                / np.sqrt(np.square(y2 - y1) / np.square(x1 - h - x2) + 1)
            )
            d = 1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 - (y1 + h)) / (x1 - x2) - y1 - h + x1 * (y1 + h - y2) / (x1 - x2)) \
                / np.sqrt(np.square(y2 - (y1 + h)) / np.square(x1 - x2) + 1) \
                - np.abs(y0 + x0 * (y2 - (y1 - h)) / (x1 - x2) - y1 + h + x1 * (y1 - h - y2) / (x1 - x2)) \
                / np.sqrt(np.square(y2 - (y1 - h)) / np.square(x1 - x2) + 1)
            )
            e = 1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 - y1) / (x1 - (x2 + h)) - y1 + x1 * (y1 - y2) / (x1 - (x2 + h))) \
                / np.sqrt(np.square(y2 - y1) / np.square(x1 - (x2 + h)) + 1) \
                - np.abs(y0 + x0 * (y2 - y1) / (x1 - (x2 - h)) - y1 + x1 * (y1 - y2) / (x1 - (x2 - h))) \
                / np.sqrt(np.square(y2 - y1) / np.square(x1 - (x2 - h)) + 1)
            )
            f = 1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 + h - y1) / (x1 - x2) - y1 + x1 * (y1 - (y2 + h)) / (x1 - x2)) \
                / np.sqrt(np.square(y2 + h - y1) / np.square(x1 - x2) + 1) \
                - np.abs(y0 + x0 * (y2 - h - y1) / (x1 - x2) - y1 + x1 * (y1 - (y2 - h)) / (x1 - x2)) \
                / np.sqrt(np.square(y2 - h - y1) / np.square(x1 - x2) + 1)
            )

            # Filling the matrix with the partial derivatives
            m[offset + i][n[0]] = u
            m[offset + i][n[0] + 1] = v
            m[offset + i][n[1]] = w
            m[offset + i][n[1] + 1] = d
            m[offset + i][n[2]] = e
            m[offset + i][n[2] + 1] = f

        return m

    def __generate_point_list(self):
        """
        Generate a list of points composed by a tuple of coordinates (x, y), the id of the object, the id of the shape of the object, and the id of the point inside the shape.
        Generate a nested list of objects composed of all points it is made of, sorted by shapes.
        Generate lists of points that will accept constraints
        """

        def add_point_constraint(pid, weights):
            for cname, clist in self.__constraints.items():
                if cname in weights.keys():
                    w = np.inf if weights[cname] == -1 else weights[cname]
                    self.__constraints[cname].append((pid, w))

        unique_points = []
        point_id = 0
        for oid, shapes in enumerate(self.__OBJECTS):
            weights = self.__WEIGHTS[oid]
            o = []
            for sid, shape in shapes.iterrows():
                s = []
                geomtype = shape.geometry.geom_type
                points = self.__get_coordinates(shape.geometry)
                for pid, p in enumerate(points):
                    if geomtype == 'Polygon' and pid >= (len(points) - 1):
                        continue
                    else:
                        unique_points.append(p)
                        self.__points.append((oid, sid, pid))
                        s.append(point_id)
                        add_point_constraint(point_id, weights)
                        point_id += 1
                o.append(s)
            self.__shapes.append(o)

        return np.array(unique_points)

    def __calculate_spatial_conflicts(self, points):
        """
        Retrieve conflicting pairs of nodes and nodes, or nodes and links.
        """

        def add_node_to_node(c):
            add = True
            for n in self.__nodes:
                if (c[0] == n[0] and c[1] == n[1]) or (c[0] == n[1] and c[1] == n[0]):
                    add = False
            if add:
                self.__nodes.append(c)

        def retrieve_nodes_links(shape, geomtype, p, min_dist, weight):
            point = shapely.Point(points[p])
            # Loop through all the points of the LineString or Polygon shape
            for pid1, p1 in enumerate(shape):
                notlast = True
                # Check if it's not the last point of the shape
                if pid1 < (len(shape) - 1):
                    p2 = shape[pid1 + 1]
                else:
                    notlast = False
                    if geomtype == 'Polygon':
                        p2 = shape[0]
                    else:
                        continue
                # Retrieve the current and next point, calculate the line between them
                point1 = shapely.Point(points[p1])
                point2 = shapely.Point(points[p2])
                line = shapely.LineString([point1, point2])
                # Calculate the distance between the point and point 1 and 2, and with the line
                pdist1 = shapely.distance(point, point1)
                pdist2 = shapely.distance(point, point2)
                ldist = shapely.distance(point, line)
                # First, re-checks if one of the three distance is below the minimum distance threshold
                if pdist1 < min_dist or pdist2 < min_dist or ldist < min_dist:
                    # If the line distance is the smallest, insert a node to link spatial conflict
                    if ldist < pdist1 and ldist < pdist2:
                        self.__links.append([p, p1, p2, min_dist, ldist, weight])
                    else:
                        # Determine which node is the closest
                        if pdist1 < pdist2:
                            add_node_to_node([p, p1, min_dist, pdist1, weight])
                        else:
                            if notlast:
                                add_node_to_node([p, p2, min_dist, pdist2, weight])

        # Check if the considered point is spatially conflicting with other points
        def check_conflicts(p):
            point = shapely.Point(points[p])
            oid, sid, pid = self.__points[p][0], self.__points[p][1], self.__points[p][2]
            dmin = self.__DISTANCES[oid]
            weight = self.__WEIGHTS[oid]['spatial']
            # Loop through all objects
            for oid1, o in enumerate(self.__shapes):
                # Loop through all shapes
                for sid1, s in enumerate(o):
                    # Checks if the shape is not the same as the concerned point
                    if oid1 == oid and sid1 == sid:
                        continue
                    else:
                        # Retrieve the geometry of the shape
                        shape = self.__OBJECTS[oid1].geometry[sid1]
                        # Checks if that shape is within the minimum distance before retrieving pairs of nodes and links
                        if shapely.dwithin(point, shape, dmin):
                            # Stores the geometry type
                            geomtype = shape.geom_type
                            # Checks if the shape is a point
                            if geomtype == 'Point':
                                distance = shapely.distance(point, shape)
                                self.__nodes.append([p, s[0], dmin, distance, weight])
                            # If it's not, checking if it's closer to a node or a segment
                            else:
                                retrieve_nodes_links(s, geomtype, p, dmin, weight)

        # Loop through all objects
        for o in self.__shapes:
            # Loop through all shapes
            for s in o:
                # Loop through all points
                for p in s:
                    # Check conflicts with other points
                    check_conflicts(p)

    def __calculate_stiffness(self, point, points):
        """
        Estimate the stiffness by calculating the amount of movement between following point on a shape.
        TODO: Only working for buildings
        """
        current_id = point[0]
        position = self.__points[current_id]
        oid, sid, pid = position[0], position[1], position[2]
        shape = self.__shapes[oid][sid]
        pid1 = pid + 1
        x0, y0 = points[current_id][0], points[current_id][1]

        if pid < (len(shape) - 1):
            pid1 = pid + 1
        else:
            pid1 = 0

        next_id = shape[pid1] 
        x1, y1 = points[next_id][0], points[next_id][1]
        dx = np.abs(x1 - x0)
        dy = np.abs(y1 - y0)

        return current_id, dx, dy

    def __update_distances(self, points):
        """
        Update the distance between pairs of nodes and pairs of nodes and links
        """
        for i, n in enumerate(self.__nodes):
            n1, n2 = points[n[0]], points[n[1]]
            x1, x2, y1, y2 = n1[0], n2[0], n1[1], n2[1]
            v = np.array([x2 - x1, y2 - y1])
            norm = np.linalg.norm(v)
            self.__nodes[i][3] = norm
            #print('{0} to {1} = {2}'.format(n[0], n[1], norm))

        for i, n in enumerate(self.__links):
            n0 = np.array(points[n[0]])
            n1 = np.array(points[n[1]])
            n2 = np.array(points[n[2]])
            d = np.linalg.norm(np.cross(n2 - n1, n1 - n0)) / np.linalg.norm(n2 - n1)
            self.__links[i][4] = d
            #print('{0} to {1} and {2} = {3}'.format(n[0], n[1], n[2], d))

    def __get_coordinates(self, shape):
        """
        Returns a list of coordinates from a shapely geometry.constraints
        """
        if shape.geom_type == 'Polygon':
            return shape.exterior.coords
        else:
            return shape.coords

    def __reconstruct_geometries(self, points):
        """
        Reconstruct Geoseries with the generalized geometries
        """
        # Loop through the results objects
        for oid, o in enumerate(self.__RESULTS):
            # Loop through the shapes of each objects
            for sid, s in o.iterrows():
                # Retrieve geomtype and shape from the list of points
                geomtype = s.geometry.geom_type
                shape = self.__shapes[oid][sid]
                # If the shape is a simple point, updating it
                if geomtype == 'Point':
                    self.__RESULTS[oid].loc[sid, 'geometry'] = shapely.Point(points[shape[0]])
                else:
                    # If the shape is a line or polygon
                    sha = []
                    # Loop through all the points of the shape and append them to a list
                    for p in shape:
                        sha.append(shapely.Point(points[p]))
                    # Create the geometry from the list depending on the geometry type
                    if geomtype == 'LineString':
                        geometry = shapely.LineString(sha)
                    else:
                        # If the geometry is a polygon, add the first point as the last to close it
                        sha.append(sha[0])
                        geometry = shapely.Polygon(sha)
                    # Update the geometry of the shape
                    self.__RESULTS[oid].loc[sid, 'geometry'] = geometry
        return self.__OBJECTS

    def __export_spatial_conflicts(self, points):
        pointslist = []
        shape_id = 0
        for o in self.__shapes:
            for s in o:
                for p in s:
                    pointslist.append({
                    'id' : p,
                    'object' : shape_id,
                    'geometry' : shapely.Point(points[p])
                }) 
                shape_id += 1         

        points_gdf = geopandas.GeoDataFrame(pointslist, crs=3857)

        nodes = []
        for i, node in enumerate(self.__nodes):
            nodes.append({
                'nid0': node[0],
                'nid1': node[1],
                'group': i,
                'geometry': shapely.LineString([points[node[0]], points[node[1]]])
            })

        if len(nodes) > 0:
            nodes_gdf = geopandas.GeoDataFrame(nodes, crs=3857)
            nodes_gdf.to_file("cartagen4py/data/data_bourbonnaise/nodes.geojson", driver="GeoJSON")

        links = []
        for lid, link in enumerate(self.__links):
            p0 = shapely.Point(points[link[0]])
            p1 = shapely.Point(points[link[1]])
            p2 = shapely.Point(points[link[2]])
            line = shapely.LineString([p1, p2])
            dist = line.project(p0)
            projected = line.interpolate(dist)
            links.append({
                'lid': lid,
                'geometry': shapely.LineString([p0, projected])
            })

        if len(links) > 0:
            links_gdf = geopandas.GeoDataFrame(links, crs=3857)
            links_gdf.to_file("cartagen4py/data/data_bourbonnaise/links.geojson", driver="GeoJSON")
        
        points_gdf.to_file("cartagen4py/data/data_bourbonnaise/points_nodes.geojson", driver="GeoJSON")
        
        