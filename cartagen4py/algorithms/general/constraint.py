# This is an implementation of the least squares method proposed by Lars E. Harrie (1999)

import shapely, pprint, geopandas
import numpy as np
from cartagen4py.utils.partitioning.network import network_partition
np.set_printoptions(suppress=True)

class ConstraintMethod:
    """
    Initialize constraint method object
    Parameters
    ----------
    max_iteration : int optional
        This is the maximum number of iteration before breaking the loop. If constraints and weights are correctly set, the norm tolerance threshold should be reached before the maximum number of iteration.
        Default value is set to 1000.
    norm_tolerance : float optional
        The threshold below which the norm of the resulting point matrix is acceptable enough to break the iteration loop.
        The default value is set to 0.05.
    verbose : boolean optional
        For debugging purposes, choose to print some key values while the constraint method is calculated.
        Default set to False.
    """
    def __init__(self, max_iteration=1000, norm_tolerance=0.05, verbose=False):
        self.MAX_ITER = max_iteration
        self.NORM_TOLERANCE = norm_tolerance
        self.VERBOSE = verbose

        self.__EXPORT_SPATIAL_DEBUG = False

        self.__ALLOWED_CONSTRAINTS = {
            'Point': ['movement'],
            'LineString': ['movement', 'stiffness', 'curvature'],
            'Polygon': ['movement', 'stiffness', 'curvature']
        }

        # Create dummy arrays for distances and conflicts weights
        self.__DISTANCES = np.ndarray((0, 0))
        self.__CONFLICTS = np.ndarray((0, 0))
        
        # List of future objects added to the generalization
        self.__OBJECTS = []
        # The same list of objects once the generalization has been made
        self.__RESULTS = []
        # Stores weights as entered by the user with the add() method
        self.__WEIGHTS = []
        # Stores the number of different shapes
        self.__SHAPE_COUNT = 0
        
        # Stores all points of all objects
        self.__points = []
        # Stores all points of all objects in nested lists of shapes
        self.__shapes = []

        # Stores points influenced by specific constraints
        self.__constraints = {
            'movement': [],
            'stiffness': [],
            'curvature': []
        }

        # Nodes and links will be respectively populated with node to node and node to link spatial conflicts
        self.__nodes = []
        self.__links = []


    def add(self, *objects, **weights):
        """
        Add one or multiple geographic object for preparation before the generalization. If multiple objects are provided, they must have the same geometry type.
        Parameters
        ----------
        object : **Geopandas**, *GeoSerie*.
            One or multiple GeoSerie of geographic objects, can be points, lines or polygons (if multigeometry are provided, they will be exploded).
            If multiple objects are provided, they must be the same geometry type because the same constraints will be applied.
        weights : **int**, *optional*.
            Specify a weight value for a specific constraint. Possible weights:
            movement : **int**.
                If an object is able to move, specify the weight of the movement constraint.
            stiffness : **int**.
                If a polygon or a line is able to move but is not flexible, specify the weight of the stiffness constraint. A point cannot have a stiffness constraint.
            curvature : **int**.
                If a polygon or line is flexible, specify the weight of the curvature constraint. A point cannot have a curvature constraint.
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
            
            # Increment the number of shapes
            self.__SHAPE_COUNT += len(o)

            # Add the object for safe keeping
            self.__OBJECTS.append(o)
            # Add the object for generalization
            self.__RESULTS.append(o)

    def __reconstruct_weights(self, geomtype, **weights):
        """
        Reconstruct weights as a dict. Check whether weights are correcty set
        """
        w = {}

        # Retrieve allowed weights for the given geometry type
        allowed = self.__ALLOWED_CONSTRAINTS[geomtype]
        for name, value in weights.items():
            # If weight is not allowed, raise an error
            if name not in allowed:
                raise Exception('{0} weight does not exists or cannot be applied to {1}.\n'.format(name, geomtype) +
                    'Available weights for {0}: {1}'.format(geomtype, ', '.join(allowed)))
            # Check if the weight is an integer
            if isinstance(value, int):
                # If it's a negative integer other than -1, raise an error
                if value < 0:
                    raise Exception('Provided weight ({0}) is negative.'.format(value))
            # Raise an error if the weight is not an integer
            else:
                raise Exception('Provided weight ({0}) is not an integer.'.format(value))

            # Add the weight to the dict
            w[name] = value

        # Check if movement weight is provided and raise an error if not
        if 'movement' not in w.keys():
            raise Exception('You must provide a movement weight.')

        return w

    def add_spatial_conflicts(self, distances, spatial_weights):
        """
        Defines spatial conflict between pairs of geographic objects previously added.
        Parameters
        ----------
        distances : **Numpy**, *ndarray*.
            A numpy ndarray of dimension (x, x) where x is the number of objects added to the constraint method object before generalization.
            Value at index (i, j) must be the same as the value at index (j, i) and represents the distance between geographic objects added
            in the order they were added. 
        spatial_weights : **Numpy**, *ndarray*.
            Same as the distances matrix but the values represent the weight of the spatial conflict of the two considered objects.
        """

        # Retrieve the number of objects
        nb_obj = len(self.__OBJECTS)

        # Retrieve both array dimensions
        sd = distances.shape
        sw = spatial_weights.shape

        # Retrieve the number of row and colums of each array
        xd, yd, xw, yw = sd[0], sd[1], sw[0], sw[1]

        # Checks if both array have the same dimension
        if xd != xw or yd != yw:
            raise Exception("The distances matrix and the spatial weights matrix must have the same dimensions.")
        else:
            # Checks if matrices are square
            if xd != yd:
                raise Exception("The spatial conflict matrix must be square.")
            else:
                dimension_error = False
                # Checks if the dimensions of the matrices are the same as the number of objects
                if xd != nb_obj:
                    dimension_error = True
                if yd != nb_obj:
                    dimension_error = True
                if dimension_error:
                    raise Exception("The provided matrix must have both dimensions equal to the number of provided objects." + \
                        "Please follow the documentation to set up spatial conflicts.")

        # If nothing has failed, stores both matrices
        self.__DISTANCES = distances
        self.__CONFLICTS = spatial_weights


    def generalize(self, network_partitioning=False):
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
        if self.__EXPORT_SPATIAL_DEBUG:
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

            # Recalculate distances between conflicting nodes after the first iteration
            self.__update_distances(points)

        # Recreate objects geometries
        result = self.__reconstruct_geometries(points)

        return result

    def __build_Y(self, points):
        """
        Build the observation vector.
        """
        stiffness = self.__constraints['stiffness']
        curvature = self.__constraints['curvature']
        nodes = self.__nodes
        links = self.__links

        # Define the size of the vector
        size = 2 * len(points) +  2 * len(stiffness) + 3 * len(curvature) + len(nodes) + len(links)
    
        # Create a zero filled np array of the wanted shape
        Y = np.zeros(size)

        # Add x and y value to the begining of the vector representing the movement constraint
        for i, p in enumerate(points):
            Y[2 * i] = p[0]
            Y[2 * i + 1] = p[1]

        offset = 2 * len(points)
        for i, p in enumerate(stiffness):
            current_id = p[0]
            # Retrieve the next point in the shape
            apply, next_id = self.__get_next_point_in_shape(current_id, points)
            # The apply boolean is set to False if the current point is the last of the shape in a LineString
            if apply:
                dx, dy = self.__calculate_stiffness(current_id, next_id, points)
                # Difference between the current and the following along the x axis (x - x+1)
                Y[offset + 2 * i] = dx
                # Difference between the current and the following along the y axis (y - y+1)
                Y[offset + 2 * i + 1] = dy
        
        offset += 2 * len(stiffness)
        for i, p in enumerate(curvature):
            current_id = p[0]
            # Retrieve the previous and the next point in the shape
            apply, previous_id, next_id = self.__get_surrounding_points_in_shape(current_id, points)
            # The apply boolean is set to False if the current point is the first or the last of the shape in a LineString
            if apply:
                alpha, normp, normn = self.__calculate_curvature(previous_id, current_id, next_id, points)
                # First value of the constraint is the angle formed by the previous, the current and the next point
                Y[offset + i] = alpha
                # Second value is the norm of the vector formed by the previous and the current point
                Y[offset + i + 1] = normp
                # Third value is the norm of the vector formed by the current and the next point
                Y[offset + i + 2] = normn

        # Keep 0 for the spatial conflicts (nodes and links)

        self.__Y = Y

    def __build_B(self, points):
        """
        Build the matrix B -> Y - S(X).
        """
        stiffness = self.__constraints['stiffness']
        curvature = self.__constraints['curvature']
        nodes = self.__nodes
        links = self.__links

        S = np.zeros(2 * len(points) + 2 * len(stiffness) + 3 * len(curvature) + len(nodes) + len(links))

        for i, p in enumerate(points):
            S[2 * i] = p[0]
            S[2 * i + 1] = p[1]

        offset = 2 * len(points)
        for i, p in enumerate(stiffness):
            current_id = p[0]
            # Retrieve the next point in the shape
            apply, next_id = self.__get_next_point_in_shape(current_id, points)
            # The apply boolean is set to False if the current point is the last of the shape in a LineString
            if apply:
                dx, dy = self.__calculate_stiffness(current_id, next_id, points)
                # Difference between the current and the following along the x axis (x - x+1)
                S[offset + 2 * i] = dx
                # Difference between the current and the following along the y axis (y - y+1)
                S[offset + 2 * i + 1] = dy

        offset += 2 * len(stiffness)
        for i, p in enumerate(curvature):
            current_id = p[0]
            # Retrieve the previous and the next point in the shape
            apply, previous_id, next_id = self.__get_surrounding_points_in_shape(current_id, points)
            # The apply boolean is set to False if the current point is the first or the last of the shape in a LineString
            if apply:
                alpha, normp, normn = self.__calculate_curvature(previous_id, current_id, next_id, points)
                # First value of the constraint is the angle formed by the previous, the current and the next point
                S[offset + i] = alpha
                # Second value is the norm of the vector formed by the previous and the current point
                S[offset + i + 1] = normp
                # Third value is the norm of the vector formed by the current and the next point
                S[offset + i + 2] = normn

        offset += 3 * len(curvature)
        for i, n in enumerate(nodes):
            if n[3] > n[2]:
                # Add 0 if the actual distance is higher than the minimum distance
                S[offset + i] = 0
            else:
                # Add dist_min - dist if the actual distance is below the minimum distance
                S[offset + i] = n[2] - n[3]

        offset += len(nodes)
        for i, l in enumerate(links):
            if l[4] > l[3]:
                # Same as nodes
                S[offset + i] = 0
            else:
                # Same as nodes
                S[offset + i] = l[3] - l[4]

        return self.__Y - S

    def __build_W(self, points):
        """
        Build the weighting matrix.
        """
        constraints = []

        movement = []
        # Loop through moving points
        for i, p in enumerate(points):
            w = self.__constraints['movement'][i][1]
            # Add the weight value to the movement list two times for x and y
            movement.extend([w, w])
        # Create a full matrix 
        constraints.append(np.full(len(movement), movement))

        stiffness = []
        # Loop through the stiffness constraint
        for s in self.__constraints['stiffness']:
            stiffness.extend([s[1], s[1]])
        constraints.append(np.full(len(stiffness), stiffness))

        curvature = [] 
        # Loop through the curvature constraint
        for c in self.__constraints['curvature']:
            curvature.extend([c[1], c[1], c[1]])
        constraints.append(np.full(len(curvature), curvature))

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
        constraints.append(np.full(len(conflicts), conflicts))

        # Create the diagonal weighting matrix by concatenating all constraints matrices
        self.__W = np.diag(np.concatenate((constraints)))

    def __compute_dx(self, points):
        """
        Apply the least square method to the model.
        """
        # Build the jacobian matrix A
        A = self.__build_A(points)
        # Build the matrix B which is S - Y
        B = self.__build_B(points)

        atp = A.T @ self.__W
        atpa = atp @ A
        atpb = atp @ B

        # Solves the equation
        dx = np.linalg.solve(atpa, atpb)

        return dx

    def __build_A(self, points):
        """
        Build the jacobian matrix of the model.
        """
        # Build an identity matrix for the movement constraint
        identity = np.identity(2 * len(points))

        # Calulcate the stiffness matrix and stacking it on A
        stiffness = self.__build_stiffness(points)
        A = np.vstack((identity, stiffness))

        # Calulcate the stiffness matrix and stacking it on A
        curvature = self.__build_curvature(points)
        A = np.vstack((A, curvature))

        # Calulcate the spatial conflicts matrix and stacking it on A
        spatial = self.__build_spatial(points)
        A = np.vstack((A, spatial))

        return A

    def __build_stiffness(self, points):
        """
        Create a matrix for the stiffness constraint.
        """
        stiffness = self.__constraints['stiffness']
        m = np.zeros((2 * len(stiffness), 2 * len(points)))

        for i, p in enumerate(stiffness):
            current_id = p[0]
            # Retrieve the next point in the shape
            apply, next_id = self.__get_next_point_in_shape(current_id, points)
            # The apply boolean is set to False if the current point is the first or the last of the shape in a LineString
            if apply:
                # x for the current point
                m[2 * i][2 * current_id] = 1
                # x for the next point
                m[2 * i][2 * next_id] = -1
                # y for the current point
                m[2 * i + 1][2 * current_id + 1] = 1
                # y for the next point
                m[2 * i + 1][2 * next_id + 1] = -1

        return m

    def __build_curvature(self, points):
        """
        Create a matrix for the curvature constraint.
        """
        curvature = self.__constraints['curvature']
        m = np.zeros((3 * len(curvature), 2 * len(points)))

        for i, p in enumerate(curvature):
            p = p[0]
            # Retrieve the previous and the next point in the shape
            apply, pp, pn = self.__get_surrounding_points_in_shape(p, points)
            # The apply boolean is set to False if the current point is the first or the last of the shape in a LineString
            if apply:
                xp, yp = points[pp][0], points[pp][1]
                x, y = points[p][0], points[p][1]
                xn, yn = points[pn][0], points[pn][1]

                x_xp = x - xp
                y_yp = y - yp
                x_xn = x - xn
                y_yn = y - yn
                b = ((x_xp) * (x_xp) + (y_yp) * (y_yp))
                d = ((x_xn) * (x_xn) + (y_yn) * (y_yn))
                c = ((x_xp) * (y_yn) - (x_xn) * (y_yp))
                bd = b * d
                a = (bd) ** (-0.5)
                
                # df/dxi-1
                m[3 * i][2 * pp] = a * (-x_xp * c + y_yn * b) / b
                # df/dyi-1
                m[3 * i][2 * pp + 1] = a * (-x_xn * b - y_yp * c) / b
                # df/dxi
                m[3 * i][2 * p] = a * ((-yp + yn) * bd + c * ((x_xp) * d + (x_xn) * b)) / bd
                # df/dyi
                m[3 * i][2 * p + 1] = a * ((xp - xn) * bd + c * ((y_yp) * d + (y_yn) * b)) / bd
                # df/dxi+1
                m[3 * i][2 * pn] = a * (-x_xn * c - y_yp * d) / d
                # df/dyi+1
                m[3 * i][2 * pn + 1] = a * (x_xp * d - y_yn * c) / d

                # Calculate length of the previous and next segment
                lp = np.sqrt((x - xp) * (x - xp) + (y - yp) * (y - yp))
                ln = np.sqrt((xn - x) * (xn - x) + (yn - y) * (yn - y))

                # Influence of the length of previous segment on the current point
                # x
                m[3 * i + 1][2 * p] = (x - xp) / lp
                # y
                m[3 * i + 1][2 * p + 1] = (y - yp) / lp
                # Influence of the length of previous segment on the previous point (p)
                # xp
                m[3 * i + 1][2 * pp] = -((x - xp) / lp)
                # yp
                m[3 * i + 1][2 * pp + 1] = -((y - yp) / lp)

                # Influence of the length of following segment on the next point (n)
                # xn
                m[3 * i + 2][2 * pn] = (xn - x) / ln
                # yn
                m[3 * i + 2][2 * pn + 1] = (yn - y) / ln
                # Influence of the length of following segment on the current point
                # x
                m[3 * i + 2][2 * p] = -((xn - x) / ln)
                # y
                m[3 * i + 2][2 * p + 1] = -((yn - y) / ln)

        return m

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

            # Calculate equations factors
            a = - ((x1 - x2) / norm)
            b = - ((y1 - y2) / norm)
            c = - ((x2 - x1) / norm)
            d = - ((y2 - y1) / norm)

            # Filling the matrix
            m[i][2 * n[0]] = a
            m[i][2 * n[0] + 1] = b
            m[i][2 * n[1]] = c
            m[i][2 * n[1] + 1] = d

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
            u = - ((np.abs(a * (x0 + h) + b * y0 + c) - np.abs(a * (x0 - h) + b * y0 + c)) \
                / np.sqrt(a * a + b * b) * 2 * h)
            v = - ((np.abs(a * x0 + b * (y0 + h) + c) - np.abs(a * x0 + b * (y0 - h) + c)) \
                / np.sqrt(a * a + b * b) * 2 * h)

            w = - (1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 - y1) / (x1 + h - x2) - y1 + (x1 + h) * (y1 - y2) / (x1 + h - x2)) \
                / np.sqrt((y2 - y1) * (y2 - y1) / ((x1 + h - x2) * (x1 + h - x2)) + 1) \
                - np.abs(y0 + x0 * (y2 - y1) / (x1 - h - x2) - y1 + (x1 - h) * (y1 - y2) / (x1 - h - x2)) \
                / np.sqrt((y2 - y1) * (y2 - y1) / ((x1 - h - x2) * (x1 - h - x2)) + 1)
            ))

            d = - (1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 - (y1 + h)) / (x1 - x2) - y1 - h + x1 * (y1 + h - y2) / (x1 - x2)) \
                / np.sqrt((y2 - (y1 + h)) * (y2 - (y1 + h)) / ((x1 - x2) * (x1 - x2)) + 1) \
                - np.abs(y0 + x0 * (y2 - (y1 - h)) / (x1 - x2) - y1 + h + x1 * (y1 - h - y2) / (x1 - x2)) \
                / np.sqrt((y2 - (y1 - h)) * (y2  - (y1 - h)) / ((x1 - x2) * (x1 - x2)) + 1)
            ))

            e = - (1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 - y1) / (x1 - (x2 + h))- y1 + x1 * (y1 - y2) / (x1 - (x2 + h))) \
                / np.sqrt((y2 - y1)* (y2 - y1) / ((x1 - (x2 + h)) * (x1 - (x2 + h))) + 1) \
                - np.abs(y0 + x0 * (y2 - y1) / (x1 - (x2 - h)) - y1 + x1 * (y1 - y2) / (x1 - (x2 - h))) \
                / np.sqrt((y2 - y1) * (y2 - y1) / ((x1 - (x2 - h)) * (x1 - (x2 - h))) + 1)))

            f = - (1 / (2 * h) * (
                np.abs(y0 + x0 * (y2 + h - y1) / (x1 - x2) - y1 + x1 * (y1 - (y2 + h)) / (x1 - x2)) \
                / np.sqrt((y2 + h - y1) * (y2 + h - y1) / ((x1 - x2) * (x1 - x2)) + 1) \
                - np.abs(y0 + x0 * (y2 - h - y1) / (x1 - x2) - y1 + x1 * (y1 - (y2 - h)) / (x1 - x2)) \
                / np.sqrt((y2 - h - y1) * (y2 - y1 - h) / ((x1 - x2) * (x1 - x2)) + 1)))

            # Filling the matrix with the partial derivatives
            # For point 0, the node in the node-to-link conflict
            # x0
            m[offset + i][2 * n[0]] = u
            # y0
            m[offset + i][2 * n[0] + 1] = v
            # For point 1 of the line (1, 2)
            # x1
            m[offset + i][2 * n[1]] = w
            # y1
            m[offset + i][2 * n[1] + 1] = d
            # For point 2 of the line (1, 2)
            # x2
            m[offset + i][2 * n[2]] = e
            # y2
            m[offset + i][2 * n[2] + 1] = f

        return m

    def __calculate_stiffness(self, current_id, next_id, points):
        """
        Estimate the stiffness by calculating the amount of movement between following point on a shape.
        """
        x0, y0 = points[current_id][0], points[current_id][1]
        x1, y1 = points[next_id][0], points[next_id][1]
        return x0 - x1, y0 - y1

    def __calculate_curvature(self, previous_id, current_id, next_id, points):
        """
        Estimate the curvature by calculating the angle formed by the point, its previous and its following point.
        """
        pp = np.array((points[previous_id][0], points[previous_id][1]))
        pc = np.array((points[current_id][0], points[current_id][1]))
        pn = np.array((points[next_id][0], points[next_id][1]))

        u = (pc - pp)
        normu = np.linalg.norm(u)
        v = (pn - pc)
        normv = np.linalg.norm(v)

        return np.cross(u / normu, v / normv), normu, normv

    def __calculate_spatial_conflicts(self, points):
        """
        Retrieve conflicting pairs of nodes and nodes, or nodes and links.
        """

        # Check if the node to link conflict concern segments of different objects crossing themselves
        def line_crossing(p, point, line, points):
            # Retrieve previous and following points of the shape
            a, previous_id, next_id = self.__get_surrounding_points_in_shape(p, points)
            # Loop through both points
            for i in [previous_id, next_id]:
                # If point exists
                if i is not None:
                    # Create the geometry
                    ip = shapely.Point(points[i])
                    # Create the line formed with the current point
                    crossline = shapely.LineString([point, ip])
                    # If that line crosses the concerned line, return True
                    if shapely.crosses(line, crossline):
                        return True
            # If none of the two lines cross the concerned line, return False
            return False

        # Check if a node to node conflict already exists before adding one
        def add_node_to_node(c):
            add = True
            for n in self.__nodes:
                if (c[0] == n[0] and c[1] == n[1]):
                    add = False
                if (c[0] == n[1] and c[1] == n[0]):
                    add = False
            if add:
                self.__nodes.append(c)

        def retrieve_nodes_links(shape, geomtype, p, conflict_dist, min_dist, weight):
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
                if pdist1 < conflict_dist or pdist2 < conflict_dist or ldist < conflict_dist:
                    # If the line distance is the smallest, insert a node to link spatial conflict
                    if ldist < pdist1 and ldist < pdist2:
                        if line_crossing(p, point, line, points) == False:
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
            # Loop through all objects
            for oid1, o in enumerate(self.__shapes):
                # Loop through all shapes
                for sid1, s in enumerate(o):
                    # Checks if the shape is not the same as the concerned point
                    if oid1 == oid and sid1 == sid:
                        continue
                    else:
                        # Getting the weight of the spatial conflict constraint from the matrix
                        weight = self.__CONFLICTS[oid][oid1]
                        # Getting the distance value from the distances matrix
                        min_dist = self.__DISTANCES[oid][oid1]
                        # Setting a distance equal to 1.5 times the min distance to retrieve conflicting objects
                        # Setting conflict_dist = min_dist doesn't take the 1.5 time conflicting entities
                        conflict_dist = 1.5 * min_dist
                        
                        # Retrieve the geometry of the shape
                        shape = self.__OBJECTS[oid1].geometry[sid1]
                        # Checks if that shape is within the minimum distance before retrieving pairs of nodes and links
                        if shapely.dwithin(point, shape, conflict_dist):
                            # Stores the geometry type
                            geomtype = shape.geom_type
                            # Checks if the shape is a point
                            if geomtype == 'Point':
                                distance = shapely.distance(point, shape)
                                add_node_to_node([p, s[0], min_dist, distance, weight])
                            # If it's not, checking if it's closer to a node or a segment
                            else:
                                retrieve_nodes_links(s, geomtype, p, conflict_dist, min_dist, weight)

        # Loop through all objects
        for o in self.__shapes:
            # Loop through all shapes
            for s in o:
                # Loop through all points
                for p in s:
                    # Check conflicts with other points
                    check_conflicts(p)

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
        # Loop through each objects
        for oid, shapes in enumerate(self.__OBJECTS):
            # retrieve the weights of the concerned object
            weights = self.__WEIGHTS[oid]
            o = []
            # Loop through each shape of the object
            for sid, shape in shapes.iterrows():
                s = []
                # Retrieve the geometry type
                geomtype = shape.geometry.geom_type
                # Retrieve x, y of the points of the shape
                points = self.__get_coordinates(shape.geometry)
                # Loop through each points in the shape
                for pid, p in enumerate(points):
                    # Skipping the last point of a polygon as it is the same as the first one
                    if geomtype == 'Polygon' and pid >= (len(points) - 1):
                        continue
                    else:
                        # Append point to the full list of points
                        unique_points.append(p)
                        # Append the position of the point within its object and shape
                        self.__points.append((oid, sid, pid))
                        # Append the point to its shape
                        s.append(point_id)
                        # Add the constraint associated with the point
                        add_point_constraint(point_id, weights)
                        point_id += 1
                # Add the list of shapes to the object
                o.append(s)
            # Add the list of objects to the list
            self.__shapes.append(o)

        # Return an array of unique points
        return np.array(unique_points)

    def __get_next_point_in_shape(self, current_id, points):
        """
        Look for the next point in the shape.
        Return a boolean set to True if a point was found, i.e. the current point is not the last of a LineString.
        If the current point is the last point of a polygon, it returns the first point.
        """
        # Get the position of the current point
        position = self.__points[current_id]
        oid, sid, pid = position[0], position[1], position[2]

        # Get its shape
        shape = self.__shapes[oid][sid]
        # Get the geometry type
        geomtype = self.__OBJECTS[oid].geometry[0].geom_type
        
        apply = True
        next_id = None

        if pid < (len(shape) - 1):
            # Return the real next id if it's not the last point of the shape
            next_id = shape[pid + 1]
        else:
            if geomtype == 'Polygon':
                # If the current point id the last of the shape and it's a polygon, return the first point
                next_id = shape[0]
            else:
                # Set the boolean to false for the last point of a linestring
                apply = False
        
        return apply, next_id

    def __get_surrounding_points_in_shape(self, current_id, points):
        """
        Look for the previous and the next point in the shape.
        Return a boolean set to True if both points were found, i.e. the current point is not the first or the last of a LineString.
        If the current point is first or the last point of the polygon, it returns correct surrounding points.
        """
        # Get the position of the current point
        position = self.__points[current_id]
        oid, sid, pid = position[0], position[1], position[2]

        # Get its shape
        shape = self.__shapes[oid][sid]
        # Get the geometry type
        geomtype = self.__OBJECTS[oid].geometry[0].geom_type
        
        apply = True
        previous_id = None
        next_id = None

        # If the current point if not the first or the last
        if pid > 0 and pid < (len(shape) - 1):
            # Set previous and next as its natural previous and next point
            previous_id = shape[pid - 1]
            next_id = shape[pid + 1]
        # If the current point either the first or the last
        else:
            # Check if the shape is a polygon
            if geomtype == 'Polygon':
                # If the point is the first of the shape
                if pid == 0:
                    # Return the previous as the last point and the next as its natural next
                    previous_id = shape[-1]
                    next_id = shape[pid + 1]
                # If the point is the last of the shape
                elif pid == (len(shape) - 1):
                    # Return its natural previous and the next as the first of the shape
                    previous_id = shape[pid - 1]
                    next_id = shape[0]
            # Set boolean to False if the shape is a LineString
            else:
                apply = False
        
        return apply, previous_id, next_id

    def __update_distances(self, points):
        """
        Update the distance between pairs of nodes and pairs of nodes and links
        """
        # Loop through node to node conflicts
        for i, n in enumerate(self.__nodes):
            # Retrieve coordinates of the two nodes
            n1, n2 = points[n[0]], points[n[1]]
            x1, x2, y1, y2 = n1[0], n2[0], n1[1], n2[1]
            # Calculate the vector formed by those two points
            v = np.array([x2 - x1, y2 - y1])
            # Calculate its norm, i.e. distance and updating it
            nodedistance = np.linalg.norm(v)
            self.__nodes[i][3] = nodedistance

        # Loop through node to link conflicts
        for i, l in enumerate(self.__links):
            # Retrieve coordinates of the three nodes as array
            n0 = np.array(points[l[0]])
            n1 = np.array(points[l[1]])
            n2 = np.array(points[l[2]])
            # Calculate the distance between the first node and the line formed by the other two
            linkdistance = np.linalg.norm(np.cross(n2 - n1, n1 - n0)) / np.linalg.norm(n2 - n1)
            # Updating the distance
            self.__links[i][4] = linkdistance

    def __get_coordinates(self, shape):
        """
        Returns a list of coordinates from a shapely geometry.constraints
        """
        if shape.geom_type == 'Polygon':
            # For Polygon type
            return shape.exterior.coords
        else:
            # For LineString type
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
        """
        Exports spatial conflicts as links between pairs of nodes and pairs of nodes and links.
        This function is for DEBUGGING PURPOSES.
        """
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
        
    def get_objects_number(self):
        """
        Return the number of objects added to the generalization algorithm.
        """
        return len(self.__OBJECTS)