from collections import Counter
import shapely
import geopandas as gpd
import numpy as np
from cartagen.utils.partitioning.network import partition_networks

class LeastSquaresMethod:
    """
    Generalise cartographic features using the method of least squares.

    This method was proposed by :footcite:p:`harrie:1999` and uses the method
    of least squares to move the vertexes of every features (points, lines and polygons)
    keeping the provided constraints as satisfied as possible.

    Parameters
    ----------
    max_iteration : int, optional
        This is the maximum number of iteration before breaking the loop. If constraints and weights are correctly set,
        the norm tolerance threshold should be reached before the maximum number of iteration.
    norm_tolerance : float, optional
        The threshold below which the norm of the resulting matrix is acceptable enough to break the iteration loop.
    intra_conflicts : bool, optional
        If set to True, the features of the same GeoDataFrame can be spatially conflicting.
    crossing_node_threshold : int, optional
        When two vertexes from different line intersects, define the number of vertexes connected on both sides from the
        intersection on wich spatial conflicts between those two lines doesn't apply. Setting this value to 0 will not
        preserve line intersections.

    References
    ----------
    .. footbibliography::

    Notes
    -----
    A more in-depth overview of this method of generalisation is available inside
    the user manual of the documentation.
    """
    def __init__(self, max_iteration=1000, norm_tolerance=0.05, intra_conflicts=True, crossing_node_threshold=7):
        self.MAX_ITER = max_iteration
        self.NORM_TOLERANCE = norm_tolerance
        self.INTRA_CONFLICTS = intra_conflicts
        self.CROSSING_NODE_THRESHOLD = crossing_node_threshold

        self.__ALLOWED_CONSTRAINTS = {
            'Point': ['movement'],
            'LineString': ['movement', 'stiffness', 'curvature'],
            'Polygon': ['movement', 'stiffness', 'curvature']
        }

        # Create dummy arrays for distances and conflicts weights
        self.distances = np.ndarray((0, 0))
        self.spatial_weights = np.ndarray((0, 0))
        
        # List of future objects added to the generalisation
        self.__OBJECTS = []
        # The same list of objects once the generalisation has been made
        self.__RESULTS = []
        # Stores weights as entered by the user with the add() method
        self.__WEIGHTS = []
        # Stores the number of different shapes
        self.__SHAPE_COUNT = 0
        
        # Stores all points of all objects in nested lists of shapes
        self.__shapes = []
        self.__points = []

        # Stores points influenced by specific constraints
        self.__constraints = {
            'movement': [],
            'stiffness': [],
            'curvature': [],
            # Nodes and links will be respectively populated with node to node and node to link spatial conflicts
            'nodes': [],
            'links': []
        }

    def add(self, *objects, **weights):
        """
        Add features before generalisation. This method adds one or multiple objects
        to the generalisation along with their associated weights.
        If multiple objects are provided, they must have the same geometry type.

        Parameters
        ----------
        object : GeoDataFrame
            One or multiple GeoDataFrame of geographic objects, can be Point, LineString or Polygon (if MultiGeometry are provided, they will be exploded).
            If multiple objects are provided, they must be the same geometry type because the same constraints will be applied.
        weights : int, optional
            Specify a constraint with its weight. Possible weights:

            - movement : Specify the weight of the movement constraint.
            - stiffness : Specify the weight of the stiffness constraint.
            - curvature : Specify the weight of the curvature constraint.

            .. list-table::
                :widths: 20 30 50
                :header-rows: 1

                * - Constraint
                  - Geometry type
                  - Impact
                * - movement
                  - Point, LineString, Polygon
                  - The object should move as little as possible
                * - stiffness
                  - LineString, Polygon
                  - The internal geometry should be invariant, `i.e.` the vertexes movement within the same object will try not to move closer or away from each other.
                * - curvature
                  - LineString, Polygon
                  - The curvature of a line or a polygon border should not change, `i.e.` the angle formed by two connected segments will try not to change.
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
            # Add the object for generalisation
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
        Defines spatial conflict between pairs of geographic objects.

        Once geographic objects has been added, this method is used to set up
        spatial conflicts between each of them.

        Parameters
        ----------
        distances : ndarray
            A numpy ndarray of dimension (n, n) where n is the number of
            objects added to the constraint method in
            the same order they has been added.
            For example, if 3 objects are provided, a numpy ndarray
            of a (3, 3) shape must be provided, such as:

            .. code-block:: Python

                [[25. 28. 30.]
                 [28. 15. 20.]
                 [30. 20. 18.]]

        spatial_weights : ndarray
            Same as the distances matrix but the values represent the weight of the spatial conflict of the two objects.
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
        self.distances = distances
        self.spatial_weights = spatial_weights

    def generalize(self, *network_partitioning):
        """
        Launch the constraint generalisation on the objects.

        Parameters
        ----------
        network_partitioning : GeoDataFrame of LineString
            One or multiple line that represents the network used
            to partition the data.
        
        Returns
        -------
        generalized : list of GeoDataFrame
            A list containing the generalized objects.
            This list contains the same number of provided objects in order.
        """
        # Checks if objects are present
        if len(self.__OBJECTS) < 1:
            raise Exception('No objects provided, cannot generalize.')
        
        # Prepare geometries before launching the constraint method
        points = self.__prepare_geometries()

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
            # if self.VERBOSE:
            #     print(norm)
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
        movement = self.__constraints['movement']
        stiffness = self.__constraints['stiffness']
        curvature = self.__constraints['curvature']
        nodes = self.__constraints['nodes']
        links = self.__constraints['links']

        # Define the size of the vector
        size = 2 * len(movement) +  2 * len(stiffness) + 3 * len(curvature) + len(nodes) + len(links)
    
        # Create a zero filled np array of the wanted shape
        Y = np.zeros(size)

        # Add x and y value to the begining of the vector representing the movement constraint
        for i, p in enumerate(movement):
            coordinates = points[p[0]]
            Y[2 * i] = coordinates[0]
            Y[2 * i + 1] = coordinates[1]

        offset = 2 * len(movement)
        for i, p in enumerate(stiffness):
            dx, dy = self.__calculate_stiffness(p, points)
            # Difference between the current and the following along the x axis (x - x+1)
            Y[offset + 2 * i] = dx
            # Difference between the current and the following along the y axis (y - y+1)
            Y[offset + 2 * i + 1] = dy
        
        offset += 2 * len(stiffness)
        for i, p in enumerate(curvature):
            alpha, normp, normn = self.__calculate_curvature(p, points)
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
        movement = self.__constraints['movement']
        stiffness = self.__constraints['stiffness']
        curvature = self.__constraints['curvature']
        nodes = self.__constraints['nodes']
        links = self.__constraints['links']

        S = np.zeros(2 * len(movement) + 2 * len(stiffness) + 3 * len(curvature) + len(nodes) + len(links))

        for i, p in enumerate(movement):
            coordinates = points[p[0]]
            S[2 * i] = coordinates[0]
            S[2 * i + 1] = coordinates[1]

        offset = 2 * len(movement)
        for i, p in enumerate(stiffness):
            dx, dy = self.__calculate_stiffness(p, points)
            # Difference between the current and the following along the x axis (x - x+1)
            S[offset + 2 * i] = dx
            # Difference between the current and the following along the y axis (y - y+1)
            S[offset + 2 * i + 1] = dy

        offset += 2 * len(stiffness)
        for i, p in enumerate(curvature):
            alpha, normp, normn = self.__calculate_curvature(p, points)
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
        # Stores arrays inside a list
        constraints = []

        movement = []
        # Loop through moving points
        for m in self.__constraints['movement']:
            w = m[1]
            # Add the weight value to the movement list two times for x and y
            movement.extend([w, w])
        constraints.append(np.full(len(movement), movement))

        stiffness = []
        # Loop through stiffness constrained points
        for s in self.__constraints['stiffness']:
            w = s[2]
            stiffness.extend([w, w])
        constraints.append(np.full(len(stiffness), stiffness))

        curvature = [] 
        # Loop through curvature constrained points
        for c in self.__constraints['curvature']:
            w = c[3]
            curvature.extend([w, w, w])
        constraints.append(np.full(len(curvature), curvature))

        spatial = []
        # Loop through node to node conflicts
        for n in self.__constraints['nodes']:
            w = n[4]
            # Add the weight to the list
            spatial.append(w)
        # Loop through node to link conflicts
        for l in self.__constraints['links']:
            w = l[5]
            # Add the weight to the list
            spatial.append(w)
        # Create the spatial conflict weight matrix
        constraints.append(np.full(len(spatial), spatial))

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

        # Create the matrix
        m = np.zeros((2 * len(stiffness), 2 * len(points)))

        for i, s in enumerate(stiffness):
            # x for the current point
            m[2 * i][2 * s[0]] = 1
            # x for the next point
            m[2 * i][2 * s[1]] = -1
            # y for the current point
            m[2 * i + 1][2 * s[0] + 1] = 1
            # y for the next point
            m[2 * i + 1][2 * s[1] + 1] = -1

        return m

    def __build_curvature(self, points):
        """
        Create a matrix for the curvature constraint.
        """
        curvature = self.__constraints['curvature']

        # Create the matrix
        m = np.zeros((3 * len(curvature), 2 * len(points)))

        for i, c in enumerate(curvature):
            pp, p, pn = c[0], c[1], c[2]

            # Coordinates of the previous point
            xp, yp = points[pp][0], points[pp][1]
            # Coordinates of the current point
            x, y = points[p][0], points[p][1]
            # Coordinates of the next point
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
        nodes = self.__constraints['nodes']
        links = self.__constraints['links']

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

    def __calculate_stiffness(self, stiffness, points):
        """
        Estimate the stiffness by calculating the amount of movement between following point on a shape.
        """
        current_id, next_id = stiffness[0], stiffness[1]

        x0, y0 = points[current_id][0], points[current_id][1]
        x1, y1 = points[next_id][0], points[next_id][1]

        return x0 - x1, y0 - y1

    def __calculate_curvature(self, curvature, points):
        """
        Estimate the curvature by calculating the angle formed by the point, its previous and its following point.
        """
        previous_id, current_id, next_id = curvature[0], curvature[1], curvature[2]

        pp = np.array((points[previous_id][0], points[previous_id][1]))
        pc = np.array((points[current_id][0], points[current_id][1]))
        pn = np.array((points[next_id][0], points[next_id][1]))

        u = (pc - pp)
        normu = np.linalg.norm(u)
        v = (pn - pc)
        normv = np.linalg.norm(v)

        return np.cross(u / normu, v / normv), normu, normv

    def __prepare_geometries(self):
        """
        Prepare geometries and stores them on adequate shapes to be handle by the constraint method.
        """
        # Prepare points and shapes using their index
        points = self.__generate_point_list()

        # Setup the constraints
        self.__setup_constraints(points)

        # Flag spatially conflicting points
        self.__calculate_spatial_conflicts(points)

        return np.array(points)

    def __generate_point_list(self):
        """
        Generate a list of points composed by a tuple of coordinates (x, y), the id of the object, the id of the shape of the object, and the id of the point inside the shape.
        Generate a nested list of objects composed of all points it is made of, sorted by shapes.
        Generate lists of points that will accept constraints
        """
        unique_points = []

        def get_coordinates(self, shape):
            """
            Returns a list of coordinates from a shapely geometry.constraints
            """
            if shape.geom_type == 'Polygon':
                # For Polygon type
                return shape.exterior.coords
            else:
                # For LineString type
                return shape.coords

        def get_pid_if_duplicate(coordinates):
            """
            Return True or False if those coordinates already exists.
            If they exist, return the index of the already existing point.
            """
            for pid, c in enumerate(unique_points):
                if c == coordinates:
                    return True, pid
            return False, None

        point_id = 0
        # Loop through each objects
        for oid, shapes in enumerate(self.__OBJECTS):
            o = []
            # Loop through each shape of the object
            for sid, shape in shapes.iterrows():
                geom = shape.geometry

                s = []
                # Retrieve the geometry type
                geomtype = geom.geom_type

                points = []
                # Retrieve x, y of the points of the shape
                if geomtype == 'Polygon':
                    # For Polygon type
                    points = geom.exterior.coords
                else:
                    # For LineString type
                    points = geom.coords

                # Retrieve first and last index of the shape
                enclosing = [0, len(points) - 1]

                # Loop through each points in the shape
                for pid, p in enumerate(points):
                    # Skipping the last point of a polygon as it is the same as the first one
                    if geomtype == 'Polygon' and pid == enclosing[1]:
                        continue

                    # Stores the index of the point to apply constraints
                    cid = point_id

                    duplicate, dpid = get_pid_if_duplicate(p)
                    if duplicate:
                        # Append the already existing point index to its shape
                        s.append(dpid)
                        # Change the point index to the already existing one to add constraints
                        cid = dpid
                    else:
                        # Append the point index to its shape
                        s.append(point_id)
                        # Append the coordinates to the full list of points
                        unique_points.append(p)
                        # Increment the point index if it's a new point
                        point_id += 1
                    
                # Add the list of shapes to the object
                o.append(s)
            # Add the list of objects to the list
            self.__shapes.append(o)

        # Stores unique points
        self.__points = unique_points

        # Return an array of unique points
        return np.array(unique_points)

    def __setup_constraints(self, points):
        """
        Create constraint and their related properties to easily create matrices afterwards.
        """

        def add_movement(obj, value):
            """
            Add the movement constraint as a list with :
                - Index of the point
                - Weight of the constraint
            If the point as already a movement constraint, it replaces the current weight if the new is higher.
            """
            for s in obj:
                for pid in s:
                    add = True
                    for eid, existing in enumerate(self.__constraints['movement']):
                        if existing[0] == pid:
                            self.__constraints['movement'][eid][1] = max(value, existing[1])
                            add = False
                    if add:
                        self.__constraints['movement'].append([pid, value])
        
        def add_stiffness(obj, value, geomtype):
            """
            Add the stiffness constraint as a list with :
                - Index of the current point
                - Index of the next point in the shape
                - Weight of the constraint
            """
            for shape in obj:
                for pid, point in enumerate(shape):
                    next_id = None
                    # If it's not the last point of the shape
                    if pid < (len(shape) - 1):
                        # Return the real next id
                        next_id = shape[pid + 1]
                    # If the current point id is the last of the shape
                    else:
                        # If it's a polygon, return the first point
                        if geomtype == 'Polygon':
                            next_id = shape[0]

                    if next_id is not None:
                        self.__constraints['stiffness'].append([point, next_id, value])

        def add_curvature(obj, value, geomtype):
            """
            Add the curvature constraint as a list with :
                - Index of the previous point in the shape
                - Index of the current point
                - Index of the next point in the shape
                - Weight of the constraint
            """
            for shape in obj:
                for pid, point in enumerate(shape):
                    previous_id = None
                    next_id = None
                    # If the current point if not the first or the last
                    if pid > 0 and pid < (len(shape) - 1):
                        # Set previous and next as its natural previous and next point
                        previous_id = shape[pid - 1]
                        next_id = shape[pid + 1]
                    # If the current point is either the first or the last
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
                    if previous_id is not None and next_id is not None:
                        self.__constraints['curvature'].append([previous_id, point, next_id, value])
                            
        for oid, o in enumerate(self.__shapes):
            # Get the geometry type
            geomtype = self.__OBJECTS[oid].geometry[0].geom_type
            # Get the constraints weights
            weights = self.__WEIGHTS[oid]
            for name, value in weights.items():
                if name == 'movement':
                    add_movement(o, value)
                elif name == 'stiffness':
                    add_stiffness(o, value, geomtype)
                elif name == 'curvature':
                    add_curvature(o, value, geomtype)

    def __calculate_spatial_conflicts(self, points):
        """
        Retrieve conflicting pairs of nodes and nodes, or nodes and links.
        """

        node_crossing = []
        edge_crossing = []

        def remove_conflicts(remove, line1, line2):
            def remove_node_to_node(n1, n2):
                indexes = []
                for i, c in enumerate(self.__constraints['nodes']):
                    c1, c2 = c[0], c[1]
                    remove = False
                    if c1 == n1 and c2 == n2:
                        remove = True
                    if c1 == n2 and c2 == n1:
                        remove = True
                    if remove:
                        indexes.append(i)
                for r in indexes:
                    self.__constraints['nodes'].pop(r)
            
            def remove_node_to_link(n, n1, n2):
                indexes = []
                for i, c in enumerate(self.__constraints['links']):
                    c, c1, c2 = c[0], c[1], c[2]
                    remove = False
                    if c == n and c1 == n1 and c2 == n2:
                        remove = True
                    if remove:
                        indexes.append(i)
                for r in indexes:
                    self.__constraints['links'].pop(r)
            
            for n1 in remove:
                for n2 in line1:
                    remove_node_to_node(n1, n2)
            for n1 in remove:
                for n2 in line2:
                    remove_node_to_node(n1, n2)
            for n in line1:
                for n1 in remove:
                    for n2 in remove:
                        remove_node_to_link(n, n1, n2)
            for n in line2:
                for n1 in remove:
                    for n2 in remove:
                        remove_node_to_link(n, n1, n2)


        # Get the surrounding points in a shape with a given threshold
        def get_surrounding_in_shape(shape, value, threshold, previous_only=False, next_only=False):
            idx = shape.index(value)
            result = shape[:]
            result.remove(value)

            if next_only:
                left = idx
            else:
                left = idx - threshold
                left = min(max(0, left), len(result) - (2 * threshold))
            if previous_only:
                right = idx
            else:
                if next_only:
                    t = threshold
                else:
                    t = 2 * threshold
                right = left + t

            return result[left:right]

        # Check if two LineStrings intersect
        def flag_crossing(n, n1, n2, shape, shape1):
            # If the lines cross at one node
            if n == n1:
                add = True
                for nc in node_crossing:
                    if nc[0] == n:
                        add = False
                if add:
                    node_crossing.append([n, shape, shape1])
            else:
                # Create the line with the two nodes
                line = shapely.LineString([points[n1], points[n2]])
                # Retrieve the index of the node in the shape
                nid = shape.index(n)

                crossing = None
                # If it's not the first node of the shape
                if nid > 0:
                    # Create the line between the node and the previous node
                    pp = shape[nid - 1]
                    linep = shapely.LineString([points[n], points[pp]])
                    # Check if both lines crosses
                    if shapely.crosses(line, linep):
                        crossing = [n, pp, n1, n2]
                # If it's not the last node of the shape
                if nid < (len(shape) - 1) and crossing is None:
                    # Create the line between the node and the next node
                    pn = shape[nid + 1]
                    linen = shapely.LineString([points[n], points[pn]])
                    # Check if both lines crosses
                    if shapely.crosses(line, linen):
                        crossing = [n, pn, n1, n2]

                # If a crossing has been found
                if crossing is not None:
                    add = True
                    # Check if it's not already present before adding to the list
                    for ec in edge_crossing:
                        if Counter(ec[0]) == Counter(crossing):
                            add = False
                    if add:
                        edge_crossing.append([crossing, shape, shape1])


        # Check if a node to link conflict already exists before adding one
        def add_node_to_link(c):
            add = True
            for n in self.__constraints['links']:
                if (c[0] == n[0] and c[1] == n[1] and c[2] == n[2]):
                    add = False
                if (c[0] == n[0] and c[1] == n[2] and c[2] == n[1]):
                    add = False
            if add:
                self.__constraints['links'].append(c)

        # Check if a node to node conflict already exists before adding one
        def add_node_to_node(c):
            add = False
            if c[0] != c[1]:
                add = True
                for n in self.__constraints['nodes']:
                    if (c[0] == n[0] and c[1] == n[1]):
                        add = False
                    if (c[0] == n[1] and c[1] == n[0]):
                        add = False
            if add:
                self.__constraints['nodes'].append(c)

        def retrieve_nodes_links(shape, shape1, geomtype, geomtype1, p, conflict_dist, min_dist, weight):
            point = shapely.Point(points[p])
            # Loop through all the points of the LineString or Polygon shape
            for pid1, p1 in enumerate(shape1):
                notlast = True
                # Check if it's not the last point of the shape
                if pid1 < (len(shape1) - 1):
                    p2 = shape1[pid1 + 1]
                else:
                    notlast = False
                    if geomtype1 == 'Polygon':
                        p2 = shape1[0]
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
                # First, checks if one of the three distance is below the minimum distance threshold
                if pdist1 < conflict_dist or pdist2 < conflict_dist or ldist < conflict_dist:
                    # For LineString, check if the two lines intersects
                    if geomtype == geomtype1 == 'LineString':
                        flag_crossing(p, p1, p2, shape, shape1)

                    # If the line distance is the smallest, insert a node to link spatial conflict
                    if ldist < pdist1 and ldist < pdist2:
                        # if line_crossing(p, point, line, points) == False:
                        add_node_to_link([p, p1, p2, min_dist, ldist, weight])
                    else:                      
                        # Determine which node is the closest
                        if pdist1 < pdist2:
                            add_node_to_node([p, p1, min_dist, pdist1, weight])
                        else:
                            if notlast:
                                add_node_to_node([p, p2, min_dist, pdist2, weight])

        # Check if the considered point is spatially conflicting with other points
        def check_conflicts(shape, p, oid, sid):
            geomtype = self.__OBJECTS[oid].geometry[sid].geom_type
            point = shapely.Point(points[p])
            # Loop through all objects
            for oid1, o in enumerate(self.__shapes):
                # Loop through all shapes
                for sid1, s in enumerate(o):
                    # Retrieve the geometry of the shape
                    shape1 = self.__OBJECTS[oid1].geometry[sid1]
                    # Stores the geometry type
                    geomtype1 = shape1.geom_type

                    # Checks if it's the same object
                    if oid1 == oid:
                        if geomtype1 == 'LineString':
                            continue
                        if sid1 == sid:
                            continue
                        if self.INTRA_CONFLICTS == False:
                            continue
                    # Getting the weight of the spatial conflict constraint from the matrix
                    weight = self.spatial_weights[oid][oid1]
                    # Getting the distance value from the distances matrix
                    min_dist = self.distances[oid][oid1]
                    # Setting a distance equal to 1.5 times the min distance to retrieve conflicting objects
                    # Setting conflict_dist = min_dist doesn't take the 1.5 time conflicting entities
                    conflict_dist = 1.5 * min_dist
                    
                    # Checks if that shape is within the minimum distance before retrieving pairs of nodes and links
                    if shapely.dwithin(point, shape1, conflict_dist):
                        # Checks if the shape is a point
                        if geomtype1 == 'Point':
                            distance = shapely.distance(point, shape1)
                            add_node_to_node([p, s[0], min_dist, distance, weight])
                        # If it's not, checking if it's closer to a node or a segment
                        else:
                            retrieve_nodes_links(shape, s, geomtype, geomtype1, p, conflict_dist, min_dist, weight)

        # Loop through all objects
        for oid, o in enumerate(self.__shapes):
            # Loop through all shapes
            for sid, s in enumerate(o):
                # Loop through all points
                for p in s:
                    # Check conflicts with other points
                    check_conflicts(s, p, oid, sid)

        threshold = self.CROSSING_NODE_THRESHOLD

        # Once every spatial conflicts have been created, removing conflicts around crossing lines
        for nc in node_crossing:
            node, line1, line2 = nc[0], nc[1], nc[2]
            remove = [node]
            remove.extend(get_surrounding_in_shape(line1, node, threshold))
            remove.extend(get_surrounding_in_shape(line2, node, threshold))
            remove_conflicts(remove, line1, line2)

        if threshold > 0:
            threshold -= 1
            for ec in edge_crossing:
                p1, p2, p3, p4 = ec[0][0], ec[0][1], ec[0][2], ec[0][3]
                line1, line2 = ec[1], ec[2]
                remove = [p1, p2, p3, p4]
                remove.extend(get_surrounding_in_shape(line1, p1, threshold, previous_only=True))
                remove.extend(get_surrounding_in_shape(line1, p2, threshold, next_only=True))
                remove.extend(get_surrounding_in_shape(line2, p3, threshold, previous_only=True))
                remove.extend(get_surrounding_in_shape(line2, p4, threshold, next_only=True))
                remove_conflicts(remove, line1, line2)

    def __update_distances(self, points):
        """
        Update the distance between pairs of nodes and pairs of nodes and links
        """
        # Loop through node to node conflicts
        for i, n in enumerate(self.__constraints['nodes']):
            # Retrieve coordinates of the two nodes
            n1, n2 = points[n[0]], points[n[1]]
            x1, x2, y1, y2 = n1[0], n2[0], n1[1], n2[1]
            # Calculate the vector formed by those two points
            v = np.array([x2 - x1, y2 - y1])
            # Calculate its norm, i.e. distance and updating it
            nodedistance = np.linalg.norm(v)
            self.__constraints['nodes'][i][3] = nodedistance

        # Loop through node to link conflicts
        for i, l in enumerate(self.__constraints['links']):
            # Retrieve coordinates of the three nodes as array
            n0 = np.array(points[l[0]])
            n1 = np.array(points[l[1]])
            n2 = np.array(points[l[2]])
            # Calculate the distance between the first node and the line formed by the other two
            linkdistance = np.linalg.norm(np.cross(n2 - n1, n1 - n0)) / np.linalg.norm(n2 - n1)
            # Updating the distance
            self.__constraints['links'][i][4] = linkdistance

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

    def get_vertexes(self):
        """
        Returns the full list of vertexes.

        Returns
        -------
        GeoDataFrame of Point
        """
        pointslist = []
        shape_id = 0
        for o in self.__shapes:
            for s in o:
                for p in s:
                    pointslist.append({
                        'id' : p,
                        'object' : shape_id,
                        'geometry' : shapely.Point(self.__points[p])
                    }) 
                shape_id += 1         
        return gpd.GeoDataFrame(pointslist)

    def get_vertexes_conflicts(self):
        """
        Returns all the vertex to vertex conflicts as a
        LineString between the two nodes.
        
        Returns
        -------
        GeoDataFrame of LineString
        """
        lines = []
        for i, node in enumerate(self.__constraints['nodes']):
            lines.append({
                'nid0': node[0],
                'nid1': node[1],
                'group': i,
                'geometry': shapely.LineString([self.__points[node[0]], self.__points[node[1]]])
            })

        if len(lines) > 0:
            return gpd.GeoDataFrame(lines, crs=3857)
        else:
            return None

    def get_links_conflicts(self):
        """
        Returns all the vertex to links conflicts as a LineString between
        the vertex and the projection of this vertex on the conflicting segment.

        Returns
        -------
        GeoDataFrame of LineString
        """
        lines = []
        for lid, link in enumerate(self.__constraints['links']):
            p0 = shapely.Point(self.__points[link[0]])
            p1 = shapely.Point(self.__points[link[1]])
            p2 = shapely.Point(self.__points[link[2]])
            line = shapely.LineString([p1, p2])
            dist = line.project(p0)
            projected = line.interpolate(dist)
            lines.append({
                'lid': lid,
                'geometry': shapely.LineString([p0, projected])
            })

        if len(lines) > 0:
            return gpd.GeoDataFrame(lines, crs=3857)
        else:
            return None
        
    def get_objects_number(self):
        """
        Returns the number of objects
        added to the generalisation algorithm.

        Returns
        -------
        int
        """
        return len(self.__OBJECTS)