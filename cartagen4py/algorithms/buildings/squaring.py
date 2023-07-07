# This is an implementation of the least squares based squaring algorithm proposed by Lokhat & Touya (https://hal.science/hal-02147792)
import numpy as np
from shapely.geometry import Polygon

# this function encapsulates the above class methods to square a collection of polygons and return the squared geometries.
def square_polygon(polygons,max_iteration=1000, norm_tolerance=0.05,
            right_tolerance=10, flat_tolerance=10,
            fixed_weight=5, right_weight=100, flat_weight=50):
    squarer = Squarer(max_iteration, norm_tolerance,
            right_tolerance, flat_tolerance, 7,
            fixed_weight, right_weight, flat_weight, 10)
    new_points = squarer.square(polygons)
    list_vertices = squarer.get_shapes_from_new_points(polygons, new_points)
    new_geoms = []
    for vertices in list_vertices:
        new_geoms.append(Polygon(vertices))
    return new_geoms

class Squarer:
    """Initialize squaring object, with default weights and tolerance set in the constructor
    """
    def __init__(
            self, max_iteration=1000, norm_tolerance=0.05,
            right_tolerance=10, flat_tolerance=10, half_right_tolerance=7,
            fixed_weight=5, right_weight=100, flat_weight=50, half_right_weight=10, switch_new=False
        ):

        self.SWITCH_NEW = switch_new
        self.MAX_ITERATION = max_iteration
        self.NORM_TOLERANCE = norm_tolerance
        self.right_tolerance = right_tolerance # 10 90° angles tolerance
        self.flat_tolerance = flat_tolerance # 10 flat angles tolerance
        self.half_right_tolerance = half_right_tolerance # 7 45/135° angles tolerance
        self.fixed_weight = fixed_weight #5
        self.right_weight = right_weight #100
        self.flat_weight = flat_weight #50
        self.half_right_weight = half_right_weight #10

        self.point_shapes = {}
        self.__lines_pindex = []
        self.indicesRight, self.indicesFlat, self.indicesHrAig, self.indicesHrObt = [], [], [], []

    # returns a list of coordinates from a shapely linestring, multilinestring or polygon 
    def __get_coords(self, shape):
        if shape.geom_type == 'MultiLineString':
            return shape[0].coords 
        elif shape.geom_type == 'Polygon':
            return shape.exterior.coords
        elif shape.geom_type == 'LineString':
            return shape.coords
        return []

    # rank of key p in a dict
    def __get_rank(self, dict_p, p):
        for i, e in enumerate(dict_p):
            if e == p:
                return i
        return -1

    # get a dict where the key is the tuple of a point and the associated value is 
    # the table of indices of the geometries the point belongs to
    def __build_dict_of_unique_points(self, shapes):
        for i, s in enumerate(shapes):
            #coords = s[0].coords if shapes[0].geom_type == 'MultiLineString' else s.exterior.coords
            coords = self.__get_coords(s)
            for p in coords:
                #print(p) # tuple
                if p in self.point_shapes:
                    if i not in self.point_shapes[p]:
                        self.point_shapes[p].append(i)
                else:
                    self.point_shapes[p] = [i]


    # get a list where each geometry is made of a list of the indices of the points forming the geometry.
    def __build_pindex_for_shapes(self, shapes):
        for s in shapes:
            index_points = []
            #coords = s[0].coords if shapes[0].geom_type == 'MultiLineString' else s.exterior.coords
            coords = self.__get_coords(s)
            for p in coords:
                index_points.append(self.__get_rank(self.point_shapes, p))
            self.__lines_pindex.append(index_points)

    # get the rank of the point of index idx_p in the shape of index idx_s
    def __get_rank_point_in_shape(self, idx_p, idx_s):
        for i, idx_pp in enumerate(self.__lines_pindex[idx_s]):
            if idx_pp == idx_p: return i
        return -1

    # get the list of potential angles around the point of index idx_p as a triplet of indices [[idx_prev, idx_p, idx_next]...]
    def __get_angle_triplets(self, idx_p, unik_points):
        p = unik_points[idx_p]
        lines_containing_p = self.point_shapes[p]
        if len(lines_containing_p) == 1:
            idx_l = lines_containing_p[0]
            r = self.__get_rank_point_in_shape(idx_p, idx_l)
            # single node, no angle
            if r == 0 or r == len(self.__lines_pindex[idx_l]) - 1 : 
                return []
            # interior point of a line
            idx_prev = self.__lines_pindex[idx_l][r - 1]
            idx_next = self.__lines_pindex[idx_l][r + 1]
            #print(f'POINT({p[0]} {p[1]})')
            return [[idx_prev, idx_p, idx_next]]
        # points intersecting multiple lines
        else:
            triplets = []
            #print(f'POINT({p[0]} {p[1]})')
            for i in range(len(lines_containing_p) - 1):
                idx_l1 = lines_containing_p[i]
                r = self.__get_rank_point_in_shape(idx_p, idx_l1)
                idx_prev = self.__lines_pindex[idx_l1][r - 1] if r > 0 else self.__lines_pindex[idx_l1][r + 1]
                for j in range(i + 1, len(lines_containing_p)):
                    idx_l2 = lines_containing_p[j]
                    r = self.__get_rank_point_in_shape(idx_p, idx_l2)
                    idx_next = self.__lines_pindex[idx_l2][r - 1] if r > 0 else self.__lines_pindex[idx_l2][r + 1]
                    triplets.append([idx_prev, idx_p, idx_next])
            return triplets

    def __get_vecs_around(self, t, unik_points) : # t = [idx_prec, idx_p, idx_suiv]
        """ return vectors formed by a triplet of indexes
        """
        pr, p, s = unik_points[t[0]], unik_points[t[1]], unik_points[t[2]]
        if self.SWITCH_NEW:
            v1 = np.array([pr[0] - p[0], pr[1] - p[1]])
        else:
            v1 = np.array([p[0] - pr[0], p[1] - pr[1]])
        v2 = np.array([s[0] - p[0], s[1] - p[1]])
        return v1, v2

    def __idx_angles_remarquables(self, unik_points):
        #unik_points = list(self.point_shapes)
        rTol = np.cos((np.pi / 2) - self.right_tolerance * np.pi / 180)
        fTol = np.deg2rad(self.flat_tolerance)
        hrTol1 = np.cos((np.pi / 4) - self.half_right_tolerance * np.pi / 180)
        hrTol2 = np.cos((np.pi / 4) + self.half_right_tolerance * np.pi / 180)
        for idx_p in range(len(unik_points)):
            triplets = self.__get_angle_triplets(idx_p, unik_points)
            for t in triplets:
                v1, v2 = self.__get_vecs_around(t, unik_points)
                n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
                v1n = v1 / n1 if n1 != 0. else np.array([0.,0.]) #n1
                v2n = v2 / n2 if n2 != 0. else np.array([0.,0.]) #n2
                dot = v1n.dot(v2n)
                cross = np.cross(v1n, v2n).item(0)
                if (np.abs(dot) <= rTol):
                    self.indicesRight.append(t)
                elif (cross <= fTol):
                    self.indicesFlat.append(t)
                #elif (dot <= hrTol1 and dot >= hrTol2):
                #    self.indicesHrAig.append(t)
                #elif (dot >= -hrTol1 and dot <= -hrTol2):
                #    self.indicesHrObt.append(t)
        print(f'potential angles -- R: {len(self.indicesRight)} - F: {len(self.indicesFlat)}')
        #print(f'potential angles -- R: {len(self.indicesRight)} - F: {len(self.indicesFlat)} - HRa: {len(self.indicesHrAig)} - HRo: {len(self.indicesHrObt)}')


    def __get_Y(self, unik_points):
        """ Observation vector
        """
        nb_points = len(unik_points)
        self.Y = np.zeros(2 * nb_points + len(self.indicesRight) + len(self.indicesFlat))
        #self.Y = np.zeros(2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrObt) + len(self.indicesHrAig))
        for i, p in enumerate(unik_points):
            self.Y[2*i] = p[0]
            self.Y[2*i+1] = p[1]
        #offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat)
        #for i, t in enumerate(self.indicesHrAig):
        #    v1, v2 = self.__get_vecs_around(t, unik_points)
        #    d = np.linalg.norm(v1) * np.linalg.norm(v2) * np.cos(np.pi / 4)
        #    self.Y[offset + i] = d
        #offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrAig)
        #for i, t in enumerate(self.indicesHrObt):
        #    v1, v2 = self.__get_vecs_around(t, unik_points)
        #    d = np.linalg.norm(v1) * np.linalg.norm(v2) * np.cos(3 * np.pi / 4)
        #    self.Y[offset + i] = d

    # B = Y - S(Xcourant)
    def __get_B(self, points):
        nb_points = len(points)
        S = np.zeros(2 * nb_points + len(self.indicesRight) + len(self.indicesFlat))
        #S = np.zeros(2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrObt) + len(self.indicesHrAig))
        for i, p in enumerate(points):
            S[2*i] = p[0]
            S[2*i+1] = p[1]
        offset = 2 * nb_points
        for i, t in enumerate(self.indicesRight):
            v1, v2 = self.__get_vecs_around(t, points)
            d = v1.dot(v2)
            S[offset + i] = d
        offset = 2 * nb_points + len(self.indicesRight)
        for i, t in enumerate(self.indicesFlat):
            v1, v2 = self.__get_vecs_around(t, points)
            d = np.cross(v1, v2).item(0)
            S[offset + i] = d
        offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat)
        #for i, t in enumerate(self.indicesHrAig):
        #    v1, v2 = self.__get_vecs_around(t, points)
        #    d = v1.dot(v2) 
        #    S[offset + i] = d
        #offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrAig)
        #for i, t in enumerate(self.indicesHrObt):
        #    v1, v2 = self.__get_vecs_around(t, points)
        #    d = v1.dot(v2) 
        #    S[offset + i] = d
        return self.Y - S

    # Weight Matrix
    # n = 2 * nb_points + indicesRight.size() + indicesFlat.size() + indicesHrAig.size() + indicesHrObt.size()
    def __get_P(self):
        nb_points = len(self.point_shapes)
        nb_rights, nb_flats =  len(self.indicesRight), len(self.indicesFlat)
        #nb_half_rights = len(self.indicesHrAig) + len(self.indicesHrObt)
        wfix = np.full(2*nb_points, self.fixed_weight)
        wRight = np.full(nb_rights, self.right_weight)
        wFlat = np.full(nb_flats, self.flat_weight)
        #wHr = np.full(nb_half_rights, self.half_right_weight)
        self.P = np.diag(np.concatenate((wfix, wRight, wFlat)))

    ## new vectors
    def __partial_derivatives_dotp(self, points, indices):
        nb_points = len(points)
        nb_indices = len(indices)
        m = np.zeros((nb_indices, 2*nb_points))
        for i, t in enumerate(indices):
            idx_prec, idx, idx_suiv = t[0], t[1], t[2]
            pr, p, s = points[t[0]], points[t[1]], points[t[2]]
            # df en Xi-1, Yi-1
            dfx = p[0] - s[0]
            dfy = p[1] - s[1]
            m[i][2*idx_prec] = dfx
            m[i][2*idx_prec + 1] = dfy
            # df en Xi, Yi
            dfx = s[0] - 2*p[0] + pr[0]
            dfy = s[1] - 2*p[1] + pr[1]
            m[i][2*idx] = dfx
            m[i][2*idx + 1] = dfy
            # df en Xi+1, Yi+1
            dfx = p[0] - pr[0]
            dfy = p[1] - pr[1]
            m[i][2*idx_suiv] = dfx
            m[i][2*idx_suiv + 1] = dfy
        return m

    def __partial_derivatives_cross(self, points, indices):
        nb_points = len(points) #- 1
        nb_indices = len(indices)
        m = np.zeros((nb_indices, 2*nb_points))
        for i, t in enumerate(indices):
            idx_prec, idx, idx_suiv = t[0], t[1], t[2]
            pr, p, s = points[t[0]], points[t[1]], points[t[2]]
            # df en Xi-1, Yi-1
            dfx = p[1] - s[1]
            dfy = -p[0] + s[0]
            m[i][2*idx_prec] = dfx
            m[i][2*idx_prec + 1] = dfy
            # df en Xi, Yi
            dfx = s[1] - pr[1]
            dfy = -s[0] + pr[0]
            m[i][2*idx] = dfx
            m[i][2*idx + 1] = dfy
            # df en Xi+1, Yi+1
            dfx = -p[1] + pr[1]
            dfy = p[0] - pr[0]
            m[i][2*idx_suiv] = dfx
            m[i][2*idx_suiv + 1] = dfy
        return m


    def __get_A(self, points):
        nb_points = len(points) #- 1
        id = np.identity(2 * nb_points)
        partialR = self.__partial_derivatives_dotp(points, self.indicesRight)
        partialCross = self.__partial_derivatives_cross(points, self.indicesFlat)
        #partialHr1 = self.__partial_derivatives_dotp(points, self.indicesHrAig)
        #partialHr2 = self.__partial_derivatives_dotp(points, self.indicesHrObt)
        a = np.vstack((id, partialR))
        if isinstance(partialCross, np.ndarray):
            a = np.vstack((a, partialCross))
        # if isinstance(partialHr1, np.ndarray):
        #    a = np.vstack((a, partialHr1))
        #if isinstance(partialHr1, np.ndarray):
        #    a = np.vstack((a, partialHr2))
        return a

    def __compute_dx(self, points):
        A = self.__get_A(points)
        B = self.__get_B(points)
        atp = A.T @ self.P 
        atpa = atp @ A
        atpb = atp @ B
        #dx = np.linalg.lstsq(atpa, atpb)
        #dx = np.linalg.inv(atpa) @ atpb
        dx = np.linalg.solve(atpa, atpb)
        return dx
    
    def __prepare_square(self, shapes):
        if len(shapes) == 0:
            return np.array([])
        self.geom_type = shapes[0].geom_type
        self.__build_dict_of_unique_points(shapes)
        self.__build_pindex_for_shapes(shapes)
        unik_points = list(self.point_shapes)
        self.__idx_angles_remarquables(unik_points)
        self.__get_Y(unik_points)
        self.__get_P()
        return np.array(unik_points)
    
    def square(self, shapes):
        """squares a collection of shapely multilinestrings or polygons
        returns a numpy array of the points after the least square process
        """
        points = self.__prepare_square(shapes)
        nb_points = len(points)
        for i in range(self.MAX_ITERATION):
            dx = self.__compute_dx(points)
            points += dx.reshape((nb_points, 2))
            print(i, np.linalg.norm(dx, ord=np.inf))
            if np.linalg.norm(dx, ord=np.inf) < self.NORM_TOLERANCE:
                break
        self.nb_iters = i
        return points


    # rebuild shapes with updated points from least square process
    def get_shapes_from_new_points(self, original_shapes, new_points):
        """rebuild a list of coordinates from the original collection  of shapes
        and the points obtained from the square process
        """
        unik_points = list(self.point_shapes)
        new_s = []
        for l in original_shapes:
            #coords, is_poly = (l[0].coords, False) if l.geom_type == 'MultiLineString' else (l.exterior.coords, True)
            coords = self.__get_coords(l)
            is_poly = True if self.geom_type == 'Polygon' else False
            size = len(coords)
            new_s.append(np.zeros((size, 2)) )
        for idx_p, p in enumerate(unik_points):
            index_of_lines = self.point_shapes[p] # point_lines_idx[p]
            #print(index_of_lines)
            for idx_l in index_of_lines:
                r = self.__get_rank_point_in_shape(idx_p, idx_l)
                #print(idx_l, r, new_points[idx_p])
                new_s[idx_l][r] = np.array(new_points[idx_p])
                if r == 0 and is_poly:
                    new_s[idx_l][-1] = np.array(new_points[idx_p])
        return new_s
