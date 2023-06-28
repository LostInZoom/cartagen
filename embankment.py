from timeit import default_timer as timer
import logging
import numpy as np
import pygeos

#log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
log_format = '%(message)s'
formatter = logging.Formatter(log_format)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

loglsd = logging.getLogger(__name__)
loglsd.addHandler(stream_handler)
loglsd.setLevel(logging.INFO)

class LSDisplacer:
    MAX_ITER = 250
    NORM_DX = 0.3 
    DIST = 'MIN'
    KKT = False
    # constraints
    ID_CONST = True
    ANGLES_CONST = True
    EDGES_CONST = True
    DIST_CONST = True
    # weights and differentiation step
    H = 2.0 
    PFix = 8.
    PEdges_int = 1
    PEdges_int_non_seg = 1
    PEdges_ext = 15
    Pedges_ext_far = 0.5
    PAngles = 8
    PDistRoads = 200
    FLOATING_NORM = True

    def __init__(self, points_talus, roads_wkts_and_buffers, talus_lengths, edges, edges_dist_min=10, edges_dist_max=30):
        self.points_talus = points_talus
        self.roads_shapes =  pygeos.from_wkt([r[0] for r in roads_wkts_and_buffers])
        self.buffers = np.array([r[1] for r in roads_wkts_and_buffers])
        self.talus_lengths = talus_lengths
        self.edges = edges
        self.edges_dist_min = edges_dist_min
        self.edges_dist_max = edges_dist_max
        self.x_courant = self.points_talus.copy()
        self.nb_vars = len(self.x_courant)
        self.e_lengths_ori = self._get_edges_lengths(self.edges)
        self.angles_ori = self.angles_crossprod()
        self.P = self.get_P()
        self.meta = {"nb_iters": -1, "time_s": -1, "dx_reached": -1}


    def set_params(MAX_ITER=250, NORM_DX=0.3, H=2.0, DIST='MIN', KKT=False,
                   ID_CONST=True, ANGLES_CONST=True, EDGES_CONST=True, DIST_CONST=True,
                   PFix=1., PEdges_int=10, PEdges_ext=2, PAngles=50, Pedges_ext_far=0, PEdges_int_non_seg=5, PDistRoads=1000):
        LSDisplacer.MAX_ITER = MAX_ITER
        LSDisplacer.NORM_DX = NORM_DX
        LSDisplacer.DIST = DIST
        LSDisplacer.KKT = KKT
        LSDisplacer.ID_CONST = ID_CONST
        LSDisplacer.ANGLES_CONST = ANGLES_CONST
        LSDisplacer.EDGES_CONST = EDGES_CONST
        LSDisplacer.DIST_CONST = DIST_CONST
        LSDisplacer.H = H
        LSDisplacer.PFix = PFix
        LSDisplacer.PEdges_int = PEdges_int
        LSDisplacer.PEdges_int_non_seg = PEdges_int_non_seg
        LSDisplacer.PEdges_ext = PEdges_ext
        LSDisplacer.Pedges_ext_far = Pedges_ext_far
        LSDisplacer.PAngles = PAngles
        LSDisplacer.PDistRoads = PDistRoads


    # length of edge referenced by index idx (pt_ini, pt_fin) in pts
    def _edge_length(self, idx ): #, pts):
        xa, ya, xb, yb = self.x_courant[idx[0]*2], self.x_courant[idx[0]*2 + 1], self.x_courant[idx[1]*2], self.x_courant[idx[1]*2 + 1]
        return ((xa - xb)**2 + (ya - yb)**2)**0.5

    # talus' number for point at index idx 
    def _num_talus(self, idx):
        i = 0
        s = self.talus_lengths[0]
        while idx >= s:
            i += 1
            s += self.talus_lengths[i]
        return i

    # pygeos, returns an array of linestrings from the points of all talus lines
    def _lines_from_points(points, talus_lengths):
        offset = 0
        lines = []
        coords = points.reshape(-1, 2)
        for size in talus_lengths:
            l = pygeos.creation.linestrings(coords[offset: offset+size])
            lines.append(l)
            offset += size
        return np.array(lines)

    # pygeos specific, return a multiline from all points of all talus
    def _multiline_from_points(points, talus_lengths):
        lines = LSDisplacer._lines_from_points(points, talus_lengths) 
        return pygeos.creation.multilinestrings(lines)
    

    def print_linestrings_wkts(self):
        lines = self.get_linestrings_wkts()
        for l in lines:
            print(l)
    
    def get_linestrings_wkts(self):
        points = self.x_courant.reshape(-1, 2)
        lines = LSDisplacer._lines_from_points(points, self.talus_lengths)
        lines = pygeos.to_wkt(lines, rounding_precision=-1)
        return lines


    # set edges original lengths, and set minimal distance for too small inter talus edges
    def _get_edges_lengths(self, edges):
        edges_lengths_ori = []
        for e in edges:
            el = self._edge_length(e) #, self.points_talus)
            if el < self.edges_dist_min and self._num_talus(e[0]) != self._num_talus(e[1]):
                loglsd.debug("min reached inter edges")
                el = self.edges_dist_min
            edges_lengths_ori.append(el)
        return np.array(edges_lengths_ori)

    def angles_crossprod(self):
        """ returns normalized crossproduct of all angles for "inside" points in all talus lines
        """
        offset = 0
        cross_products = []
        for size in self.talus_lengths:
            for i in range(offset + 1, (offset + size - 1)):
                prec = np.array((self.x_courant[2*i - 2], self.x_courant[2*i - 1]))
                pt = np.array((self.x_courant[2*i], self.x_courant[2*i + 1]))
                suiv = np.array((self.x_courant[2*i + 2], self.x_courant[2*i + 3]))
                u = (pt - prec) #/ np.linalg.norm(pt - prec)
                u = u / np.linalg.norm(u)
                v = (suiv - pt) #/ np.linalg.norm(suiv - pt)
                v = v / np.linalg.norm(v)
                cross_products.append(np.cross(u, v))
            offset += size
        return np.array(cross_products)
   
    # derived with sympy, micro optimized
    # raw from sympy it is :
    # # df/dxi-1
    # m[2*i - 2] = (((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))**(-0.5)*(-1.0*(x - xpp)*((x - xpp)*(y - yss) - (x - xss)*(y - ypp)) + (y - yss)*((x - xpp)**2 + (y - ypp)**2))/((x - xpp)**2 + (y - ypp)**2)
    # # df/dyi-1
    # m[2*i - 1] = (((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))**(-0.5)*((-x + xss)*((x - xpp)**2 + (y - ypp)**2) - 1.0*(y - ypp)*((x - xpp)*(y - yss) - (x - xss)*(y - ypp)))/((x - xpp)**2 + (y - ypp)**2)
    # # df/dxi
    # m[2*i] = (((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))**(-0.5)*((-ypp + yss)*((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2) + ((x - xpp)*(y - yss) - (x - xss)*(y - ypp))*((x - xpp)*((x - xss)**2 + (y - yss)**2) + (x - xss)*((x - xpp)**2 + (y - ypp)**2)))/(((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))
    # # df/dyi
    # m[2*i + 1] = (((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))**(-0.5)*((xpp - xss)*((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2) + ((x - xpp)*(y - yss) - (x - xss)*(y - ypp))*((y - ypp)*((x - xss)**2 + (y - yss)**2) + (y - yss)*((x - xpp)**2 + (y - ypp)**2)))/(((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))
    # # df/dxi+1
    # m[2*i + 2] = (((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))**(-0.5)*(-1.0*(x - xss)*((x - xpp)*(y - yss) - (x - xss)*(y - ypp)) + (-y + ypp)*((x - xss)**2 + (y - yss)**2))/((x - xss)**2 + (y - yss)**2)
    # # df/dyi+1
    # m[2*i + 3] = (((x - xpp)**2 + (y - ypp)**2)*((x - xss)**2 + (y - yss)**2))**(-0.5)*((x - xpp)*((x - xss)**2 + (y - yss)**2) - 1.0*(y - yss)*((x - xpp)*(y - yss) - (x - xss)*(y - ypp)))/((x - xss)**2 + (y - yss)**2)                
    def cross_norm_diff(self):
        offset = 0
        cross_products = []
        for size in self.talus_lengths:
            m = np.zeros(self.nb_vars)
            for i in range(offset + 1, offset + size - 1):
                # xi, yi => x, y | xi-1, yi-1 => xpp, ypp | xi+1, yi+1 => xss, yss
                u = self.x_courant[2 * i - 2:2 * i + 4]
                xpp, ypp, x, y, xss, yss = u[0], u[1], u[2], u[3], u[4], u[5]
                x_xpp = x - xpp
                y_ypp = y - ypp
                x_xss = x - xss
                y_yss = y - yss
                b = ((x_xpp)**2 + (y_ypp)**2)
                d = ((x_xss)**2 + (y_yss)**2)
                c = ((x_xpp)*(y_yss) - (x_xss)*(y_ypp))
                bd = b*d
                a = (bd)**(-0.5)
                # df/dxi-1
                m[2*i - 2] = a*(-x_xpp*c + y_yss*b)/b
                # df/dyi-1
                m[2*i - 1] = a*(-x_xss*b - y_ypp*c)/b
                # df/dxi
                m[2*i] = a*((-ypp + yss)*bd + c*((x_xpp)*d + (x_xss)*b))/bd
                # df/dyi
                m[2*i + 1] = a*((xpp - xss)*bd + c*((y_ypp)*d + (y_yss)*b))/bd
                # df/dxi+1
                m[2*i + 2] = a*(-x_xss*c - y_ypp*d)/d
                # df/dyi+1
                m[2*i + 3] = a*(x_xpp*d - y_yss*c)/d
                cross_products.append(m)
            offset += size
        return np.array(cross_products)
      
    # each line of points_array contains points for multiple lines, offset and size is deduced from self.talus_lengths
    # returns an array of distances from road, either min or mean
    def dist_F_vectorized(self, road, i, points_array):
        ml = []
        for c in points_array:
            if LSDisplacer.DIST == 'MIN':
                m = LSDisplacer._multiline_from_points(c, self.talus_lengths)
            else:
                m = LSDisplacer._lines_from_points(c, self.talus_lengths) 
            ml.append(m)
        ml = np.array(ml)
        dists = pygeos.distance(road, ml)
        if LSDisplacer.DIST != 'MIN':
            dists = dists.mean(axis=1)
        dists = np.where(dists > self.buffers[i], 0., self.buffers[i] - dists)
        return dists


    def dist_F_diff(self, road, i):
        # diagonal matrix with H on diagonal
        h = np.eye(self.nb_vars) * self.H
        coords_plus_H = self.x_courant + h
        coords_minus_H = self.x_courant - h
        # seems a bit faster to have 2 np arrays instead of the same one splitted 
        d_plus = self.dist_F_vectorized(road, i, coords_plus_H)
        d_min = self.dist_F_vectorized(road, i, coords_minus_H)
        ds = (d_plus - d_min) / (2* self.H)
        return ds

    # idx : edge index in pts [idx_p1, idx_p2]
    def edge_length_diff(self, idx):
        xa, ya, xb, yb = self.x_courant[idx[0]*2], self.x_courant[idx[0]*2 + 1], self.x_courant[idx[1]*2], self.x_courant[idx[1]*2 + 1]
        ab = ((xa - xb)**2 + (ya - yb)**2) ** 0.5
        l = np.zeros(self.nb_vars)
        l[idx[0] * 2] = (xa - xb) / ab
        l[idx[0] * 2 + 1] = (ya - yb) / ab
        l[idx[1] * 2] = -l[idx[0] * 2] # -(xa - xb) / ab
        l[idx[1] * 2 + 1] = -l[idx[0] * 2 + 1] #- (ya - yb) / ab
        return l

    def get_P(self):
        weights = []
        if LSDisplacer.ID_CONST:
            wfix = np.full(self.nb_vars, LSDisplacer.PFix)
            weights.append(wfix)
        if LSDisplacer.ANGLES_CONST:
            wAngles = np.full(len(self.angles_ori), LSDisplacer.PAngles)
            weights.append(wAngles)
        if LSDisplacer.EDGES_CONST:
            wEdges = []
            for i, e in enumerate(self.edges):
                same_talus = self._num_talus(e[0]) == self._num_talus(e[1])
                non_consecutive_points = abs(e[0] - e[1]) != 1
                if same_talus:
                    if non_consecutive_points:
                        loglsd.debug("**** non intra edge segment : limiting weight")
                        wEdges.append(LSDisplacer.PEdges_int_non_seg)
                    else:
                        loglsd.debug("**** intra edge segment")
                        wEdges.append(LSDisplacer.PEdges_int)
                else:
                    if self._edge_length(e) >= self.edges_dist_max: 
                        loglsd.debug("**** max inter edges threshold reached : minimalizing weight for this edge")
                        wEdges.append(LSDisplacer.Pedges_ext_far)
                    else:
                        loglsd.debug("**** inter edges segment")
                        wEdges.append(LSDisplacer.PEdges_ext)
            wEdges = np.array(wEdges)
            weights.append(wEdges)
        if LSDisplacer.DIST_CONST and not LSDisplacer.KKT:
            wRoads = np.full(len(self.roads_shapes), LSDisplacer.PDistRoads)
            weights.append(wRoads)
        return np.diag(np.concatenate(weights))


    # B = Y - S(Xcourant)
    def get_B(self): 
        b = None
        # inertia
        if LSDisplacer.ID_CONST:
            b = self.points_talus - self.x_courant
        # cross prod angles
        if LSDisplacer.ANGLES_CONST:
            angles_cross = self.angles_crossprod()
            if len(angles_cross) != 0:
                if b is None:
                    b = self.angles_ori - angles_cross
                else:
                    b = np.concatenate((b, (self.angles_ori - angles_cross)))
        # triangulation edges
        if LSDisplacer.EDGES_CONST:
            e_lengths = []
            for e in self.edges:
                e_lengths.append(self._edge_length(e)) #, self.x_courant))
            ee = self.e_lengths_ori - np.array(e_lengths)
            if b is None:
                b = ee
            else:
                b = np.concatenate((b, ee))
        # distance from roads
        if LSDisplacer.DIST_CONST:
            r_dists = []
            for i, r in enumerate(self.roads_shapes):
                fk = - self.dist_F_vectorized(r, i, self.x_courant[np.newaxis,:])
                r_dists.append(fk.item())
            if b is None:
                b = np.array(r_dists)
            else:
                b = np.concatenate((b, r_dists))
        return b


    def get_A(self): 
        a = None
        # inertia
        if LSDisplacer.ID_CONST:
            a = np.identity(self.nb_vars)
        # cross prod angles
        if LSDisplacer.ANGLES_CONST:
            angles = self.cross_norm_diff()
            if len(angles) != 0:
                if a is None:
                    a = angles
                else:
                    a = np.vstack((a, angles))
        # triangulation edges
        if LSDisplacer.EDGES_CONST:
            e_lens=[]
            for e in self.edges:
                ek = self.edge_length_diff(e)
                e_lens.append(ek)
            if a is None:
                a = np.array(e_lens)
            else:
                a = np.vstack((a, e_lens))
        # distance from roads
        if LSDisplacer.DIST_CONST:
            r_dists = []
            for i, r in enumerate(self.roads_shapes):
                fk = self.dist_F_diff(r, i)
                r_dists.append(fk)
            if a is None :
                a = np.array(r_dists)
            elif len(r_dists) != 0:
                a = np.vstack((a, r_dists))
        return a

    def compute_dx(self):
        if LSDisplacer.KKT:
            idsave, angsave, edgessave = LSDisplacer.ID_CONST, LSDisplacer.ANGLES_CONST, LSDisplacer.EDGES_CONST
            LSDisplacer.DIST_CONST = False
            A = self.get_A()
            B = self.get_B()
            LSDisplacer.DIST_CONST = True
            LSDisplacer.ID_CONST, LSDisplacer.ANGLES_CONST, LSDisplacer.EDGES_CONST = False, False, False
            C = self.get_A()
            D = self.get_B()
            LSDisplacer.ID_CONST, LSDisplacer.ANGLES_CONST, LSDisplacer.EDGES_CONST = idsave, angsave, edgessave
            loglsd.debug(f"A{A.shape} B{B.shape} C{C.shape} D{D.shape} P{self.P.shape}")
            atp = A.T @ self.P 
            atpa = atp @ A
            kkt = np.vstack((2 * atpa, C))
            kkt = np.hstack((kkt, np.vstack((C.T, np.zeros((len(self.roads_shapes), len(self.roads_shapes)))))))
            atpb = atp @ B
            kkt_b = np.concatenate((2 * atpb, D))
            dx = np.linalg.lstsq(kkt, kkt_b, rcond=None)
            return dx
        A = self.get_A()
        B = self.get_B()
        loglsd.debug(f"A{A.shape} B{B.shape} P{self.P.shape}")
        atp = A.T @ self.P 
        atpa = atp @ A
        atpb = atp @ B
        dx = np.linalg.lstsq(atpa, atpb, rcond=None)
        return dx

    def square(self):
        """
        returns a numpy array of the displaced points after the least square process
        """
        if LSDisplacer.KKT:
            loglsd.info("mode KKT")
            LSDisplacer.DIST_CONST = False
        alpha = 0.1
        ro = 0.1
        min_dx = np.inf
        norm_float = LSDisplacer.NORM_DX
        start_loop = timer()
        i = 0
        for i in range(LSDisplacer.MAX_ITER):
            start = timer()
            dx = self.compute_dx()
            if LSDisplacer.KKT:
                self.x_courant += alpha * dx[0][:self.nb_vars]
            else :
                self.x_courant += alpha * dx[0] #dx
            end = timer()
            normdx = np.linalg.norm(dx[0], ord=np.inf)
            alpha =  (LSDisplacer.H * ro) / (2**0.5 * normdx) if normdx != 0 else 0.1
            min_dx = normdx if normdx < min_dx else min_dx
            loglsd.info(f'iter {i}/ |dx| : {normdx:.4f} -- NORM_DXf: {norm_float} -- mean(dx): {np.mean(dx[0]):.2f} -- {(end - start):.2f}s per step')
            if normdx < norm_float : #NORM_DX :
                break
            if LSDisplacer.FLOATING_NORM:
                norm_float = LSDisplacer.NORM_DX if i < 100 else (LSDisplacer.NORM_DX + 2 * min_dx) / 3
        end_loop = timer()
        self.meta['nb_iters'], self.meta['time_s'], self.meta['dx_reached'] = i, (end_loop - start_loop), min_dx
        loglsd.warning(f'nb iterations: {i + 1} -- min |dx| reached: {min_dx} -- NORM_DXf: {norm_float} -- {(self.meta["time_s"]):.2f}s ')
        return self.x_courant