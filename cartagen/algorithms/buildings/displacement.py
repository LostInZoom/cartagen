# This is an implementation of the random building displacement algorithm
import random, math
import geopandas as gpd
from shapely.affinity import translate
from shapely import geometry, intersects, intersection
from shapely.strtree import STRtree
from cartagen.utils.partitioning.network import partition_networks

def random_displacement(
        polygons, networks=None, polygon_distance=10.0, network_distance=10.0,
        max_trials=25, max_displacement=10.0, network_partitioning=None
    ):
    """
    Iteratively displace polygons overlapping each other and the provided network.

    Displace the provided buildings if they overlap each other
    or are closer than the width value to the provided networks.
    This method is not deterministic, at each iteration, the algorithm
    select a random angle and a random distance to displace the polygon
    and see if it improves the overlapping issues. Running twice the
    algorithm returns two different solutions.
    This algorithm is succinctly presented in Wabiński *et al.*
    :footcite:p:`wabinski:2022`

    Parameters
    ----------
    polygons : GeoDataFrame of Polygon
        The buildings to displace.
    networks : list of GeoDataFrame of LineString, optional
        A list of networks the polygons need to be moved away from.
        If left to None, polygons will only be moved away from each other.
    polygon_distance : float, optional
        The minimum acceptable distance between polygons.
    network_distance : float, optional
        The minimum acceptable distance between the polygons
        and the provided networks.
    max_trials : int, optional
        The maximum number of trials before stopping the iteration.
        A trial represent the movement of one polygon that did not
        lower the mean overlapping area between polygons and networks.
    max_displacement : float, optional
        The maximum allowed distance of displacement per iteration.
    network_partioning : GeoDataFrame of LineString, optional
        The network to partition the data with. If provided, each network
        face is treated individually, thus improving performance on larger dataset.

    Returns
    -------
    GeoDataFrame

    See Also
    --------
    network_faces:
        Calculates the faces of one or multiple networks.

    References
    ----------
    .. footbibliography::
    """

    POLYGONS = polygons.to_dict('records')

    def __loop(indexes, network):
        """
        Launch a loop to iteratively displace polygons randomly
        and checks the congestion before validating
        """
        # Calculating the mean rate of overlapping polygons with other polygons and the provided network
        # It represents the mean congestion of polygons within the polygon block
        rate_mean = __get_buildings_overlapping_rate_mean(indexes, network)

        # Starting trial count
        trial = 0

        # Launching the loop which will displace buildings randomly as
        # long as the rate mean is above 0 and the max trial count is not exceeded
        while rate_mean > 0 and trial <= max_trials:
            # Selecting a random building index
            random_index = random.randint(0, len(indexes) - 1)
            random_building = indexes[random_index]

            overlap = __get_building_overlap(random_building, indexes, network)

            # Checking if that building is overlapping
            if overlap > 0:
                # Selecting a random angle (0-360)
                random_angle = random.uniform(0, 360)
                # Selecting a random length (0-max displacement variable)
                random_length = random.uniform(0, max_displacement)
                # Calculate displacement for x and y
                dx = math.cos(random_angle) * random_length
                dy = math.sin(random_angle) * random_length

                # Translating the building with the random values
                untranslated = POLYGONS[random_building]['geometry']
                translated = translate(untranslated, dx, dy)
                POLYGONS[random_building]['geometry'] = translated

                # Calulcating the new rate mean
                new_rate_mean = __get_buildings_overlapping_rate_mean(indexes, network)

                # If the new rate mean is equal or higher, we cancel the translation
                if (new_rate_mean >= rate_mean):
                    POLYGONS[random_building]['geometry'] = untranslated
                    trial += 1
                # Else, resetting the trial count and updating the rate mean
                else:
                    rate_mean = new_rate_mean
                    trial = 0

    def __get_buildings_overlapping_rate_mean(indexes, network):
        """
        Calculate the overlapping mean rate.
        """
        # If no polygons are provided, returns 0
        if len(indexes) == 0:
            return 0

        total = 0
        nb = 0
        # For each polygon, calculating the area overlapping other polygons,
        # and the network depending on a distance value
        for i in indexes:
            total += __get_building_overlap(i, indexes, network)
            nb += 1
        # Calculating the mean area
        if nb > 0:
            mean = total/nb

        return mean

    def __get_building_overlap(index, indexes, network):
        """
        Calculate the area of overlapping between polygons and the network
        """
        overlap = None
        b = POLYGONS[index]['geometry']
        buffered = b.buffer(polygon_distance)

        # For each buildings...
        for i in indexes:
            # Checking if it's not the same building
            if i != index:
                overlap = __get_overlapping_geometries(buffered, POLYGONS[i]['geometry'], overlap)
            
        # For each network section
        for n in network:
            overlap = __get_overlapping_geometries(b, n, overlap)

        # Returning the area of the geometry if it exists
        if (overlap is None) or (overlap.is_empty == True) or (overlap.area == 0):
            return 0
        else:
            return overlap.area

    def __get_overlapping_geometries(geom1, geom2, overlap):
        """
        Calculate the geometry of the intersection between a geographic object and a building
        """
        # If the building intersects the object
        if intersects(geom1, geom2):
            # Creating the intersection between the building and the object
            intersection = intersection(geom1, geom2)
            # If the geometry is empty, return the intersection...
            if overlap is None:
                return intersection
            # Else, returning the union between the intersection and the existing geometry
            else:
                return overlap.union(intersection)
        else:
            return overlap

    # Retrieve crs and convert polygons to list of dict records
    crs = polygons.crs

    network = []
    # Unpack the provided networks
    if networks is not None:
        for line in networks:
            # Append the buffered lines to the list
            for l in line.geometry:
                network.append(l.buffer(network_distance))

    if network_partitioning is not None:
        # Create the partitions -> tuple ([polygon index], [partition polygon geometry])
        partitions = partition_networks(polygons, *network_partitioning)

        # Loop through each partition
        for i, partition in enumerate(partitions[0]):
            partition_network = []

            if networks is not None:
                # Create the spatial index for the buffered network
                btree = STRtree(network)
                # Retrieve network sections that intersects the considered network face
                intersects = list(btree.query(partitions[1][i], predicate='intersects'))                
                # Append the buffered networks only if it truly intersects
                for inter in intersects:
                    if intersects(network[inter], partitions[1][i]):
                        partition_network.append(network[inter])

            # Launch the random displacement for the current partition
            __loop(partition, partition_network)
    else:
        # Launch the random displacement with all polygons and all buffered networks
        __loop(list(range(0, len(POLYGONS))), network)

    return gpd.GeoDataFrame(POLYGONS, crs=crs)


def displacement_buildings_simulated_annealing(
        polygons, networks=None, d_min1=10.0, d_min2=10.0,
        ppcost=1, plcost=5, max_d=5.0, lambda_=0.9, omega=40, zeta=20, psi=50, k=29
    ):
    """
    Iteratively displace polygons overlapping (or too close to) each other and the provided network.
    The displacement is optimised using a simulated annealing process, which is a probabilistic technique 
    for approximating the global optimum of a given function. The algorithm is based on the work of 
    Ware *et al.* :footcite:p:`ware_automated_2003`.

    Parameters
    ----------
    polygons : GeoDataFrame of Polygon
        The buildings to displace.
    networks : list of LineString, optional
        A list of networks the polygons need to be moved away from.
        If left to None, polygons will only be moved away from each other.
    d_min1 : float, optional
        The minimum acceptable distance between polygons.
    d_min2 : float, optional
        The minimum acceptable distance between polygons and networks.
    ppcost : int, optional
        The cost of overlapping polygons in the simulated annealing process. The higher the value, the more the algorithm will try to lower the mean overlapping area between polygons.
        lower the mean overlapping area between polygons and networks.
    plcost : int, optional
        The cost of overlapping polygons and networks in the simulated annealing process. The higher the value
    max_d : float, optional
        The maximum allowed distance of displacement per iteration.
    lambda_ : float, optional
        The cooling rate for the simulated annealing process.
    omega : int, optional
        The number of iterations for each temperature level.
    zeta : int, optional
        The number of times to restart the algorithm.
    psi : int, optional
        The number of times to accept a worse solution.
    k : int, optional
        The number of trial positions to consider at each iteration.

    Returns
    -------
    GeoDataFrame

    See Also
    --------
    random_displacement:
        Displace polygons randomly until no more improvements are found.

    Notes
    -----
    This algorithm is an optimisation process that tries many random displacements of the polygons before finding an optimal solution. It can very slow compared to other displacement algorithms, but it is more likely to find a better solution. 

    References
    ----------
    .. footbibliography::
    """
    class WareSimulatedAnnealing:
        def __init__(self, polygons, lines, params=None):
            self.M = polygons.copy(deep=True)  # Objets modifiables [cite: 12]
            self.F = lines     # Objets fixes [cite: 12]
            
            # Paramètres exacts de l'article (Section 2.5 & 2.6)
            self.d_min1 = 7.5   # Seuil PP [cite: 149]
            self.d_min2 = 7.5   # Seuil PL [cite: 149]
            self.ppcost = 1     # [cite: 149]
            self.plcost = 5    # [cite: 149]
            self.max_d = 5    # Distance d [cite: 149]
            
            # Paramètres de refroidissement (Table 1 & Section 2.6)
            self.lam = 0.9      # lambda [cite: 177]
            self.omega = 40     # alpha [cite: 177]
            self.zeta = 20      # [cite: 177]
            self.psi = 50       # [cite: 177]
            self.k = 29         # Nombre de positions d'essai (tp) [cite: 154]

            if params:
                self.d_min1 = params.get('d_min1')
                self.d_min2 = params.get('d_min2')
                self.ppcost = params.get('ppcost')
                self.plcost = params.get('plcost')
                self.max_d = params.get('d_max')
                self.lam = params.get('lambda')
                self.omega = params.get('omega')
                self.zeta = params.get('zeta')
                self.psi = params.get('psi')
                self.k = params.get('k')
            
            # État actuel : liste de [géométrie_translatée, index_tp_actuel]
            self.current_positions = [0] * len(self.M) 
            self.trial_vectors = self._build_trial_vectors(self.k) # k=29 [cite: 154]

        def _build_trial_vectors(self, k):
            """Génère les vecteurs tp selon la Figure 2[cite: 46, 47]."""
            vectors = [(0, 0)] # tp1: position originale [cite: 15]
            for i in range(k - 1):
                angle = 2 * math.pi * i / (k - 1)
                # On distribue les positions radialement jusqu'à d [cite: 16, 18]
                dist = self.max_d if i % 2 == 0 else self.max_d * 0.5
                vectors.append((dist * math.cos(angle), dist * math.sin(angle)))
            return vectors

        def get_obj_cost(self, obj_idx, tp_idx, all_positions):
            cost = 0
            # On génère la géométrie candidate
            poly_i = translate(self.M.geometry.iloc[obj_idx], *self.trial_vectors[tp_idx])
            
            # --- Conflit PP (Polygone-Polygone) ---
            for j, pos_j_idx in enumerate(all_positions):
                if obj_idx == j: continue
                poly_j = translate(self.M.geometry.iloc[j], *self.trial_vectors[pos_j_idx])
                
                dist = poly_i.distance(poly_j)
                if dist < self.d_min1:
                    # Si intersection réelle : coût très élevé (non explicite mais nécessaire)
                    if poly_i.intersects(poly_j):
                        cost += self.ppcost * 50 
                    else:
                        # Si simple proximité : coût standard [cite: 24]
                        cost += self.ppcost 
            
            # --- Conflit PL (Polygone-Ligne) ---
            for line in self.F:
                dist_line = poly_i.distance(line)
                if dist_line < self.d_min2:
                    if poly_i.intersects(line):
                        cost += self.plcost * 50 # On ne traverse pas la route !
                    else:
                        cost += self.plcost # [cite: 49]
                        
            return cost

        def run(self):
            n = len(self.M) # Nombre de bâtiments (16 dans votre cas)
            D_curr = list(self.current_positions)
            
            # 1. Calibration (la boucle de 500 n'est QUE pour T_initial) [cite: 172]
            T = self._calibrate_T(D_curr)
            
            # 2. Boucle de refroidissement (psi = 50 paliers max) 
            for stage in range(self.psi):
                moves_attempted = 0
                successes = 0
                
                # On teste jusqu'à (omega * n) positions à cette température [cite: 144, 177]
                # Pour 16 bâtiments, omega=100 -> 1600 essais par palier.
                max_moves = self.omega * n 
                target_successes = self.zeta * n # ex: 30 * 16 
                
                while moves_attempted < max_moves:
                    idx = random.randrange(n)
                    new_tp = random.randrange(len(self.trial_vectors))
                    
                    # Calcul du Delta E (variation de conflit) [cite: 113]
                    old_cost = self.get_obj_cost(idx, D_curr[idx], D_curr)
                    new_cost = self.get_obj_cost(idx, new_tp, D_curr)
                    delta_E = old_cost - new_cost 
                    
                    # Critère d'acceptation de Metropolis [cite: 115, 118]
                    if delta_E > 0 or random.random() < math.exp(delta_E / T):
                        D_curr[idx] = new_tp
                        successes += 1
                        
                        # Si on a trouvé beaucoup de bonnes positions, on refroidit plus vite [cite: 145]
                        if successes >= target_successes:
                            break
                    
                    moves_attempted += 1
                    
                # 3. Mise à jour de la température (Géométrique) 
                T *= self.lam 
                
                # Si aucun mouvement n'a été accepté, l'algorithme a convergé (ou est bloqué) [cite: 146]
                if successes == 0:
                    break
            
            # Update internal state with final solution
            self.current_positions = D_curr
            
            # Create a copy of the GeoDataFrame and update geometries with translations
            result_gdf = self.M.copy()
            translated_geometries = [translate(self.M.geometry.iloc[i], *self.trial_vectors[D_curr[i]]) for i in range(n)]
            result_gdf.geometry = translated_geometries
            
            return result_gdf

        def _calibrate_T(self, state):
            """Méthode de la section 2.6 pour tau[cite: 173, 175]."""
            deltas = []
            for _ in range(500):
                i = random.randrange(len(self.M))
                ntp = random.randrange(len(self.trial_vectors))
                diff = self.get_obj_cost(i, ntp, state) - self.get_obj_cost(i, state[i], state)
                if diff > 0: deltas.append(diff)
            if not deltas: return 1.0
            m_delta = sum(deltas) / len(deltas)
            return m_delta / math.log(3) # [cite: 175]
    
    params = {
        'd_min1': d_min1,
        'd_min2': d_min2,
        'ppcost': ppcost,
        'plcost': plcost,
        'd_max': max_d,
        'lambda': lambda_,
        'omega': omega,
        'zeta': zeta,
        'psi': psi,
        'k': k
    }

    generalizer = WareSimulatedAnnealing(polygons, networks, params)
    output_dataframe = generalizer.run()

    return output_dataframe