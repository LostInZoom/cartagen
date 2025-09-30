import numpy as np
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import nearest_points
import math
import matplotlib.pyplot as plt
from itertools import combinations

class ElasticBeamPolyline:
    """
    Implements the elastic beams algorithm for polyline displacement.
    Based on the paper "Cartographic Displacement in Generalization: Introducing Elastic Beams"
    by M. Bader and M. Barrault.
    """

    def __init__(self, polyline: LineString, E=1.0, A=1.0, I=1.0):
        """
        Initializes the beam structure from a Shapely polyline.

        Args:
            polyline: The Shapely polyline to be displaced.
            E: Modulus of elasticity (material constant).
            A: Cross-sectional area of the beam.
            I: Moment of inertia of the beam.
        """
        self.original_coords = np.array(polyline.coords)
        self.current_coords = self.original_coords.copy()
        self.num_nodes = len(self.original_coords)
        self.E = E
        self.A = A
        self.I = I

        # Each node has 3 degrees of freedom: x displacement, y displacement, and rotation.
        self.dof = 3 * self.num_nodes
        self.K_global = np.zeros((self.dof, self.dof))

    def _get_local_stiffness_matrix(self, p1, p2):
        """
        Calculates the local stiffness matrix for a single beam element.
        Based on equation (11) from the paper.
        """
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        L = np.sqrt(dx**2 + dy**2)

        if L == 0:
            return np.zeros((6, 6))

        c = dx / L
        s = dy / L
        
        # Coefficients to simplify the calculation
        C1 = self.A * c**2 + (12 * self.I / L**2) * s**2
        C2 = self.A * s**2 + (12 * self.I / L**2) * c**2
        C3 = (self.A - 12 * self.I / L**2) * c * s
        C4 = (6 * self.I / L) * s
        C5 = (6 * self.I / L) * c
        
        # Local stiffness matrix in a global coordinate system (Eq. 11)
        K_local = np.array([
            [C1, C3, -C4, -C1, -C3, -C4],
            [C3, C2, C5, -C3, -C2, C5],
            [-C4, C5, 4 * self.I, C4, -C5, 2 * self.I],
            [-C1, -C3, C4, C1, C3, C4],
            [-C3, -C2, -C5, C3, C2, -C5],
            [-C4, C5, 2 * self.I, C4, -C5, 4 * self.I]
        ])
        
        return (self.E / L) * K_local

    def _assemble_global_matrix(self):
        """
        Assembles the global stiffness matrix from all beam elements.
        """
        self.K_global = np.zeros((self.dof, self.dof))
        for i in range(self.num_nodes - 1):
            p1 = self.current_coords[i]
            p2 = self.current_coords[i+1]
            
            K_local = self._get_local_stiffness_matrix(p1, p2)
            
            # Degrees of freedom (DOF) indices for the current element
            dof_indices = [3*i, 3*i+1, 3*i+2, 3*(i+1), 3*(i+1)+1, 3*(i+1)+2]
            
            # Assembly of the global matrix (K_global)
            for r in range(6):
                for c in range(6):
                    self.K_global[dof_indices[r], dof_indices[c]] += K_local[r, c]

    def _solve_for_displacement(self, f, gamma=0.5, mu=1.0, fix_boundaries=True, iterations=100, verbose=False):
        """
        Iteratively solves the system to calculate displacement.
        Based on equation (13) from the paper.

        Args:
            f: Vector of external forces.
            gamma: Time step / inertia factor.
            mu: Force multiplier.
            fix_boundaries: True to fix the first and last node.
            iterations: Number of iterations for convergence.
            verbose: If True, prints iteration details.
        """
        if verbose:
            print(f"Maximum force norm: {np.max(np.abs(f))}")
        
        # Boundary conditions: fix the first and last node.
        # We fix the first 6 and last 6 DOFs to 0.
        fixed_dof = {0, 1, 2, self.dof - 3, self.dof - 2, self.dof - 1}

        # Initial displacement vector
        d_prev = np.zeros(self.dof)

        for t in range(iterations):
            self._assemble_global_matrix()
            
            # Create the system matrix (I + gamma*K)
            system_matrix = np.eye(self.dof) + gamma * self.K_global
            
            # Right-hand side vector (d_prev + gamma * mu * f)
            rhs_vector = d_prev + gamma * mu * f
            
            # Apply boundary conditions
            if fix_boundaries:
                for idx in fixed_dof:
                    system_matrix[idx, :] = 0
                    system_matrix[idx, idx] = 1
                    rhs_vector[idx] = 0
                
            # Solve the system of equations
            try:
                d_curr = np.linalg.solve(system_matrix, rhs_vector)
            except np.linalg.LinAlgError:
                print("Matrix singularity error. Stopping iteration.")
                break
            
            # Update node positions
            for i in range(self.num_nodes):
                self.current_coords[i, 0] += d_curr[3*i]
                self.current_coords[i, 1] += d_curr[3*i+1]
            
            d_prev = d_curr
            if verbose:
                print(f"Iteration step {t+1}/{iterations}: Maximum displacement: {np.max(np.abs(d_curr))}")

        return self.current_coords
    

def _calculate_proximity_forces(polyline1: LineString, polyline2: LineString, min_dist: float):
    """
    Calcule les vecteurs de forces de proximité pour deux polylignes.
    Basé sur l'équation (16) de l'article Bader & Barrault (2001).

    Args:
        polyline1: La première polyligne Shapely.
        polyline2: La deuxième polyligne Shapely.
        min_dist: the minimum distance between the two lines.

    Returns:
        Deux vecteurs NumPy : (f1, f2), où f1 et f2 sont les vecteurs de forces
        (en x, y et rotation) pour chaque nœud de chaque polyligne.
    """
    coords1 = np.array(polyline1.coords)
    coords2 = np.array(polyline2.coords)
    
    # Chaque nœud a 3 degrés de liberté (x, y, rotation).
    num_nodes1 = len(coords1)
    num_nodes2 = len(coords2)
    
    f1 = np.zeros(3 * num_nodes1)
    f2 = np.zeros(3 * num_nodes2)

    # Calcul des forces de la polyligne 2 sur la polyligne 1
    for i in range(num_nodes1):
        node1_geom = Point(coords1[i])

        # Trouver le point le plus proche sur la deuxième polyligne
        point1, point2 = nearest_points(node1_geom, polyline2)
        
        # Calcul de la distance
        distance = point1.distance(point2)
        
        # Vecteur pointant du point2 vers le point1 (force de répulsion)
        if distance < min_dist and distance > 1e-9:
            normal_vector = np.array([point1.x - point2.x, point1.y - point2.y])
            force_magnitude = (min_dist - min(distance, min_dist)) / (distance * min_dist)
            force_vector = normal_vector * force_magnitude
        else:
            force_vector = np.array([0.0, 0.0])
        
        # Appliquer la force au nœud i (déplacements en x et y)
        f1[3*i] = force_vector[0]
        f1[3*i+1] = force_vector[1]

    # Calcul des forces de la polyligne 1 sur la polyligne 2 (réciproque)
    for i in range(num_nodes2):
        node2_geom = Point(coords2[i])
        
        # Trouver le point le plus proche sur la première polyligne
        point2, point1 = nearest_points(node2_geom, polyline1)
        
        distance = point2.distance(point1)
        
        # Vecteur pointant du point1 vers le point2 (force de répulsion)
        if distance < min_dist and distance > 1e-9:
            normal_vector = np.array([point2.x - point1.x, point2.y - point1.y])
            force_magnitude = (min_dist - min(distance, min_dist)) / (distance * min_dist)
            force_vector = normal_vector * force_magnitude
        else:
            force_vector = np.array([0.0, 0.0])
        
        f2[3*i] = force_vector[0]
        f2[3*i+1] = force_vector[1]

    return f1, f2


def beams_displacement(lines_geodataframe, min_dist, E=1.0, A=1.0, I=1.0, gamma=0.5, mu=0.1, fix_boundaries=True, verbose=False, iterations=100):
    """
    Applies the elastic beams algorithm to a geodataframe of polyline features. The forces are computed based on proximity between those line features.

    Args:
        polylines: List of Shapely polylines to be displaced.
        min_dist: The minimum distance that should separate each polyline.
        E: Modulus of elasticity (material constant).
        A: Cross-sectional area of the beam.
        I: Moment of inertia of the beam.
        gamma: Time step / inertia factor.
        mu: Force multiplier.
        fix_boundaries: True to fix the first and last node of each polyline.
        verbose: If True, prints iteration details.
        iterations: Number of iterations for convergence.

    Returns:
        A new GeoDataFrame of displaced polylines.
    """
    displaced_lines = lines_geodataframe.copy()

    # create a dictionary to hold the forces
    forces_dict = {idx: np.zeros(3 * len(line.geometry.coords)) for idx, line in displaced_lines.iterrows()}
    
    for line1, line2 in combinations(displaced_lines.iterrows(), 2):
        # Calculate proximity forces between the two polylines
        id1 = line1[0]
        id2 = line2[0]
        f1, f2 = _calculate_proximity_forces(line1[1]['geometry'], line2[1]['geometry'], min_dist)

        # Accumulate forces for each polyline
        if verbose:
            print("Forces on line", id1, ":", f1)
            print("Forces on line", id2, ":", f2)
        forces_dict[id1] += f1
        forces_dict[id2] += f2

    for idx, line in displaced_lines.iterrows():    
        f = forces_dict[idx]
        if np.all(f == 0):
            continue  # No forces to apply, skip displacement

        if verbose:
            print(f"Maximum force norm: {np.max(np.abs(f))}")
        # Initialize elastic beam objects for the polyline
        beam = ElasticBeamPolyline(line.geometry, E=E, A=A, I=I)
        
        # Solve for displacements
        displaced_coords = beam._solve_for_displacement(f, gamma=gamma, mu=mu, fix_boundaries=fix_boundaries, iterations=iterations, verbose = verbose)
        
        # Update the geometries in the GeoDataFrame
        displaced_lines.at[idx, 'geometry'] = LineString(displaced_coords)

    return displaced_lines


def _test_beams():
    # Create a simple polyline for testing
    coords = [(0, 0), (10, 0), (20, 5), (30, 5)]
    polyline = LineString(coords)
    
    # Initialize the elastic beam object
    beam = ElasticBeamPolyline(polyline, E=1e6, A=10, I=1)
    
    # Define external forces.
    # In this case, we apply a repulsive force to an internal node
    # to simulate a conflict. Each node has 3 DOFs (x, y, rotation).
    f = np.zeros(beam.dof)
    # Apply an upward force (positive in y) on the third node (index 2)
    # The index 3*2+1 corresponds to the y displacement of the third node.
    force_magnitude = 50.0
    f[3*2+1] = force_magnitude 
    
    # Solve the system iteratively
    displaced_coords = beam._solve_for_displacement(f, gamma=0.1, mu=1.0, iterations=50)

    # Convert the displaced coordinates to a Shapely polyline
    displaced_polyline = LineString(displaced_coords)
    
    # Visualize the results
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    # Display the original polyline
    x_orig, y_orig = zip(*beam.original_coords)
    ax.plot(x_orig, y_orig, 'b-o', label='Original Polyline')
    
    # Display the displaced polyline
    x_disp, y_disp = zip(*displaced_coords)
    ax.plot(x_disp, y_disp, 'r-o', label='Displaced Polyline')
    
    # Display the force vector
    node_index = 2
    ax.arrow(coords[node_index][0], coords[node_index][1],
             0, force_magnitude/10,  # The 10 is just a scaling factor for visualization
             head_width=0.5, head_length=0.5, fc='k', ec='k', label='Applied Force')
    
    ax.set_title("Polyline Displacement with Elastic Beams Algorithm")
    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")
    ax.legend()
    ax.grid(True)
    ax.set_aspect('equal', adjustable='box')
    plt.show()

def _test_forces():
    # Créer deux polylignes proches
    line_A = LineString([(0, 0), (10, 0), (20, 5)])
    line_B = LineString([(0.5, 2), (10.5, 2), (20.5, 7)])
    
    # Définir les paramètres des forces de proximité
    threshold_dist = 3.0 # d0
    
    # Calculer les forces
    forces_A, forces_B = _calculate_proximity_forces(line_A, line_B, threshold_dist)
    
    print("Forces calculées pour la ligne A :")
    print(forces_A)
    print("\nForces calculées pour la ligne B :")
    print(forces_B)
    
    # Visualisation des forces (optionnel, pour montrer la direction)
    fig, ax = plt.subplots()
    
    x_a, y_a = line_A.xy
    x_b, y_b = line_B.xy
    
    ax.plot(x_a, y_a, 'b-o', label='Ligne A')
    ax.plot(x_b, y_b, 'r-o', label='Ligne B')

    # Afficher les vecteurs de force
    for i in range(len(line_A.coords)):
        fx, fy = forces_A[3*i], forces_A[3*i+1]
        ax.arrow(line_A.coords[i][0], line_A.coords[i][1], fx, fy, head_width=0.2, color='blue')
    
    for i in range(len(line_B.coords)):
        fx, fy = forces_B[3*i], forces_B[3*i+1]
        ax.arrow(line_B.coords[i][0], line_B.coords[i][1], fx, fy, head_width=0.2, color='red')
        
    ax.set_title("Forces de Proximité entre deux Polylignes")
    ax.set_aspect('equal', adjustable='box')
    plt.legend()
    plt.grid(True)
    plt.show()
 

# Example of use
if __name__ == '__main__':
    _test_forces()