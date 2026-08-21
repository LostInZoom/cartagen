import numpy as np
from scipy.ndimage import convolve

def compute_curvature_and_slope(dem, cellsize):
    """
    Computes approximations of curvature (profile/plan) and slope.
    Based on finite differences (Wilson & Gallant, 2000).
    """
    # Gradients horizontaux et verticaux (premières dérivées)
    zy, zx = np.gradient(dem, cellsize)
    
    # Secondes dérivées
    zyy, zxy = np.gradient(zy, cellsize)
    zxx, _   = np.gradient(zx, cellsize)
    
    pente = np.sqrt(zx**2 + zy**2)
    
    # Pour éviter les divisions par zéro sur les zones parfaitement plates
    vis_pente = pente.copy()
    vis_pente[vis_pente == 0] = 1e-5
    
    # Approximation de la courbure maximale/minimale (Laplacien pour simplifier les crêtes/vallées)
    # Un Laplacien négatif indique une crête (convexe), positif une vallée (concave)
    laplacien = zxx + zyy
    
    return laplacien, pente

def calculate_hillshade(dem, azimuth, altitude, cell_size=1.0):
    """
    Calcule l'ombrage analytique d'un MNT en utilisant des opérations Numpy pures.
    Formule standard de Horn (utilisée par ArcGIS/QGIS).
    """
    # Conversion des angles en radians
    azimuth_rad = np.radians(360 - azimuth + 90) # Ajustement de la direction du soleil
    altitude_rad = np.radians(altitude)
    zenith_rad = np.pi / 2.0 - altitude_rad

    # Noyaux de convolution de Horn pour les gradients X et Y
    kernel_x = np.array([[-1, 0, 1],
                         [-2, 0, 2],
                         [-1, 0, 1]]) / (8.0 * cell_size)
    
    kernel_y = np.array([[ 1,  2,  1],
                         [ 0,  0,  0],
                         [-1, -2, -1]]) / (8.0 * cell_size)

    # Calcul des gradients verticaux (dz/dx) et horizontaux (dz/dy)
    dz_dx = convolve(dem, kernel_x, mode='reflect')
    dz_dy = convolve(dem, kernel_y, mode='reflect')

    # Calcul de la pente (slope) et de l'exposition (aspect) en radians
    slope = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
    
    # Calcul de l'aspect (gestion du diviseur par zéro avec np.arctan2)
    aspect = np.arctan2(dz_dy, -dz_dx)
    # Ajustement de l'aspect pour correspondre aux conventions géographiques
    aspect = np.where(aspect < 0, aspect + 2 * np.pi, aspect)

    # Formule mathématique de l'ombrage (Hillshade)
    shaded = (np.cos(zenith_rad) * np.cos(slope)) + \
             (np.sin(zenith_rad) * np.sin(slope) * np.cos(azimuth_rad - aspect))
             
    # Renvoyer des valeurs normalisées (0 à 1)
    return np.clip(shaded, 0, 1)