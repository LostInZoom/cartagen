import numpy as np
import rasterio
from scipy.ndimage import uniform_filter, gaussian_filter
from cartagen.utils.raster.DTM import compute_curvature_and_slope

def sculpt_dtm(input_path, output_path, 
                    base_filter_size=15, 
                    curvature_filter_size=9, 
                    exaggeration_factor=1.5, 
                    excavation_factor=1.5, 
                    slope_threshold=0.05):
    """
    Implementation of the "Terrain Sculptor" algorithm from Leonowicz *et al.* :footcite:p:`leonowicz:2010`.
    This algorithm generalizes a Digital Terrain Model (DTM) by smoothing the terrain while preserving significant 
    features such as ridges and valleys. It uses curvature and slope information to create a more visually appealing 
    representation of the terrain.
    
    Parameters
    ----------
    input_path : str
        The path to the input DEM file.
    output_path : str
        The path where the generalized DEM will be saved.
    base_filter_size : int, optional
        The size of the filter for the base surface smoothing.
    curvature_filter_size : int, optional
        The size of the filter for curvature detection.
    exaggeration_factor : float, optional
        The factor for exaggerating the terrain features.
    excavation_factor : float, optional
        The factor for deepening the terrain features.
    slope_threshold : float, optional
        The threshold for slope-based classification.
    interval : float, optional
        Interval used for the interpolation during the initial displacement vector calculation.
    
    Returns
    -------
    None :
        The function modifies the input DEM in place and saves the result to the specified output path.

    Notes
    -------
    The terrain sculptor algorithm is also available as a free and open Java tool :footcite:p:`leonowicz:2010`.
    
    See Also
    --------
    douglas_peucker_3d_DTM : Another method for DTM generalization using the Douglas-Peucker algorithm in 3D.

    References
    ----------
    .. footbibliography::
    """

    with rasterio.open(input_path) as src:
        dem = src.read(1).astype(np.float64)
        meta = src.meta.copy()
        cellsize = src.res[0] # Récupère la résolution du pixel (ex: 25m)
        nodata = src.nodata

    # Masque pour les valeurs NoData
    masque_valide = (dem != nodata) if nodata is not None else np.ones_like(dem, dtype=bool)
    
    # --- ÉTAPE 1 : MNT lissé (Base surface) ---
    # Filtrage par moyenne mobile (low-pass mean filter)
    dem_lisse = uniform_filter(dem, size=base_filter_size)
    
    # --- ÉTAPE 2 : Détection des structures par courbure ---
    # On lisse le MNT d'origine avant de calculer les courbures pour enlever le bruit
    dem_pre_courbure = uniform_filter(dem, size=curvature_filter_size)
    courbure, _ = compute_curvature_and_slope(dem_pre_courbure, cellsize)
    
    # --- ÉTAPE 3 : Modèles exagérés et creusés ---
    # On crée une surface accentuée vers le haut (crêtes) et vers le bas (vallées)
    dem_exagere = dem_lisse + (dem - dem_lisse) * exaggeration_factor
    dem_creuse = dem_lisse - (dem_lisse - dem) * excavation_factor
    
    # --- ÉTAPE 4 : Création des poids et des modèles Montagne / Plaine ---
    # Normalisation basique des courbures pour créer des grilles de poids [0, 1]
    # Crêtes : Courbure négative (Laplacien < 0)
    poids_cretes = np.clip(-courbure / np.std(courbure), 0, 1)
    # Vallées : Courbure positive (Laplacien > 0)
    poids_vallees = np.clip(courbure / np.std(courbure), 0, 1)
    
    # Combinaison (Step 4)
    model_montagne = poids_cretes * dem_exagere + (1 - poids_cretes) * dem_lisse
    model_plaine = poids_vallees * dem_creuse + (1 - poids_vallees) * dem_lisse
    
    # --- ÉTAPE 5 : Recombinaison finale via la Pente (Slope mask) ---
    # Calcul de la pente sur le modèle lissé pour séparer plaine et montagne
    _, pente_lisse = compute_curvature_and_slope(dem_lisse, cellsize)
    # Filtre passe-bas sur la pente pour adoucir la transition (Figure 3G)
    pente_lisse_floue = gaussian_filter(pente_lisse, sigma=5)
    
    # Poids de transition basé sur le seuil de pente (sigmoïde ou linéaire)
    # Si pente > seuil_pente -> plutôt montagne (1), si inférieure -> plutôt plaine (0)
    poids_pente = np.clip(pente_lisse_floue / slope_threshold, 0, 1)
    
    # Modèle final combiné
    mnt_final = poids_pente * model_montagne + (1 - poids_pente) * model_plaine
    
    # Rétablir les NoData d'origine si nécessaire
    if nodata is not None:
        mnt_final[~masque_valide] = nodata
        
    # --- ÉTAPE 6 : Export du MNT Généralisé ---
    meta.update(dtype=rasterio.float64)
    with rasterio.open(output_path, "w", **meta) as dst:
        dst.write(mnt_final, 1)
        
