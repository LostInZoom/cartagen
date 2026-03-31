import geopandas as gpd
from shapely.geometry import LineString
from typing import Optional, List

from cartagen.algorithms.lines.bends import accordion, schematization
from cartagen.algorithms.lines.breaks import max_break, min_break
from cartagen.algorithms.lines.coalescence import coalescence_splitting
from cartagen.algorithms.propagation.network import network_propagation, network_propagation_batch
from cartagen.utils.lines.smoothing.gaussian import smooth_gaussian
from cartagen.utils.lines.smoothing.platre import smooth_platre

def galbe(
    network: gpd.GeoDataFrame,
    width: float,
    smoothing_method: str = 'platre',
    max_propagation: float = None,
    damping: float = 0.4,
    sigma: float = 30,
    inplace: bool = False
) -> gpd.GeoDataFrame:
    """
    Generalize a sinuous mountain road network.

    This process proposed by Mustière :footcite:p:`mustiere:2001-a`

    Parameters
    ----------
    network : GeoDataFrame of LineString
        The road network to generalize.
    width : float
        The width of the roads symbol.

    Notes
    -----
    This algorithm works best when using a properly made road network. The network must
    avoid as much as possible having nodes of degree 2, e.g. chunks of road which are
    not divided at intersections.
    """
    if not inplace:
        network = network.copy()

    for index, row in network.iterrows():
        geom = row.geometry
        chunks = coalescence_splitting(geom, width)


def galbe_single_line(
    line: LineString, 
    width: float, 
    smoothing_method: str = 'platre',
    sigma: float = 30,
    sample: int = None
) -> LineString:
    """
    Traite une ligne unique : Détection -> Accordéon/Schématisation -> Lissage.
    """

    # 1. Détection des zones de conflit (coalescence)
    chunks = coalescence_splitting(line, tolerance=width)
    processed_chunks = []
    
    for chunk in chunks:
        geom = chunk['geometry']
        level = chunk['coalescence']
        
        # Si pas de conflit (level 0), on garde la géométrie
        if level == 0:
            processed_chunks.append(geom)
            continue
            
        # 2. Transformation locale
        # On tente l'accordéon pour élargir les virages
        candidate = accordion(geom, width=width, sigma=sigma, sample=sample)
        
        # Si l'accordéon n'a pas assez d'espace (peu d'évolution), on schématise
        if candidate.length < geom.length * 1.02:
            candidate = schematization(candidate, sigma=sigma, sample=sample)
        
        # Petit lissage gaussien local pour nettoyer le chunk
        candidate = smooth_gaussian(candidate, sigma=sigma/3, sample=sample)
        processed_chunks.append(candidate)
    
    # 3. Réassemblage
    final_coords = []
    for i, g in enumerate(processed_chunks):
        coords = list(g.coords)
        if i > 0:
            final_coords.extend(coords[1:])
        else:
            final_coords.extend(coords)
    
    reassembled_line = LineString(final_coords)
    
    # 4. Lissage final (Le "Galbe")
    if smoothing_method == 'platre':
        # Platre est excellent pour la fluidité des réseaux routiers
        return smooth_platre(reassembled_line, threshold=width/4)
    else:
        return smooth_gaussian(reassembled_line, sigma=sigma, sample=sample)

def galbe_geodataframe(
    gdf: gpd.GeoDataFrame,
    road_width: float,
    max_propagation: float = 500,
    damping: float = 0.1,
    smoothing_method: str = 'platre',
    sigma: float = 30,
    inplace: bool = False
) -> gpd.GeoDataFrame:
    """
    Applique le processus GALBE à l'ensemble d'un GeoDataFrame.
    
    Args:
        gdf: Le GeoDataFrame de LineStrings à traiter.
        road_width: Largeur du symbole (ex: 15m pour une route principale).
        max_propagation: Distance max de diffusion du mouvement dans le réseau.
        damping: Facteur d'amortissement de la propagation.
        smoothing_method: 'platre' ou 'gaussian'.
        sigma: Paramètre pour le lissage gaussien.
        inplace: Si True, modifie le GDF d'entrée.
    """
    if not inplace:
        gdf = gdf.copy()
        
    modifications = []
    
    # Parcourir toutes les lignes pour calculer les corrections nécessaires
    for idx, row in gdf.iterrows():
        original_geom = row.geometry
        
        # On applique le coeur algorithmique
        new_geom = galbe_single_line(
            original_geom, 
            width=road_width, 
            smoothing_method=smoothing_method,
            sigma=sigma
        )
        
        # Si la géométrie a été modifiée, on l'ajoute à la liste des propagations
        if not original_geom.equals(new_geom):
            modifications.append((idx, new_geom))
            
    # Lancer la propagation réseau pour assurer la continuité aux intersections
    # Cette fonction utilise la logique de "fausses bordures" définie précédemment
    result_gdf = network_propagation_batch(
        gdf,
        modifications,
        max_propagation_distance=max_propagation,
        damping_factor=damping,
        inplace=True
    )
    
    return result_gdf