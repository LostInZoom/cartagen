# import geopandas as gpd
# import numpy as np
# from shapely.geometry import LineString
# from shapely import hausdorff_distance
# from typing import Optional, List

# from cartagen.algorithms.lines.bends import accordion, schematization
# from cartagen.algorithms.lines.breaks import max_break, min_break
# from cartagen.algorithms.lines.coalescence import coalescence_splitting
# from cartagen.algorithms.propagation.network import propagation_network, propagation_network_batch
# from cartagen.utils.lines.smoothing.gaussian import smooth_gaussian
# from cartagen.utils.lines.smoothing.platre import smooth_platre

# def galbe(
#     network, width, hausdorff=100.0, max_propagation=None, damping=0.4, sigma=30,
#     break_exaggeration=1.0, accordion_exaggeration=1.0, inplace=False
#     ):
#     """
#     Generalize a sinuous mountain road network.

#     This process proposed by Mustière :footcite:p:`mustiere:2001-a`

#     Parameters
#     ----------
#     network : GeoDataFrame of LineString
#         The road network to generalize.
#     width : float
#         The width of the roads symbol.

#     Notes
#     -----
#     This algorithm works best when using a properly made road network. The network must
#     avoid as much as possible having nodes of degree 2, e.g. chunks of road which are
#     not divided at intersections.
#     """
#     if not inplace:
#         network = network.copy()

#     # Iterate over individual lines
#     for index, row in network.iterrows():
#         geom = row.geometry
#         # Calculate the coalescence of the line
#         lines, coalescences = coalescence_splitting(geom, width)

#         final = []
#         for i, line in enumerate(lines):
#             coalescence = coalescences[i]
            
#             # If no coalescence is found, apply only a gaussian smoothing
#             if coalescence == 'None':
#                 final.append(smooth_gaussian(line, sigma=sigma))
#                 continue
            
#             # If coalescence is found on one side, apply max break to exaggerate the bend
#             if coalescence == 'Left' or coalescence == 'Right':
#                 final.append(max_break(line, width, exaggeration=break_exaggeration))
#                 continue

#             # Here, coalescence is found on both sides
#             # Start by applying accordion to 
#             spread = accordion(line, width, exaggeration=accordion_exaggeration, sigma=sigma)

#             hd = hausdorff_distance(line, spread)
#             count = 0
#             max_iter = 50

#             while hd > hausdorff and count < max_iter:
#                 unbend = schematization(line, sigma=sigma)
#                 spread = accordion(line, width, exaggeration=accordion_exaggeration, sigma=sigma)
#                 count += 1
            
#             chunks, c = coalescence_splitting()
    
#     return None

# def galbe_single_line(
#     line: LineString, 
#     width: float, 
#     smoothing_method: str = 'platre',
#     sigma: float = 30,
#     sample: int = None
# ) -> LineString:
#     """
#     Traite une ligne unique : Détection -> Accordéon/Schématisation -> Lissage.
#     """

#     # 1. Détection des zones de conflit (coalescence)
#     chunks = coalescence_splitting(line, tolerance=width)
#     processed_chunks = []
    
#     for chunk in chunks:
#         geom = chunk['geometry']
#         level = chunk['coalescence']
        
#         # Si pas de conflit (level 0), on garde la géométrie
#         if level == 0:
#             processed_chunks.append(geom)
#             continue
            
#         # 2. Transformation locale
#         # On tente l'accordéon pour élargir les virages
#         candidate = accordion(geom, width=width, sigma=sigma, sample=sample)
        
#         # Si l'accordéon n'a pas assez d'espace (peu d'évolution), on schématise
#         if candidate.length < geom.length * 1.02:
#             candidate = schematization(candidate, sigma=sigma, sample=sample)
        
#         # Petit lissage gaussien local pour nettoyer le chunk
#         candidate = smooth_gaussian(candidate, sigma=sigma/3, sample=sample)
#         processed_chunks.append(candidate)
    
#     # 3. Réassemblage
#     final_coords = []
#     for i, g in enumerate(processed_chunks):
#         coords = list(g.coords)
#         if i > 0:
#             final_coords.extend(coords[1:])
#         else:
#             final_coords.extend(coords)
    
#     reassembled_line = LineString(final_coords)
    
#     # 4. Lissage final (Le "Galbe")
#     if smoothing_method == 'platre':
#         # Platre est excellent pour la fluidité des réseaux routiers
#         return smooth_platre(reassembled_line, threshold=width/4)
#     else:
#         return smooth_gaussian(reassembled_line, sigma=sigma, sample=sample)

# def galbe_geodataframe(
#     gdf: gpd.GeoDataFrame,
#     road_width: float,
#     max_propagation: float = 500,
#     damping: float = 0.1,
#     smoothing_method: str = 'platre',
#     sigma: float = 30,
#     inplace: bool = False
# ) -> gpd.GeoDataFrame:
#     """
#     Applique le processus GALBE à l'ensemble d'un GeoDataFrame.
    
#     Args:
#         gdf: Le GeoDataFrame de LineStrings à traiter.
#         road_width: Largeur du symbole (ex: 15m pour une route principale).
#         max_propagation: Distance max de diffusion du mouvement dans le réseau.
#         damping: Facteur d'amortissement de la propagation.
#         smoothing_method: 'platre' ou 'gaussian'.
#         sigma: Paramètre pour le lissage gaussien.
#         inplace: Si True, modifie le GDF d'entrée.
#     """
#     if not inplace:
#         gdf = gdf.copy()
        
#     modifications = []
    
#     # Parcourir toutes les lignes pour calculer les corrections nécessaires
#     for idx, row in gdf.iterrows():
#         original_geom = row.geometry
        
#         # On applique le coeur algorithmique
#         new_geom = galbe_single_line(
#             original_geom, 
#             width=road_width, 
#             smoothing_method=smoothing_method,
#             sigma=sigma
#         )
        
#         # Si la géométrie a été modifiée, on l'ajoute à la liste des propagations
#         if not original_geom.equals(new_geom):
#             modifications.append((idx, new_geom))
            
#     # Lancer la propagation réseau pour assurer la continuité aux intersections
#     # Cette fonction utilise la logique de "fausses bordures" définie précédemment
#     result_gdf = network_propagation_batch(
#         gdf,
#         modifications,
#         max_propagation_distance=max_propagation,
#         damping_factor=damping,
#         inplace=True
#     )
    
#     return result_gdf

import numpy as np
import geopandas as gpd
from shapely import hausdorff_distance
from shapely.geometry import LineString
from typing import Tuple, List, Optional

from cartagen.algorithms.lines.bends import accordion, schematization
from cartagen.algorithms.lines.breaks import max_break, min_break
from cartagen.algorithms.lines.coalescence import coalescence_splitting
from cartagen.algorithms.propagation.network import propagation_network, propagation_network_batch
from cartagen.utils.lines.smoothing.gaussian import smooth_gaussian
from cartagen.utils.lines.smoothing.platre import smooth_platre


def chain_line_parts(parts: List[LineString]) -> LineString:
    """
    Chain multiple line parts back into a single line.
    
    Parameters
    ----------
    parts : List[LineString]
        List of line parts to chain together.
    
    Returns
    -------
    LineString
        The chained line.
    """
    if not parts:
        return LineString()
    if len(parts) == 1:
        return parts[0]
    
    all_coords = []
    for part in parts:
        coords = list(part.coords)
        if all_coords and coords[0] == all_coords[-1]:
            # Avoid duplicate point at junction
            coords = coords[1:]
        all_coords.extend(coords)
    
    return LineString(all_coords)


def galbe(
    network: gpd.GeoDataFrame,
    width: float,
    exaggeration: float = 1.0,
    legibility_factor: float = 1.7,
    sigma_smooth: float = 5,
    sigma_platre: float = 1,
    hausdorff_threshold: Optional[float] = 100.0,
    sample: Optional[int] = None,
    propagation_distance: Optional[float] = None,
    damping: float = 0.1,
    tolerance: float = 1e-6,
    inplace: bool = False
) -> gpd.GeoDataFrame:
    """
    Generalize sinuous mountain roads.
    
    GALBE (Généralisation Adaptative du Linéaire Basée sur l'Empâtement) is an
    adaptive generalization process for linear features based on coalescence detection.
    The process was proposed by Mustière (1998) for generalizing sinuous mountain roads.
    
    The algorithm follows these steps:
    1. Detect coalescence (pastiness) on the line
    2. Split the line into homogeneous parts based on coalescence type
    3. Apply appropriate algorithms to each part:
        - No coalescence: Gaussian smoothing
        - One-sided coalescence: Min or Max break (depending on context)
        - Two-sided coalescence: Accordion (with optional Schematization)
    4. Propagate displacements after each destructive algorithm
    5. Chain parts back together
    6. Apply light Platre smoothing to remove artifacts
    7. Propagate final displacements through the network
    
    Parameters
    ----------
    network : GeoDataFrame
        The road network to generalize. Must contain LineString geometries.
    width : float
        The width in meters of the road symbol casing (used for coalescence detection).
    exaggeration : float, optional
        Multiplicator for the width parameter in break and accordion algorithms.
        Default is 1.0.
    legibility_factor : float, optional
        Factor used in coalescence detection to group nearby coalescence chunks.
        Default is 1.7.
    sigma_smooth : float, optional
        Gaussian smoothing strength. Default is 5.
    sigma_platre : float, optional
        Platre smoothing strength for final smoothing. Default is 1.
    hausdorff_threshold : float, optional
        Maximum acceptable Hausdorff distance when applying accordion.
        If exceeded, schematization is applied first. If None, no threshold is used.
        Default is 100.0.
    sample : int, optional
        Gaussian smoothing sample size.
    propagation_distance : float, optional
        Maximum distance for displacement propagation through the network.
    damping : float, optional
        Damping factor for network propagation. Default is 0.1.
    tolerance : float, optional
        Tolerance for considering roads as touching in network propagation.
        Default is 1e-6.
    inplace : bool, optional
        Whether to modify the GeoDataFrame in place. Default is False.
    
    Returns
    -------
    GeoDataFrame
        The generalized road network.
    
    References
    ----------
    Mustière, S. (1998). Processus de généralisation cartographique de routes à 
    l'aide de déformations de Moindres Carrés. Bulletin du Comité Français de 
    Cartographie, 156, 45-52.
    
    Examples
    --------
    >>> roads = gpd.read_file("mountain_roads.shp")
    >>> generalized = galbe(roads, width=10.0, exaggeration=1.5)
    """
    if not inplace:
        network = network.copy()
    
    final_modifications = []
    
    for idx, row in network.iterrows():
        line = row.geometry
        
        if not isinstance(line, LineString) or line.is_empty:
            continue
        
        original_line = line
        
        # Step 1: Detect coalescence and split the line
        parts, coalescence_types = coalescence_splitting(line, width, legibility_factor)
        
        # Step 2: Process parts based on coalescence type
        # Order: two-sided, one-sided, none
        processed_parts = list(parts)  # Start with original parts
        
        # Sort parts by priority: 'both' > 'left'/'right' > 'none'
        priority_map = {'both': 0, 'left': 1, 'right': 1, 'none': 2}
        indexed_parts = list(enumerate(zip(parts, coalescence_types)))
        indexed_parts.sort(key=lambda x: priority_map.get(x[1][1], 3))

        from cartagen.utils.debug import plot_debug
        
        for part_idx, (part, coal_type) in indexed_parts:
            original_part = part
            
            if coal_type == 'none':
                # No coalescence: apply Gaussian smoothing (non-destructive)
                processed_part = smooth_gaussian(part, sigma_smooth, sample)
            
            elif coal_type in ['left', 'right']:
                # One-sided coalescence: use max_break (destructive)
                print('max_break')
                processed_part = max_break(part, width, exaggeration)
                
                # Propagate displacement to other parts if needed
                if not processed_part.equals(original_part):
                    print('prop')
                    processed_parts = propagation_network(processed_parts, part_idx, processed_part, propagation_distance, damping, tolerance)
            
            elif coal_type == 'both':
                # Two-sided coalescence: bend series
                print('accordion')
                
                # Try accordion first (destructive)
                accordion_result = accordion(part, width, exaggeration)
                
                # Check Hausdorff distance if threshold is specified
                if hausdorff_threshold is not None:
                    hd = hausdorff_distance(part, accordion_result)
                    
                    if hd > hausdorff_threshold:
                        # Accordion degrades planimetry too much
                        # Apply schematization first, then accordion again
                        print('schematization')
                        schematized = schematization(part)
                        processed_part = accordion(schematized, width, exaggeration)
                    else:
                        processed_part = accordion_result
                else:
                    processed_part = accordion_result

                # Propagate displacement after accordion
                if not processed_part.equals(original_part):
                    print('prop')
                    processed_parts = propagation_network(processed_parts, part_idx, processed_part, propagation_distance, damping, tolerance)
                
                # Re-detect coalescence on the processed bend series
                # to handle residual coalescence at bend summits
                sub_parts, sub_types = coalescence_splitting(processed_part, width, legibility_factor)
                
                # Process sub-parts if there's still coalescence
                sub_processed = list(sub_parts)  # Start with original sub-parts
                
                for sub_idx, (sub_part, sub_type) in enumerate(zip(sub_parts, sub_types)):
                    if sub_type in ['left', 'right']:
                        # Isolated bend in series: use min_break (destructive)
                        print('min_break')
                        original_sub_part = sub_part
                        sub_processed_part = min_break(sub_part, width * exaggeration)

                        # Propagate displacement to other sub-parts
                        if not sub_processed_part.equals(original_sub_part):
                            print('prop')
                            sub_processed = propagation_network(sub_processed, sub_idx, sub_processed_part, propagation_distance, damping, tolerance)
                
                processed_part = chain_line_parts(sub_processed)
            
            else:
                # Unknown type: keep original
                processed_part = part
            
            processed_parts[part_idx] = processed_part
        
        # Step 3: Chain processed parts back together
        chained_line = chain_line_parts(processed_parts)
        
        # Step 4: Apply light Platre smoothing to remove artifacts from chaining
        final_line = smooth_platre(chained_line, sigma=sigma_platre)
        
        # Step 5: Check for self-intersection (topology error)
        if not final_line.is_simple:
            print('topology_error')
            # Topology error detected: keep original line
            final_line = original_line
        
        # Store modification for final batch propagation
        if not final_line.equals(original_line):
            final_modifications.append((idx, final_line))
            network.at[idx, 'geometry'] = final_line
    
    # # Step 6: Final propagation of displacements through the entire network
    # if final_modifications and propagation_distance is not None:
    #     network = propagation_network_batch(
    #         network,
    #         final_modifications,
    #         propagation_distance=propagation_distance,
    #         damping=damping,
    #         tolerance=tolerance,
    #         inplace=True
    #     )
    
    return network