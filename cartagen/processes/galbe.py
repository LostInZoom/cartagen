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

def __chain_line(lines):
    if len(lines) == 0:
        return []

    if len(lines) == 1:
        return lines[0]

    line = []
    for i, l in enumerate(lines):
        if i == 0 or i == len(lines) - 1:
            line.extend(l.coords)
        else:
            line.extend(l.coords[1:-1])
    
    return LineString(line)

def galbe(
    network: gpd.GeoDataFrame,
    width: float,
    exaggeration: float = 1.0,
    legibility_factor: float = 1.7,
    sigma_smooth: float = 5,
    sigma_platre: float = 1,
    curvature_platre: float = 0.05,
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
                processed_parts[part_idx] = smooth_gaussian(part, sigma_smooth, sample)
            
            elif coal_type in ['left', 'right']:
                # One-sided coalescence: use max_break (destructive)
                processed_part = max_break(part, width, exaggeration)
                
                # Propagate displacement to other parts if needed
                if not processed_part.equals(original_part):
                    print('max_break')
                    newparts = propagation_network(processed_parts, part_idx, processed_part, propagation_distance, damping, tolerance, inplace=True)
                    processed_parts = newparts

            elif coal_type == 'both':
                # Two-sided coalescence: bend series
                print('accordion')

                # Try accordion first (destructive)
                accordion_result = accordion(part, width, exaggeration)
                
                # Check Hausdorff distance
                temp = processed_parts.copy()
                propagation_network(temp, part_idx, accordion_result, propagation_distance, damping, tolerance, inplace=True)
                hd = hausdorff_distance(line, __chain_line(temp))
                
                if hd > hausdorff_threshold:
                    # Accordion degrades planimetry too much
                    # Apply schematization first, then accordion again
                    print('schematization')
                    schematized = schematization(part)
                    processed_part = accordion(schematized, width, exaggeration)
                else:
                    processed_part = accordion_result
                
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
                            propagation_network(sub_processed, sub_idx, sub_processed_part, propagation_distance, damping, tolerance, inplace=True)
                
                processed_part = __chain_line(sub_processed)
                processed_parts[part_idx] = processed_part

        # Step 3: Chain processed parts back together
        chained_line = __chain_line(processed_parts)
        
        # Step 4: Apply light Platre smoothing to remove artifacts from chaining
        final_line = smooth_platre(chained_line, sigma=sigma_platre, curvature=curvature_platre)
        
        # # Step 5: Check for self-intersection (topology error)
        # if not final_line.is_simple:
        #     print('topology_error')
        #     # Topology error detected: keep original line
        #     final_line = line
        
        # Store modification for final batch propagation
        if not final_line.equals(line):
            print('validating')
            final_modifications.append((idx, final_line))
            network.at[idx, 'geometry'] = final_line
    
    return network