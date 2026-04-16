import numpy as np
import math
import pandas as pd
import geopandas as gpd
import random

from shapely import affinity
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, LinearRing, Polygon, MultiPolygon, GeometryCollection

from cartagen.utils.lines.simplification.topographic import simplify_topographic

def generalise_area_patches(polygons, scaling, initial_scale=25000, final_scale=50000, level=1, closeness=None, max_reselect_iterations=5, area_tolerance=0.2):
    """
    Generalise area patches using enlargement and contraction. 

    This algorithm, proposed by Müller and Wang :footcite:p:`muller:1992`, simplifies
    area patches by enlarging or contracting them according to their
    relative size, then eliminating, reselecting, merging, displacing, and
    simplifying the survivors.

    This algorithm can be used for patches of vegetation for example, but
    also for islands, or any areal geographic entities.

    Parameters
    ----------
    polygons : GeoDataFrame of Polygon or MultiPolygon
        Input patches in projected coordinates (metres). MultiPolygon
        rows are supported and are reassembled after generalisation.
    scaling : float
        User-defined maximum blanket width in **mm on the source map**
        (= t* in eq. 2).  Internally converted to ground units as
        ``scaling × 1e-3 × initial_scale``.
    initial_scale : int
        Denominator of the source map scale (M_s).  Default 25 000.
    final_scale : int
        Denominator of the target map scale (M_t).  Default 50 000.
    level : int {1, 2, 3, 4}
        Radical-law selection model (I–IV in the paper). Lower values are
        more selective (fewer patches retained); higher values are less
        selective.  The reselect loop may adjust this automatically.
    closeness : float or None
        Distance threshold in ground units. Patches closer than this value
        are subject to spatial constraints (see step 6).  ``None`` disables
        the spatial-constraint step entirely.
    max_reselect_iterations : int
        Maximum number of reselect feedback iterations (default 5).
    area_tolerance : float
        Fractional tolerance on total patch area for the reselect loop
        (e.g. 0.05 = 5 %).

    Returns
    -------
    GeoDataFrame
        Generalised patches with Polygon or MultiPolygon geometries
        and an updated ``area`` column, in the same CRS as the input.

    References
    ----------
    .. footbibliography::
    """
    if level not in [1, 2, 3, 4]:
        raise ValueError("Level must be 1, 2, 3, or 4")

    # ------------------------------------------------------------------
    # 1. MultiPolygon explosion
    # ------------------------------------------------------------------
    polygons = polygons.copy().reset_index(drop=True)
    # Track which original row each exploded part came from
    polygons['_orig_idx'] = polygons.index
    has_multipolygons = polygons.geometry.geom_type.isin(['MultiPolygon']).any()
    if has_multipolygons:
        polygons = polygons.explode(index_parts=False).reset_index(drop=True)

    # ------------------------------------------------------------------
    # 2. Pre-process
    # ------------------------------------------------------------------
    # Convert t* from mm on source map to ground units
    scaling_ground = scaling * 1e-3 * initial_scale

    polygons['area'] = polygons.geometry.area
    polygons = polygons.sort_values(by='area', ascending=False).reset_index(drop=True)

    original_total_area = polygons['area'].sum()

    # Eliminate threshold (eq. 6): ω̄_s = ω̄_t × M_t²  (ground units)
    omega_t = 0.5e-6  # 0.5 mm² expressed in m²
    ground_threshold_source = omega_t * (final_scale ** 2)

    # ------------------------------------------------------------------
    # 3–9. Radical-law selection + reselect feedback loop
    # ------------------------------------------------------------------
    Ns = len(polygons)
    max_area = polygons['area'].iloc[0]
    scale_ratio = initial_scale / final_scale
    powers = {1: 2.0, 2: 1.5, 3: 1.0, 4: 0.5}

    current_level = level
    for _reselect_iter in range(max_reselect_iterations):

        # 3. Radical-law threshold T
        Nt = Ns * (scale_ratio ** powers[current_level])
        idx_T = int(np.clip(Nt, 0, Ns - 1))
        T = polygons['area'].iloc[idx_T]

        # 4. Compactness
        polys_work = polygons.copy()
        polys_work['compactness'] = (
            (4 * np.pi * polys_work['area']) / (polys_work.geometry.length ** 2)
        ).clip(upper=1.0)

        # 5. Blanket-width assignment (eqs. 2–3)
        MAX = max_area
        max_diff = MAX - T
        if max_diff <= 0:
            polys_work['blank_w'] = 0.0
        else:
            K = scaling_ground / math.sqrt(max_diff)
            polys_work['blank_w'] = (
                polys_work['compactness'] * K *
                np.sqrt(np.abs(polys_work['area'] - T))
            )
            polys_work.loc[polys_work['area'] < T, 'blank_w'] *= -1

        # 6. Spatial constraints
        if closeness is not None and closeness > 0:
            sindex = polys_work.sindex
            active_mask = polys_work['blank_w'] != 0
            to_check = polys_work[active_mask]

            if not to_check.empty:
                indices, distances = sindex.nearest(
                    to_check.geometry, exclusive=True, return_distance=True
                )
                min_dists = np.full(len(polys_work), np.inf)
                query_pos = indices[0]
                to_check_positions = np.where(active_mask)[0]
                for qp, dist in zip(query_pos, distances):
                    orig_pos = to_check_positions[qp]
                    min_dists[orig_pos] = min(min_dists[orig_pos], dist)

                close_mask = (polys_work['blank_w'] > 0) & (min_dists < closeness)
                polys_work.loc[close_mask, 'blank_w'] = 0.0

                remote_mask = (polys_work['blank_w'] < 0) & (min_dists > closeness)
                polys_work.loc[remote_mask, 'blank_w'] = 0.0

        # 7. Buffer application
        polys_work['geometry'] = polys_work.geometry.buffer(polys_work['blank_w'])

        # 8. Eliminate
        survivors = polys_work[
            (~polys_work.geometry.is_empty) &
            (polys_work.geometry.area >= ground_threshold_source)
        ].copy()

        eliminated = polys_work[
            polys_work.geometry.is_empty |
            (polys_work.geometry.area < ground_threshold_source)
        ].copy()

        # 9. Reselect
        reselected_geoms = gpd.GeoDataFrame(columns=eliminated.columns, crs=polygons.crs)

        if not eliminated.empty and not survivors.empty:
            eps_ground = 0.4e-3 * final_scale

            eliminated['dist_to_survivor'] = eliminated.geometry.apply(
                lambda x: survivors.distance(x).min()
            )

            isolated = eliminated[eliminated['dist_to_survivor'] > eps_ground].copy()
            cluster = eliminated[eliminated['dist_to_survivor'] <= eps_ground].copy()

            dfs_to_concat = []
            if not isolated.empty:
                dfs_to_concat.append(isolated)
            if not cluster.empty:
                reduction_factor = final_scale / initial_scale
                reselection_ratio = 1.0 / reduction_factor
                num_to_keep = max(1, int(len(cluster) * reselection_ratio))
                cluster_reselected = cluster.sample(n=min(num_to_keep, len(cluster)))
                if not cluster_reselected.empty:
                    dfs_to_concat.append(cluster_reselected)

            if dfs_to_concat:
                reselected_geoms = pd.concat(dfs_to_concat)

        if not reselected_geoms.empty:
            new_radius = math.sqrt(ground_threshold_source / math.pi)
            reselected_geoms = reselected_geoms.copy()
            reselected_geoms['geometry'] = (
                reselected_geoms.geometry.centroid.buffer(new_radius)
            )
            valid_dfs = [df for df in [survivors, reselected_geoms] if not df.empty]
            result_polys = pd.concat(valid_dfs)
        else:
            result_polys = survivors

        if result_polys.empty:
            break

        # Area feedback check
        new_total_area = result_polys.geometry.area.sum()
        relative_diff = abs(new_total_area - original_total_area) / max(original_total_area, 1e-12)

        if relative_diff <= area_tolerance:
            break

        if new_total_area < original_total_area:
            if current_level < 4:
                current_level += 1
            else:
                break
        else:
            if current_level > 1:
                current_level -= 1
            else:
                break

    polygons = gpd.GeoDataFrame(result_polys, crs=polygons.crs).reset_index(drop=True)
    polygons['area'] = polygons.geometry.area

    # ------------------------------------------------------------------
    # 10. Merge
    # ------------------------------------------------------------------
    merged_geometry = polygons.geometry.unary_union

    if merged_geometry is None or merged_geometry.is_empty:
        return gpd.GeoDataFrame(geometry=[], crs=polygons.crs)

    if merged_geometry.geom_type == 'MultiPolygon':
        distinct_patches = list(merged_geometry.geoms)
    else:
        distinct_patches = [merged_geometry]

    polygons = gpd.GeoDataFrame(
        geometry=distinct_patches, crs=polygons.crs
    ).reset_index(drop=True)
    polygons['area'] = polygons.geometry.area

    # ------------------------------------------------------------------
    # 11. Displace
    # ------------------------------------------------------------------
    eps_ground = 0.4e-3 * final_scale
    n_patches = len(polygons)
    displacements = np.zeros((n_patches, 2))
    sindex = polygons.sindex

    for i in range(n_patches):
        geom_i = polygons.geometry.iloc[i]
        area_i = polygons['area'].iloc[i]
        centroid_i = geom_i.centroid

        potential_neighbors = sindex.query(geom_i.buffer(eps_ground))

        for j in potential_neighbors:
            if i >= j:
                continue

            geom_j = polygons.geometry.iloc[j]
            dist = geom_i.distance(geom_j)

            if 0 < dist < eps_ground:
                area_j = polygons['area'].iloc[j]
                centroid_j = geom_j.centroid

                dx = centroid_i.x - centroid_j.x
                dy = centroid_i.y - centroid_j.y
                vec_len = math.sqrt(dx ** 2 + dy ** 2)

                if vec_len == 0:
                    continue

                ux, uy = dx / vec_len, dy / vec_len
                gap_to_fill = eps_ground - dist
                total_area = area_i + area_j

                move_i = gap_to_fill * (area_j / total_area)
                move_j = gap_to_fill * (area_i / total_area)

                displacements[i] += [ux * move_i, uy * move_i]
                displacements[j] -= [ux * move_j, uy * move_j]

    new_geoms = []
    for i in range(n_patches):
        dx, dy = displacements[i]
        geom = polygons.geometry.iloc[i]
        if dx != 0 or dy != 0:
            new_geoms.append(affinity.translate(geom, xoff=dx, yoff=dy))
        else:
            new_geoms.append(geom)

    polygons['geometry'] = new_geoms

    # ------------------------------------------------------------------
    # 12. Topological validity check
    # ------------------------------------------------------------------
    def ensure_valid(geom):
        if not geom.is_valid:
            return geom.buffer(0)
        return geom

    polygons['geometry'] = polygons.geometry.apply(ensure_valid)
    polygons = polygons[~polygons.geometry.is_empty].copy()

    # ------------------------------------------------------------------
    # 13. Simplify contours
    # ------------------------------------------------------------------
    polygons['geometry'] = polygons.geometry.apply(simplify_topographic)
    polygons['geometry'] = polygons.geometry.apply(ensure_valid)
    polygons = polygons[~polygons.geometry.is_empty].copy()

    # ------------------------------------------------------------------
    # 14. MultiPolygon reassembly
    # ------------------------------------------------------------------
    if has_multipolygons and '_orig_idx' in polygons.columns:
        # Spatial join back to original footprints to recover _orig_idx
        original_footprints = gpd.GeoDataFrame(
            {'_orig_idx': polygons['_orig_idx'].unique()},
            geometry=[polygons.loc[polygons['_orig_idx'] == idx].geometry.unary_union for idx in polygons['_orig_idx'].unique()],
            crs=polygons.crs
        )
        # Dissolve parts sharing the same original index
        polygons = (
            polygons
            .dissolve(by='_orig_idx', as_index=False)
            .reset_index(drop=True)
        )
        polygons['area'] = polygons.geometry.area

    polygons = polygons.drop(columns=['_orig_idx'], errors='ignore')

    return polygons.reset_index(drop=True)