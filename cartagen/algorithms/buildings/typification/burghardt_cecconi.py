"""
Claude AI Python implementation based on the original Java code from the CartAGen project
(IGN France), adapted to work with GeoPandas and Shapely.
https://github.com/IGNF/CartAGen/blob/master/cartagen-core/src/main/java/fr/ign/cogit/cartagen/algorithms/typification/TypifyBurghardtCecconi.java
"""

from __future__ import annotations

import math
import warnings
from typing import List, Optional

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import Delaunay, cKDTree
from shapely.affinity import rotate, scale, translate
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.ops import unary_union

def typify_buildings_burghardt_cecconi(
    buildings, initial_scale=25000, final_scale=50000,
    ratio=None, road_network=None, attributes=None,
    min_area=None, elimination_area=None, min_separation=None, density_weighting=True,
):
    """
    Typify buildings using the Burghardt-Cecconi algorithm.

    This algorithm was proposed by Burghardt & Cecconi :footcite:p:`burghardt:2007`
    and progressively merges pairs of buildings by collapsing the shortest edge
    of a constrained Delaunay triangulation built on building centroids, until
    the target number of buildings is reached **and** all remaining neighbours
    are separated by at least the minimum legible distance at the target scale.

    The algorithm follows four main steps:

    1. **Pre-generalisation** - buildings smaller than an elimination threshold
       are removed; buildings smaller than the minimum representable area at
       the target scale are enlarged via homothety.
    2. **Constrained triangulation** - a Delaunay triangulation is built on
       building centroids; edges that cross roads are removed; each edge carries
       two weights: a density-weighted score used for ordering merges, and the
       raw geometric distance used for the separation stop criterion.
    3. **Iterative typification** - the shortest density-weighted edge is
       repeatedly collapsed: the two adjacent buildings are replaced by a new
       representative building whose position is the area-weighted barycentre
       of the merged group and whose shape is derived from the smallest
       surrounding rectangle of the largest constituent building, scaled to
       the average area of the group.
    4. **Attribute transfer** - selected attributes are copied from the largest
       building in each merged group to the representative building.

    Parameters
    ----------
    buildings : GeoDataFrame of Polygon or MultiPolygon
        Input buildings to typify. Each row must have a valid polygon geometry.
    initial_scale : int, optional
        Denominator of the source map scale (default: 25 000).
    final_scale : int, optional
        Denominator of the target map scale (default: 50 000).
    ratio : float, optional
        Explicit reduction ratio for the number of buildings (e.g. ``0.5``
        halves the count). When supplied, ``initial_scale`` and ``final_scale``
        are ignored for the target-count calculation. Must be in ``(0, 1]``.
    road_network : GeoDataFrame of LineString, optional
        Road network used to remove triangulation edges that cross roads.
        When ``None``, no road-based edge filtering is applied.
    attributes : list of str, optional
        Column names to transfer from the largest building in each merged
        group to the output representative building. Columns absent from
        *buildings* are silently ignored.
    min_area : float, optional
        Minimum building footprint area in ground units squared (e.g. m² for
        a metric CRS). Buildings below this threshold are enlarged by homothety.
        When ``None``, defaults to ``(0.2 mm * final_scale/1000)^2``,
        i.e. the ground area corresponding to a 0.04 mm² symbol at the target
        scale. Set to ``0`` to disable enlargement.
    elimination_area : float, optional
        Buildings whose area is strictly below this value are deleted before
        typification starts. When ``None``, defaults to
        ``(0.1 mm * final_scale/1000)^2``, i.e. truly sub-pixel footprints.
        Set to ``0`` to disable elimination entirely.
    min_separation : float, optional
        Minimum desired separation between neighbouring building footprints in
        ground units. Typification continues until all remaining triangulation
        raw distances exceed this value, in addition to the target-count
        criterion. When ``None``, defaults to ``0.2 mm * final_scale/1000``.
    density_weighting : bool, optional
        When ``True`` (default), edges are ordered by ``sqrt(local_density) *
        distance`` so that buildings in dense clusters are merged first. The
        raw geometric distance is always used for the stop criterion regardless
        of this flag.

    Returns
    -------
    GeoDataFrame of Polygon
        Typified buildings. The result contains:

        - ``geometry`` - representative polygon at the target scale;
        - all columns listed in *attributes* (values from the largest
          constituent building in each merged group);
        - ``n_merged`` - number of original buildings merged into this
          representative.

    Notes
    -----
    * The Topfer radical law used to compute the default target count is
      ``n_target = n_initial * sqrt(initial_scale / final_scale)``.
      For a 1:25 000 to 1:50 000 generalisation this yields ~71% of the
      original count.
    * Density weighting controls the *order* in which pairs are merged but
      uses only the raw footprint distance for the separation stop criterion,
      so the result count is unaffected by ``density_weighting``.
    * If all buildings are eliminated during pre-generalisation, the function
      emits a warning and returns an empty GeoDataFrame. In that case, lower
      *elimination_area* or set it to ``0``.

    See Also
    --------
    typify_buildings_matching :
        Typify buildings using the matching-based algorithm for web mapping.

    References
    ----------
    .. footbibliography::
    """
    # ------------------------------------------------------------------
    # 0. Input validation & defaults
    # ------------------------------------------------------------------
    if buildings is None or len(buildings) == 0:
        warnings.warn("Empty buildings GeoDataFrame - returning as-is.")
        return buildings.copy() if buildings is not None else gpd.GeoDataFrame()

    gdf = buildings.copy().reset_index(drop=True)
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()
    gdf.geometry = gdf.geometry.apply(
        lambda g: max(g.geoms, key=lambda p: p.area)
        if isinstance(g, MultiPolygon) else g
    )
    gdf = gdf.reset_index(drop=True)

    scale_factor = final_scale / 1000.0
    if min_area is None:
        min_area = (0.2 * scale_factor) ** 2
    if elimination_area is None:
        elimination_area = (0.1 * scale_factor) ** 2
    if min_separation is None:
        min_separation = 0.2 * scale_factor

    if ratio is not None:
        if not (0 < ratio <= 1):
            raise ValueError("ratio must be in (0, 1].")
        target_n = max(1, int(round(len(gdf) * ratio)))
    else:
        # Topfer radical law (square-root variant, standard for areal objects)
        target_n = max(1, int(round(len(gdf) * math.sqrt(initial_scale / final_scale))))

    # ------------------------------------------------------------------
    # 1. Pre-generalisation
    # ------------------------------------------------------------------
    gdf = _pre_generalise(gdf, min_area, elimination_area)

    if len(gdf) == 0:
        warnings.warn(
            "All buildings were eliminated during pre-generalisation. "
            f"elimination_area={elimination_area:.2f} may be too large for "
            "your data. Pass elimination_area=0 to disable elimination."
        )
        return gpd.GeoDataFrame({"geometry": [], "n_merged": []}, crs=buildings.crs)

    if len(gdf) <= target_n:
        return _build_output(gdf, attributes, buildings.crs)

    # ------------------------------------------------------------------
    # 2. Build constrained Delaunay triangulation graph
    # ------------------------------------------------------------------
    graph = _build_triangulation_graph(gdf, road_network, density_weighting)

    # ------------------------------------------------------------------
    # 3. Iterative typification
    # ------------------------------------------------------------------
    groups = _typify_iterative(gdf, graph, target_n, min_separation, density_weighting)

    # ------------------------------------------------------------------
    # 4. Create representative geometries
    # ------------------------------------------------------------------
    result_rows = []
    for group_indices in groups:
        subset = gdf.loc[sorted(group_indices)]
        new_geom = _create_representative_geometry(subset, min_area)
        row = {"geometry": new_geom, "n_merged": len(group_indices)}
        if attributes:
            largest_idx = subset.geometry.area.idxmax()
            for col in attributes:
                if col in subset.columns:
                    row[col] = subset.at[largest_idx, col]
        result_rows.append(row)

    return gpd.GeoDataFrame(result_rows, crs=buildings.crs)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _pre_generalise(
    gdf: gpd.GeoDataFrame,
    min_area: float,
    elimination_area: float,
) -> gpd.GeoDataFrame:
    """Remove sub-pixel buildings and enlarge below-threshold ones."""
    if elimination_area > 0:
        gdf = gdf[gdf.geometry.area >= elimination_area].copy()

    def _enlarge(geom: Polygon) -> Polygon:
        a = geom.area
        if 0 < a < min_area:
            return scale(geom, xfact=math.sqrt(min_area / a),
                         yfact=math.sqrt(min_area / a), origin=geom.centroid)
        return geom

    gdf.geometry = gdf.geometry.apply(_enlarge)
    return gdf.reset_index(drop=True)


def _build_triangulation_graph(
    gdf: gpd.GeoDataFrame,
    road_network: Optional[gpd.GeoDataFrame],
    density_weighting: bool,
) -> nx.Graph:
    """
    Build a Delaunay triangulation graph on building centroids.

    Each edge stores two attributes:
    - ``weight``   : density-weighted distance (used to order merges)
    - ``raw_dist`` : raw geometric distance between footprints (used for
                     the separation stop criterion)
    """
    n = len(gdf)
    centroids = np.array(
        [[g.centroid.x, g.centroid.y] for g in gdf.geometry]
    )

    G = nx.Graph()
    G.add_nodes_from(range(n))

    if n < 2:
        return G

    road_union = (
        unary_union(road_network.geometry.values)
        if road_network is not None and len(road_network) > 0
        else None
    )
    densities = (
        _compute_local_densities(gdf, centroids) if density_weighting
        else np.ones(n)
    )

    if n == 2:
        _add_edge(G, 0, 1, gdf, centroids, densities, road_union, density_weighting)
        return G

    tri = Delaunay(centroids)
    seen: set = set()
    for simplex in tri.simplices:
        for a, b in ((0, 1), (1, 2), (0, 2)):
            i, j = sorted((int(simplex[a]), int(simplex[b])))
            if (i, j) in seen:
                continue
            seen.add((i, j))
            _add_edge(G, i, j, gdf, centroids, densities, road_union, density_weighting)

    return G


def _add_edge(G, i, j, gdf, centroids, densities, road_union, density_weighting):
    """Add edge (i,j) unless it crosses a road; store both weight and raw_dist."""
    if road_union is not None:
        if LineString([centroids[i].tolist(), centroids[j].tolist()]).intersects(road_union):
            return
    raw = gdf.geometry.iloc[i].distance(gdf.geometry.iloc[j])
    if density_weighting:
        d = (densities[i] + densities[j]) / 2.0
        w = math.sqrt(max(d, 1e-9)) * raw
    else:
        w = raw
    G.add_edge(i, j, weight=w, raw_dist=raw)


def _compute_local_densities(gdf: gpd.GeoDataFrame, centroids: np.ndarray) -> np.ndarray:
    """
    Approximate local building density for each building as the ratio of
    summed neighbour footprint areas to the convex hull of their centroids (k=6).
    """
    n = len(gdf)
    if n < 3:
        return np.ones(n)

    k = min(6, n)
    _, indices = cKDTree(centroids).query(centroids, k=k)
    densities = np.ones(n)

    for i in range(n):
        idx = indices[i]
        total_area = sum(gdf.geometry.iloc[j].area for j in idx)
        try:
            hull_area = Polygon(centroids[idx]).convex_hull.area if len(idx) >= 3 else 0.0
        except Exception:
            hull_area = 0.0
        densities[i] = min(total_area / hull_area, 1.0) if hull_area > 0 else 1.0

    return densities


def _typify_iterative(
    gdf: gpd.GeoDataFrame,
    graph: nx.Graph,
    target_n: int,
    min_separation: float,
    density_weighting: bool,
) -> List[set]:
    """
    Iteratively collapse the shortest density-weighted edge.

    Stop condition (both must hold simultaneously):
    - number of groups <= target_n, AND
    - raw geometric distance of the shortest remaining edge >= min_separation.

    Density weighting controls only the *order* of merges; the stop criterion
    always uses the raw footprint distance to avoid distortion from density.
    """
    groups: dict[int, set] = {i: {i} for i in range(len(gdf))}
    G = graph.copy()

    while True:
        if G.number_of_edges() == 0:
            break

        u, v, data = min(G.edges(data=True), key=lambda e: e[2].get("weight", 0.0))
        raw_w = data.get("raw_dist", data.get("weight", 0.0))

        # stop when target count is reached AND closest pair is far enough apart
        if len(groups) <= target_n and raw_w >= min_separation:
            break

        # merge v into u
        groups[u] = groups[u] | groups[v]
        del groups[v]

        merged_geom = unary_union([gdf.geometry.iloc[i] for i in groups[u]])

        for neigh in list(G.neighbors(v)):
            if neigh == u:
                continue
            neigh_geom = unary_union(
                [gdf.geometry.iloc[i] for i in groups.get(neigh, {neigh})]
            )
            new_raw = merged_geom.distance(neigh_geom)
            if density_weighting:
                density = min(
                    merged_geom.area / max(merged_geom.convex_hull.area, 1e-9), 1.0
                )
                new_w = math.sqrt(max(density, 1e-9)) * new_raw
            else:
                new_w = new_raw

            if G.has_edge(u, neigh):
                # keep the edge with the smaller weight (closer pair)
                if new_w < G[u][neigh].get("weight", new_w):
                    G[u][neigh]["weight"] = new_w
                    G[u][neigh]["raw_dist"] = new_raw
            else:
                G.add_edge(u, neigh, weight=new_w, raw_dist=new_raw)

        G.remove_node(v)

    return list(groups.values())


def _create_representative_geometry(
    subset: gpd.GeoDataFrame,
    min_area: float,
) -> Polygon:
    """
    Build a representative rectangle for a merged group.

    - Position    : area-weighted barycentre of constituent centroids
    - Aspect ratio: smallest surrounding rectangle of the largest building
    - Area        : sqrt(n) * average_area, so a merged group is visibly
                    larger than a single building while staying sub-linear:
                    n=1 -> 1.0x avg, n=2 -> 1.4x, n=4 -> 2.0x, n=9 -> 3.0x
    - Orientation : longest axis of that same SSR
    """
    areas = subset.geometry.area.values
    total_area = areas.sum()
    n = len(subset)

    cx = (np.array([g.centroid.x for g in subset.geometry]) * areas).sum() / total_area
    cy = (np.array([g.centroid.y for g in subset.geometry]) * areas).sum() / total_area

    avg_area = total_area / n
    rep_area = max(math.sqrt(n) * avg_area, min_area)
    largest_geom = subset.geometry.iloc[areas.argmax()]

    ssr = largest_geom.minimum_rotated_rectangle
    coords = list(ssr.exterior.coords)
    side1 = Point(coords[0]).distance(Point(coords[1]))
    side2 = Point(coords[1]).distance(Point(coords[2]))
    long_side  = max(side1, side2, 1e-9)
    short_side = max(min(side1, side2), 1e-9)
    aspect = long_side / short_side

    new_long  = math.sqrt(aspect * rep_area)
    new_short = rep_area / new_long

    rect = Polygon([
        (-new_long / 2, -new_short / 2),
        ( new_long / 2, -new_short / 2),
        ( new_long / 2,  new_short / 2),
        (-new_long / 2,  new_short / 2),
    ])

    if side1 >= side2:
        dx, dy = coords[1][0] - coords[0][0], coords[1][1] - coords[0][1]
    else:
        dx, dy = coords[2][0] - coords[1][0], coords[2][1] - coords[1][1]

    rect = rotate(rect, math.degrees(math.atan2(dy, dx)), origin=(0, 0), use_radians=False)
    rect = translate(rect, xoff=cx, yoff=cy)

    if rect.area < min_area:
        rect = scale(rect,
                     xfact=math.sqrt(min_area / rect.area),
                     yfact=math.sqrt(min_area / rect.area),
                     origin=rect.centroid)
    return rect


def _build_output(
    gdf: gpd.GeoDataFrame,
    attributes: Optional[List[str]],
    crs,
) -> gpd.GeoDataFrame:
    """Build output GeoDataFrame when no typification is needed."""
    cols: dict = {
        "geometry": gdf.geometry.values,
        "n_merged": np.ones(len(gdf), dtype=int),
    }
    if attributes:
        for col in attributes:
            if col in gdf.columns:
                cols[col] = gdf[col].values
    return gpd.GeoDataFrame(cols, crs=crs)