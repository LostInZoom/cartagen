"""
Burghardt-Cecconi building typification algorithm.

Python implementation based on the original Java code from the CartAGen project
(IGN France), adapted to work with GeoPandas and Shapely.

References
----------
Burghardt, D., & Cecconi, A. (2007). Mesh simplification for building typification.
International Journal of Geographical Information Science, 21(3), 283–298.
https://doi.org/10.1080/13658810600912323
"""

from __future__ import annotations

import math
import warnings
from typing import List, Optional

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import Delaunay
from shapely.affinity import rotate, scale
from shapely.geometry import MultiPolygon, Point, Polygon
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def typify_buildings_burghardt_cecconi(
    buildings: gpd.GeoDataFrame,
    initial_scale: int = 25_000,
    final_scale: int = 50_000,
    ratio: Optional[float] = None,
    road_network: Optional[gpd.GeoDataFrame] = None,
    attributes: Optional[List[str]] = None,
    min_area: Optional[float] = None,
    elimination_area: Optional[float] = None,
    min_separation: Optional[float] = None,
    density_weighting: bool = True,
) -> gpd.GeoDataFrame:
    """
    Typify buildings using the Burghardt–Cecconi Delaunay-based algorithm.

    This algorithm was proposed by Burghardt & Cecconi :footcite:p:`burghardt:2007`
    and progressively merges pairs of buildings by collapsing the shortest edge
    of a constrained Delaunay triangulation built on building centroids until the
    target number of buildings is reached and all remaining neighbours are
    separated by at least the minimum legible distance at the target scale.

    The algorithm follows four main steps:

    1. **Pre-generalisation** – buildings smaller than an elimination threshold
       are removed; buildings smaller than the minimum area at the target scale
       are enlarged via homothety.
    2. **Constrained triangulation** – a Delaunay triangulation is built on
       building centroids; edges that cross road axes or block boundaries are
       removed; remaining edge weights are set to a density-weighted distance.
    3. **Iterative typification** – the shortest edge is repeatedly collapsed:
       the two adjacent buildings are replaced by a new representative building
       whose position is the area-weighted barycentre of the merged group and
       whose shape is derived from the smallest surrounding rectangle of the
       largest constituent building.
    4. **Attribute transfer** – selected attributes are copied from the largest
       building in each merged group to the representative building.

    Parameters
    ----------
    buildings : GeoDataFrame of Polygon or MultiPolygon
        Input buildings to typify.  Each row must have a valid polygon geometry.
        The GeoDataFrame index is preserved where possible.
    initial_scale : int, optional
        Denominator of the source map scale (default: 25 000).
    final_scale : int, optional
        Denominator of the target map scale (default: 50 000).
    ratio : float, optional
        Explicit reduction ratio for the number of buildings (e.g. ``0.5``
        halves the count).  When supplied, ``initial_scale`` and
        ``final_scale`` are ignored for the target-count calculation.
        Must be in the interval ``(0, 1]``.
    road_network : GeoDataFrame of LineString, optional
        Road network used to remove triangulation edges that cross roads.
        When ``None``, no road-based edge filtering is applied.
    attributes : list of str, optional
        Column names to transfer from the largest building in each merged
        group to the output representative building.  Columns absent from
        *buildings* are silently ignored.
    min_area : float, optional
        Minimum building footprint area **in map ground units² (e.g. m²)**.
        Buildings below this threshold at the target scale are enlarged by
        homothety.  When ``None``, defaults to
        ``(0.6 mm)² × (final_scale / 1000)²``  (i.e. 0.36 mm² symbol area
        converted to ground units).
    elimination_area : float, optional
        Buildings whose area is strictly below this value are eliminated
        before typification starts.  When ``None``, defaults to half of
        *min_area*.
    min_separation : float, optional
        Minimum desired separation between neighbouring buildings in ground
        units at the target scale.  Typification continues until all
        remaining triangulation edges exceed this value, **in addition** to
        the target-count criterion.  When ``None``, defaults to
        ``0.4 mm × (final_scale / 1000)``.
    density_weighting : bool, optional
        When ``True`` (default), edge weights are multiplied by the square
        root of the local building density (ratio of total building area to
        block area), giving priority to merging buildings in dense clusters.

    Returns
    -------
    GeoDataFrame of Polygon
        Typified buildings.  The result contains:

        - ``geometry`` – representative polygon at the target scale;
        - all columns listed in *attributes* (values from the largest
          constituent building);
        - ``n_merged`` – number of original buildings merged into this
          representative.

    Notes
    -----
    * The function works in the coordinate reference system of *buildings*
      as-is; no reprojection is performed.  Make sure the CRS uses metric
      units (e.g. a projected CRS) so that area and distance thresholds are
      meaningful.
    * Buildings that could not be assigned to any triangulation node (e.g.
      isolated buildings after road filtering) are kept unchanged.
    * The Töpfer radical law used to compute the default target count is
      ``n_target = n_initial × (initial_scale / final_scale)``.

    References
    ----------
    .. footbibliography::
    """
    # ------------------------------------------------------------------
    # 0. Input validation & defaults
    # ------------------------------------------------------------------
    if buildings is None or len(buildings) == 0:
        warnings.warn("Empty buildings GeoDataFrame – returning as-is.")
        return buildings.copy()

    # work on a clean copy with integer positional index
    gdf = buildings.copy().reset_index(drop=True)
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()
    gdf.geometry = gdf.geometry.apply(
        lambda g: g.geoms[0] if isinstance(g, MultiPolygon) and len(g.geoms) == 1
        else (g.convex_hull if isinstance(g, MultiPolygon) else g)
    )

    # scale-derived defaults (ground units, assuming CRS in metres)
    scale_factor = final_scale / 1000.0          # mm → m at target scale
    if min_area is None:
        min_area = (0.6 * scale_factor) ** 2     # ≈ 0.36 mm² → m² on ground
    if elimination_area is None:
        elimination_area = min_area * 0.5
    if min_separation is None:
        min_separation = 0.4 * scale_factor      # 0.4 mm → m on ground

    # target number of buildings (Töpfer radical law, linear variant)
    if ratio is not None:
        if not (0 < ratio <= 1):
            raise ValueError("ratio must be in (0, 1].")
        target_n = max(1, int(round(len(gdf) * ratio)))
    else:
        target_n = max(1, int(round(len(gdf) * initial_scale / final_scale)))

    # ------------------------------------------------------------------
    # 1. Pre-generalisation
    # ------------------------------------------------------------------
    gdf = _pre_generalise(gdf, min_area, elimination_area)

    if len(gdf) == 0:
        warnings.warn("All buildings eliminated during pre-generalisation.")
        return gpd.GeoDataFrame(columns=buildings.columns, crs=buildings.crs)

    if len(gdf) <= target_n:
        # nothing to do
        return _build_output(gdf, attributes, buildings.crs)

    # ------------------------------------------------------------------
    # 2. Build constrained Delaunay triangulation graph
    # ------------------------------------------------------------------
    graph = _build_triangulation_graph(
        gdf, road_network, density_weighting
    )

    # ------------------------------------------------------------------
    # 3. Iterative typification
    # ------------------------------------------------------------------
    groups = _typify_iterative(
        gdf, graph, target_n, min_separation, density_weighting
    )

    # ------------------------------------------------------------------
    # 4. Create representative geometries
    # ------------------------------------------------------------------
    result_rows = []
    for group_indices in groups:
        subset = gdf.loc[list(group_indices)]
        new_geom = _create_representative_geometry(subset, min_area, scale_factor)
        row = {"geometry": new_geom, "n_merged": len(group_indices)}
        # attribute transfer from the largest building
        if attributes:
            largest_idx = subset.geometry.area.idxmax()
            for col in attributes:
                if col in subset.columns:
                    row[col] = subset.at[largest_idx, col]
        result_rows.append(row)

    result = gpd.GeoDataFrame(result_rows, crs=buildings.crs)
    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _pre_generalise(
    gdf: gpd.GeoDataFrame,
    min_area: float,
    elimination_area: float,
) -> gpd.GeoDataFrame:
    """Remove tiny buildings and enlarge below-threshold ones."""
    # eliminate
    mask_keep = gdf.geometry.area >= elimination_area
    gdf = gdf[mask_keep].copy()

    # enlarge
    def _enlarge(geom: Polygon, area: float) -> Polygon:
        if area < min_area:
            factor = math.sqrt(min_area / area)
            centroid = geom.centroid
            return scale(geom, xfact=factor, yfact=factor, origin=centroid)
        return geom

    gdf.geometry = [
        _enlarge(row.geometry, row.geometry.area)
        for _, row in gdf.iterrows()
    ]
    return gdf.reset_index(drop=True)


def _build_triangulation_graph(
    gdf: gpd.GeoDataFrame,
    road_network: Optional[gpd.GeoDataFrame],
    density_weighting: bool,
) -> nx.Graph:
    """
    Build a Delaunay triangulation graph on building centroids.

    Edges crossing roads are removed.  Edge weights are set to the
    density-weighted distance between the two building geometries.
    """
    centroids = np.array(
        [[geom.centroid.x, geom.centroid.y] for geom in gdf.geometry]
    )

    if len(centroids) < 3:
        # cannot triangulate fewer than 3 points – return a path graph
        G = nx.Graph()
        for i in range(len(gdf)):
            G.add_node(i)
        if len(gdf) == 2:
            d = gdf.geometry.iloc[0].distance(gdf.geometry.iloc[1])
            G.add_edge(0, 1, weight=d)
        return G

    tri = Delaunay(centroids)

    G = nx.Graph()
    for i in range(len(gdf)):
        G.add_node(i)

    # build union of road geometries for fast intersection test
    road_union = None
    if road_network is not None and len(road_network) > 0:
        road_union = unary_union(road_network.geometry.values)

    # precompute local building density per building
    # (approximated as the ratio of building area to its Voronoi-like neighbourhood)
    densities = _compute_local_densities(gdf, centroids)

    # add edges from triangulation
    seen = set()
    for simplex in tri.simplices:
        for a, b in [(0, 1), (1, 2), (0, 2)]:
            i, j = simplex[a], simplex[b]
            if i > j:
                i, j = j, i
            if (i, j) in seen:
                continue
            seen.add((i, j))

            # check road intersection
            if road_union is not None:
                edge_line = _make_line(centroids[i], centroids[j])
                if edge_line.intersects(road_union):
                    continue  # skip cross-road edges

            # compute weight
            dist = gdf.geometry.iloc[i].distance(gdf.geometry.iloc[j])
            if density_weighting:
                local_density = (densities[i] + densities[j]) / 2.0
                weight = math.sqrt(max(local_density, 1e-9)) * dist
            else:
                weight = dist

            G.add_edge(i, j, weight=weight)

    return G


def _compute_local_densities(
    gdf: gpd.GeoDataFrame,
    centroids: np.ndarray,
) -> np.ndarray:
    """
    Approximate local building density for each building.

    Density is estimated as the ratio of building area to the area of the
    convex hull of its k nearest neighbours (k=5).  Returns a value in [0, 1].
    """
    n = len(gdf)
    densities = np.ones(n)
    k = min(6, n)  # number of neighbours including self

    if n < 3:
        return densities

    from scipy.spatial import cKDTree
    tree = cKDTree(centroids)
    distances, indices = tree.query(centroids, k=k)

    for i in range(n):
        neigh_idx = indices[i]
        neigh_geoms = [gdf.geometry.iloc[idx] for idx in neigh_idx]
        total_building_area = sum(g.area for g in neigh_geoms)
        pts = centroids[neigh_idx]
        if len(pts) >= 3:
            try:
                hull_area = Polygon(pts).convex_hull.area
            except Exception:
                hull_area = 0.0
        else:
            hull_area = 0.0
        if hull_area > 0:
            densities[i] = min(total_building_area / hull_area, 1.0)
        else:
            densities[i] = 1.0

    return densities


def _typify_iterative(
    gdf: gpd.GeoDataFrame,
    graph: nx.Graph,
    target_n: int,
    min_separation: float,
    density_weighting: bool,
) -> List[set]:
    """
    Iteratively merge the two nodes linked by the shortest edge.

    Returns a list of sets, each set containing the original GeoDataFrame
    indices that were merged into one representative building.
    """
    # each node starts as its own group
    groups: dict[int, set] = {i: {i} for i in range(len(gdf))}
    # representative node index for each logical group
    node_to_group: dict[int, int] = {i: i for i in range(len(gdf))}

    G = graph.copy()

    while len(groups) > target_n:
        if G.number_of_edges() == 0:
            break

        # find the shortest (minimum-weight) edge
        min_edge = min(G.edges(data=True), key=lambda e: e[2].get("weight", 0))
        u, v, data = min_edge
        w = data.get("weight", 0.0)

        # stop if minimum separation is already satisfied and target reached
        if len(groups) <= target_n and w >= min_separation:
            break

        # merge groups u and v into u
        groups[u] = groups[u] | groups[v]
        del groups[v]
        node_to_group.pop(v, None)

        # reattach all edges of v to u
        for neigh in list(G.neighbors(v)):
            if neigh == u:
                continue
            # compute new weight as distance between new merged centroid and neighbour
            group_u_indices = list(groups[u])
            group_neigh_indices = list(groups.get(neigh, {neigh}))
            merged_geom = unary_union([gdf.geometry.iloc[i] for i in group_u_indices])
            neigh_geom = unary_union([gdf.geometry.iloc[i] for i in group_neigh_indices])
            new_dist = merged_geom.distance(neigh_geom)
            if density_weighting:
                # simplified: use mean area ratio as proxy for density
                density = min(
                    merged_geom.area / max(merged_geom.convex_hull.area, 1e-9), 1.0
                )
                new_w = math.sqrt(max(density, 1e-9)) * new_dist
            else:
                new_w = new_dist

            if G.has_edge(u, neigh):
                # keep the minimum weight
                old_w = G[u][neigh].get("weight", new_w)
                G[u][neigh]["weight"] = min(old_w, new_w)
            else:
                G.add_edge(u, neigh, weight=new_w)

        G.remove_node(v)

    return list(groups.values())


def _create_representative_geometry(
    subset: gpd.GeoDataFrame,
    min_area: float,
    scale_factor: float,
) -> Polygon:
    """
    Create a representative rectangular building for a merged group.

    Position: area-weighted barycentre of constituent centroids.
    Shape: smallest surrounding rectangle of the largest building,
           scaled so that the total area equals the average area of the group.
    Orientation: general orientation of the largest building.
    """
    areas = subset.geometry.area.values
    centroids_x = np.array([g.centroid.x for g in subset.geometry])
    centroids_y = np.array([g.centroid.y for g in subset.geometry])

    total_area = areas.sum()
    cx = (centroids_x * areas).sum() / total_area
    cy = (centroids_y * areas).sum() / total_area
    barycentre = Point(cx, cy)

    # average area of the group
    avg_area = total_area / len(subset)
    avg_area = max(avg_area, min_area)

    # largest building drives shape & orientation
    largest_geom = subset.geometry.iloc[areas.argmax()]

    # smallest surrounding rectangle
    ssr = largest_geom.minimum_rotated_rectangle

    # side lengths of the SSR
    coords = list(ssr.exterior.coords)
    side1 = Point(coords[0]).distance(Point(coords[1]))
    side2 = Point(coords[1]).distance(Point(coords[2]))
    long_side = max(side1, side2)
    short_side = min(side1, side2)

    ratio_shape = long_side / max(short_side, 1e-9)

    # new dimensions preserving the shape ratio
    new_long = math.sqrt(ratio_shape * avg_area)
    new_short = avg_area / new_long

    # build a rectangle centred at origin
    half_long = new_long / 2.0
    half_short = new_short / 2.0
    rect = Polygon([
        (-half_long, -half_short),
        ( half_long, -half_short),
        ( half_long,  half_short),
        (-half_long,  half_short),
    ])

    # orientation of the largest building (angle of its SSR's longest edge)
    dx = coords[1][0] - coords[0][0]
    dy = coords[1][1] - coords[0][1]
    if side1 >= side2:
        orientation_deg = math.degrees(math.atan2(dy, dx))
    else:
        dx2 = coords[2][0] - coords[1][0]
        dy2 = coords[2][1] - coords[1][1]
        orientation_deg = math.degrees(math.atan2(dy2, dx2))

    rect = rotate(rect, orientation_deg, origin=(0, 0), use_radians=False)

    # translate to barycentre
    rect = _translate_geom(rect, cx, cy)

    # ensure minimum area
    if rect.area < min_area:
        factor = math.sqrt(min_area / rect.area)
        rect = scale(rect, xfact=factor, yfact=factor, origin=barycentre)

    return rect


def _build_output(
    gdf: gpd.GeoDataFrame,
    attributes: Optional[List[str]],
    crs,
) -> gpd.GeoDataFrame:
    """Build output GeoDataFrame when no typification is needed."""
    cols = {"geometry": gdf.geometry.values, "n_merged": np.ones(len(gdf), dtype=int)}
    if attributes:
        for col in attributes:
            if col in gdf.columns:
                cols[col] = gdf[col].values
    return gpd.GeoDataFrame(cols, crs=crs)


def _make_line(p1: np.ndarray, p2: np.ndarray):
    """Return a Shapely LineString between two coordinate arrays."""
    from shapely.geometry import LineString
    return LineString([p1.tolist(), p2.tolist()])


def _translate_geom(geom: Polygon, dx: float, dy: float) -> Polygon:
    """Translate a Shapely geometry by (dx, dy)."""
    from shapely.affinity import translate
    return translate(geom, xoff=dx, yoff=dy)