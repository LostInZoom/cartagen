from typing import List, Dict, Literal, Tuple
import numpy as np
import shapely
from shapely import Point, LineString
from shapely.ops import unary_union

from cartagen.utils.partitioning.network import network_faces

def dilate_line(line, offset, cap_style='round', quad_segs=8):
    """
    Dilate a line on both sides.

    This algorithm proposed by Mustière :footcite:p:`mustiere:2001-a` dilates
    a line on both sides by a given distance in meters. It is the basis of many
    mountain roads generalisation algorithm.
    
    Parameters
    ----------
    line : LineString
        The line to offset.
    offset : float
        The length of the offset to apply in meters.
    cap_style : str, optional
        The type of caps at the start and end of the provided line. Possible values are 'round' or 'flat'.
    quad_segs : int, optional
        The number of point allowed per circle quadrant when interpolating points using round method.

    Returns
    -------
    left, right: tuple of list of LineString
        A tuple of two list of LineString, the left and the right side of the dilation.

    See Also
    --------
    offset_line : 
        This function preserves the relationship between the provided list of points and the result.
    circle_interpolation :
        The function used to interpolate point along the quadrant of a circle.

    Warning
    -------
    This is basically a buffering algorithm rebuild from scratch. It is much slower than the
    :func:`buffer() <shapely.buffer>` algorithm proposed by shapely, so you better use this one.
    It exists only to be used alongside :func:`offset_line() <cartagen.offset_line>` to keep
    the relationship between the vertexes of the original line and the vertexes of the dilated one.

    Notes
    -----
    This algorithm returns dilation as lists of LineString. It often is composed of only one line but
    sometimes, it can creates multiple lines when the dilation creates holes.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(1, 1), (2, 3), (5, 2)])
    >>> dilate_line(line, 1.0)
    ([<LINESTRING (0.553 0.106, 0.715 0.042, 0.886 0.007, 1.06 0.002, 1.232 0.027,...>], [<LINESTRING (0.553 0.106, 0.404 0.197, 0.274 0.312, 0.165 0.449, 0.082 0.602...>])
    """

    # Check if provided cap style is allowed
    acaps = ['round', 'flat']
    if cap_style not in acaps:
        raise Exception('Offset curve cap style unknown: {0}.'.format(cap_style))

    # If offset is 0, return the provided line
    if offset == 0:
        return line

    # Calculate the offset points along the line
    groups1 = offset_line(line, offset, cap_style, quad_segs)
    groups2 = offset_line(line, -offset, cap_style, quad_segs)

    # Reconstruct the line into parts
    left, right = reconstruct_line(groups1, groups2, line, offset)

    return left[0], right[0]

def offset_line(
    line: LineString,
    offset: float,
    cap_style: Literal['round', 'flat'] = 'round',
    quad_segs: int = 8
) -> List[Dict]:
    """
    Offset the vertices of a line on one side.
 
    Offset the vertex of the line by a given distance and keeps the relationship
    between the line vertices and the result.
 
    Parameters
    ----------
    line : LineString
        The line to offset.
    offset : float
        The length of the offset to apply in meters. 
        Negative value for left-side dilation, positive for right-side.
    cap_style : {'round', 'flat'}, optional
        The type of caps at the start and end of the line.
        Default is 'round'.
    quad_segs : int, optional
        The number of points allowed per circle quadrant when interpolating 
        points using round method. Default is 8.
 
    Returns
    -------
    list of dict
        The index of each dict corresponds to the index of the provided points, 
        with the keys:
        
        - 'type': The type of point ('start', 'concave', 'convex', 'end')
        - 'original': The index of the vertex inside the input line
        - 'projected': The list of projected coordinates
 
    Raises
    ------
    ValueError
        If the line has fewer than 2 points, or if quad_segs < 1.
    
    See Also
    --------
    dilate_line : 
        This function loses the relationship between the provided line vertex 
        and the offset points.
    circle_interpolation :
        The function used to interpolate points along the quadrant of a circle.
    """
    # Validation des entrées
    if len(line.coords) < 2:
        raise ValueError("Line must have at least 2 points")
    
    if quad_segs < 1:
        raise ValueError("quad_segs must be at least 1")
    
    if abs(offset) < 1e-10:
        # Si offset est essentiellement zéro, retourner la ligne originale
        return [
            {
                'type': 'start' if i == 0 else ('end' if i == len(line.coords) - 1 else 'convex'),
                'original': i,
                'projected': [line.coords[i]]
            }
            for i in range(len(line.coords))
        ]
 
    def calculate_extremity(p1: Tuple[float, float], p2: Tuple[float, float]) -> Tuple[float, float]:
        """
        Calculate the rounded extremity of the offset curve.
        
        Parameters
        ----------
        p1 : tuple
            First point (x1, y1)
        p2 : tuple
            Second point (x2, y2)
            
        Returns
        -------
        tuple
            The extremity point (x, y)
        """
        x1, y1 = p1
        x2, y2 = p2
        
        # Distance entre les deux points
        distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        # Éviter la division par zéro
        if distance < 1e-10:
            return p2
        
        # Vecteur normalisé dans la direction p1 -> p2
        dx = (x2 - x1) / distance
        dy = (y2 - y1) / distance
        
        # Point situé à une distance 'offset' de p2 dans la direction p1 -> p2
        c_x = x2 + dx * abs(offset)
        c_y = y2 + dy * abs(offset)
        
        return (c_x, c_y)
 
    def project_point(p: Tuple[float, float], direction: np.ndarray, offset_dist: float) -> Tuple[float, float]:
        """
        Project a point perpendicular to a direction vector.
        
        Parameters
        ----------
        p : tuple
            The point to project
        direction : ndarray
            The normalized direction vector
        offset_dist : float
            The offset distance
            
        Returns
        -------
        tuple
            The projected point
        """
        # Vecteur orthogonal (rotation de 90°)
        orthogonal = np.array([-direction[1], direction[0]])
        
        # Point projeté
        projected = np.array(p) - offset_dist * orthogonal
        
        return (float(projected[0]), float(projected[1]))
 
    def snap_to_precision(point: Tuple[float, float], precision: float = 1e-9) -> Tuple[float, float]:
        """Snap coordinates to precision to avoid floating point issues."""
        return tuple(shapely.set_precision(Point(point), precision).coords[0])
 
    # Conversion des coordonnées en liste
    points = list(line.coords)
    n_points = len(points)
 
    # Stockage pour la ligne projetée
    pline = []
 
    # Boucle sur les sommets de la ligne originale
    for i in range(n_points):
        current_point = points[i]
 
        # Dernier point de la ligne
        if i >= n_points - 1:
            if cap_style == 'round':
                # Point précédent
                prev_point = points[i - 1]
                
                # Calcul de l'extrémité arrondie
                extremity = calculate_extremity(prev_point, current_point)
 
                # Interpolation entre le dernier point projeté et l'extrémité
                rotation = 'ccw' if offset < 0 else 'cw'
                
                last_projected = pline[-1]['projected'][-1]
                
                # AMÉLIORATION 4: Convertir en np.ndarray pour circle_interpolation
                try:
                    interpolation = circle_interpolation(
                        np.array(current_point, dtype=float), 
                        np.array(last_projected, dtype=float), 
                        np.array(extremity, dtype=float), 
                        rotation=rotation, 
                        quad_segs=quad_segs
                    )
                    
                    pline[-1]['type'] = 'end'
                    pline[-1]['projected'] = [snap_to_precision(tuple(p)) for p in interpolation]
                except ValueError as e:
                    # AMÉLIORATION 5: En cas d'erreur, fallback vers une simple ligne
                    print(f"Warning: circle_interpolation failed at end cap: {e}")
                    pline[-1]['type'] = 'end'
                    pline[-1]['projected'] = [
                        snap_to_precision(last_projected),
                        snap_to_precision(extremity)
                    ]
            
            break
 
        # Point suivant
        next_point = points[i + 1]
 
        # Vecteurs numpy
        v_current = np.array(current_point, dtype=float)
        v_next = np.array(next_point, dtype=float)
 
        # Vecteur de direction (du point suivant vers le point courant)
        direction = v_current - v_next
        direction_norm = np.linalg.norm(direction)
        
        # Éviter la division par zéro pour les points dupliqués
        if direction_norm < 1e-10:
            continue
        
        direction = direction / direction_norm
 
        # Projection des points
        projected_current = project_point(current_point, direction, offset)
        projected_next = project_point(next_point, direction, offset)
 
        # Premier point
        if i == 0:
            if cap_style == 'round':
                extremity = calculate_extremity(next_point, current_point)                
                rotation = 'ccw' if offset < 0 else 'cw'
                
                # AMÉLIORATION 6: Gestion d'erreur pour le premier cap
                try:
                    interpolation = circle_interpolation(
                        np.array(current_point, dtype=float), 
                        np.array(extremity, dtype=float),
                        np.array(projected_current, dtype=float), 
                        rotation=rotation, 
                        quad_segs=quad_segs
                    )
                    
                    projected_points = [snap_to_precision(tuple(p)) for p in interpolation]
                except ValueError as e:
                    print(f"Warning: circle_interpolation failed at start cap: {e}")
                    projected_points = [
                        snap_to_precision(extremity),
                        snap_to_precision(projected_current)
                    ]
            else:
                projected_points = [snap_to_precision(projected_current)]
 
            pline.append({
                'type': 'start',
                'original': i,
                'projected': projected_points
            })
 
        # Points intermédiaires
        else:
            # Récupération des derniers points projetés
            prev_projected_b = pline[-1]['projected'][-1]
            
            # Segments pour tester l'intersection
            seg_prev = LineString([pline[-2]['projected'][-1], prev_projected_b])
            seg_current = LineString([projected_current, projected_next])
 
            # Test de croisement (angle concave)
            if seg_prev.crosses(seg_current):
                point_type = 'concave'
                projected_points = [
                    snap_to_precision(prev_projected_b),
                    snap_to_precision(projected_current)
                ]
            else:
                # Angle convexe - interpolation circulaire
                point_type = 'convex'
                rotation = 'ccw' if offset < 0 else 'cw'
                
                # AMÉLIORATION 7: Gestion d'erreur pour les angles convexes
                try:
                    interpolation = circle_interpolation(
                        np.array(current_point, dtype=float), 
                        np.array(prev_projected_b, dtype=float), 
                        np.array(projected_current, dtype=float), 
                        rotation=rotation, 
                        quad_segs=quad_segs
                    )
                    projected_points = [snap_to_precision(tuple(p)) for p in interpolation]
                except ValueError as e:
                    print(f"Warning: circle_interpolation failed at vertex {i}: {e}")
                    # Fallback: simple ligne droite
                    projected_points = [
                        snap_to_precision(prev_projected_b),
                        snap_to_precision(projected_current)
                    ]
 
            # Mise à jour du dernier point
            pline[-1]['type'] = point_type
            pline[-1]['original'] = i
            pline[-1]['projected'] = projected_points
 
        # Ajout de la projection du point suivant
        pline.append({
            'type': 'end',
            'original': i + 1,
            'projected': [snap_to_precision(projected_next)]
        })
 
    return pline

def circle_interpolation(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    rotation: Literal['cw', 'ccw'] = 'cw',
    quad_segs: int = 8,
    tolerance: float = 1e-9
) -> np.ndarray:
    """
    Interpolate points along a circle arc.
    
    Given b and c, two points at equal distance from a third point a, 
    interpolates n points between b and c along the circle of center a 
    and radius ab. The number of provided points depends on the number 
    of segments allowed per circle quadrant.
    
    Parameters
    ----------
    a : ndarray, shape (2,)
        Center point coordinates
    b : ndarray, shape (2,)
        Start point coordinates
    c : ndarray, shape (2,)
        End point coordinates
    rotation : {'cw', 'ccw'}, optional
        Define the rotation direction:
        - 'cw': clockwise
        - 'ccw': counterclockwise (plane symmetry of cw using line ab)
        Default is 'cw'.
    quad_segs : int, optional
        The number of segments allowed on a quarter circle.
        This defines the resolution of interpolated points.
        Must be at least 1. Default is 8.
    tolerance : float, optional
        Numerical tolerance for distance and vector comparisons.
        Default is 1e-9.
    
    Returns
    -------
    ndarray, shape (n_points, 2)
        Array of interpolated points
    """
    # Validation
    if quad_segs < 1:
        raise ValueError(f"quad_segs must be at least 1, got {quad_segs}")
    
    if rotation not in ('cw', 'ccw'):
        raise ValueError(f"rotation must be 'cw' or 'ccw', got '{rotation}'")
    
    # Calcul des vecteurs
    ab = b - a
    ac = c - a
    
    # Calcul des rayons
    radius_ab = np.linalg.norm(ab)
    radius_ac = np.linalg.norm(ac)
    
    # Vérifications
    if radius_ab < tolerance:
        raise ValueError(f"Point b is too close to center a (distance: {radius_ab})")
    
    # AMÉLIORATION 1: Utiliser une tolérance relative pour la comparaison des rayons
    # au lieu d'une tolérance absolue
    radius_mean = (radius_ab + radius_ac) / 2
    relative_tolerance = max(tolerance, radius_mean * 1e-6)  # Au moins 1e-6 de la distance moyenne
    
    if abs(radius_ab - radius_ac) > relative_tolerance:
        # from cartagen.utils.debug import plot_debug
        # plot_debug(Point(a), Point(b), Point(c))

        raise ValueError(
            f"Points b and c are not equidistant from a: "
            f"distance(a,b)={radius_ab:.10f}, distance(a,c)={radius_ac:.10f}, "
            f"difference={abs(radius_ab - radius_ac):.10e}"
        )
    
    # AMÉLIORATION 2: Normaliser les vecteurs pour une meilleure stabilité numérique
    ab_normalized = ab / radius_ab
    ac_normalized = ac / radius_ac
    
    # Vérifier si les points sont colinéaires
    if np.allclose(ab_normalized, ac_normalized, atol=tolerance):
        return np.vstack([b, c])
    
    # Calcul de l'angle avec les vecteurs normalisés
    cos_angle = np.clip(np.dot(ab_normalized, ac_normalized), -1.0, 1.0)
    total_angle = np.arccos(cos_angle)
    
    if total_angle < tolerance:
        return np.vstack([b, c])
    
    # AMÉLIORATION 3: Utiliser le rayon moyen pour plus de stabilité
    radius = radius_mean
    
    # Nombre de points
    n_points = max(2, int(total_angle / (np.pi / 2) * quad_segs) + 2)
    
    if n_points == 2:
        return np.vstack([b, c])
    
    # Génération des angles
    angles = np.linspace(0, total_angle, n_points)
    if rotation == 'cw':
        angles = -angles
    
    # Calcul de toutes les matrices de rotation en une fois
    cos_angles = np.cos(angles[1:-1])
    sin_angles = np.sin(angles[1:-1])
    
    # Direction initiale normalisée
    direction = ab_normalized
    
    # Calcul vectorisé de tous les points intermédiaires
    dx, dy = direction
    rotated_x = cos_angles * dx - sin_angles * dy
    rotated_y = sin_angles * dx + cos_angles * dy
    
    # Points interpolés
    interpolated = a[:, np.newaxis] + radius * np.vstack([rotated_x, rotated_y])
    
    # Assemblage du résultat
    result = np.vstack([
        b,
        interpolated.T,
        c
    ])
    
    return result

def reconstruct_line(groups1, groups2, line, offset):
    """
    Reconstruct lines using provided groups of dilated points.

    This function creates buffer polygons from dilated line groups and extracts
    the boundary coordinates while preserving the relationship with original nodes.

    Parameters
    ----------
    groups1 : list of dict
        Groups of dilated points on the right side of the line.
        Each dict contains 'type', 'original', and 'projected' keys.
    groups2 : list of dict
        Groups of dilated points on the left side of the line.
    line : LineString
        The original line from which the dilation was calculated.
    offset : float
        The distance of the dilation in meters.

    Returns
    -------
    tuple of (tuple, tuple)
        ((lines_side2, breaks2), (lines_side1, breaks1))
        where lines are lists of LineString objects and breaks are lists of breakpoint dicts.
    """
    # Handle concavities by calculating intersections
    groups1 = _handle_concavity(groups1)
    groups2 = _handle_concavity(groups2)

    # Polygonize both sides to get boundary coordinates and holes
    polygon1, holes1 = _polygonize(line, groups1, offset)
    polygon2, holes2 = _polygonize(line, groups2, offset)

    # Build index maps for faster lookups
    projected1_set = _build_projected_set(groups1)
    projected2_set = _build_projected_set(groups2)
    
    iprojected1 = _build_indexed_projected(groups1)
    iprojected2 = _build_indexed_projected(groups2)

    # Extract the relevant part of polygon boundaries
    first1, last1 = _find_wrapping_indices(groups1, polygon1, iprojected1)
    last2, first2 = _find_wrapping_indices(groups2, polygon2, iprojected2)
    
    dilated1 = _extract_polygon_part(polygon1, first1, last1)
    dilated2 = _extract_polygon_part(polygon2, first2, last2)
    dilated2.reverse()

    # Handle crossing between the two sides
    dilated1, dilated2 = _handle_crossing(dilated1, dilated2)

    # Process main parts and holes
    parts1, breaks1 = _process_dilated_line(
        dilated1, groups1, iprojected1, projected1_set, holes1
    )
    parts2, breaks2 = _process_dilated_line(
        dilated2, groups2, iprojected2, projected2_set, holes2
    )
    
    return (parts2, breaks2), (parts1, breaks1)


def _handle_concavity(groups: List[Dict]) -> List[Dict]:
    """
    Handle concavities by calculating intersection points between segments.
    
    For concave vertices, find where the previous and following segments intersect
    and update the projected point to that intersection.
    
    Parameters
    ----------
    groups : list of dict
        List of projection groups.
    
    Returns
    -------
    list of dict
        Updated groups with concavities resolved.
    """
    updates = []
    
    for i, group in enumerate(groups):
        if group['type'] == 'concave':
            # Get segments before and after the concave vertex
            prev_pt = groups[i - 1]['projected'][-1]
            curr_start = group['projected'][0]
            curr_end = group['projected'][-1]
            next_pt = groups[i + 1]['projected'][0]
            
            seg1 = shapely.LineString([prev_pt, curr_start])
            seg2 = shapely.LineString([curr_end, next_pt])
            
            if shapely.crosses(seg1, seg2):
                intersection = shapely.intersection(seg1, seg2)
                updates.append((i, intersection.coords[0]))
    
    # Apply updates
    for idx, new_point in updates:
        groups[idx]['projected'] = [new_point]
    
    return groups


def _polygonize(line: LineString, groups: List[Dict], offset: float) -> Tuple[List, List]:
    """
    Create polygon from original line and projected groups, then extract boundary.
    
    Parameters
    ----------
    line : LineString
        Original line.
    groups : list of dict
        Projected point groups.
    offset : float
        Dilation offset.
    
    Returns
    -------
    tuple of (list, list)
        (boundary_coords, hole_coords_list)
    """
    coords = list(line.coords)
    reverse = [Point(c) for c in coords[::-1]]
    
    # Build vertex list: reversed original + all projected points
    vertices = reverse.copy()
    for group in groups:
        vertices.extend([Point(p) for p in group['projected']])
    vertices.append(reverse[0])
    
    # Create segments
    segments = [
        shapely.LineString([vertices[i], vertices[i + 1]])
        for i in range(len(vertices) - 1)
    ]
    
    # Polygonize the network
    polygons = network_faces(segments, convex_hull=False)
    
    # Extract holes (interior polygons beyond offset distance)
    holes = []
    abs_offset = abs(offset)
    for polygon in polygons:
        point_on_surface = polygon.point_on_surface()
        if shapely.distance(point_on_surface, line) >= abs_offset:
            coords_list = list(polygon.exterior.coords)
            # Reverse if offset is negative to maintain consistent orientation
            if offset > 0:
                holes.append(coords_list[::-1])
            else:
                holes.append(coords_list)
    
    # Merge all polygons
    merged_polygon = unary_union(polygons)
    
    # Handle MultiPolygon by keeping only the largest
    if merged_polygon.geom_type == 'MultiPolygon':
        merged_polygon = max(merged_polygon.geoms, key=lambda p: p.area)
    
    # Extract boundary (remove last duplicate vertex)
    boundary = list(merged_polygon.boundary.coords)[:-1]
    
    return boundary, holes


def _build_projected_set(groups: List[Dict]) -> set:
    """
    Build a set of all projected points for O(1) membership testing.
    
    Parameters
    ----------
    groups : list of dict
        Projected point groups.
    
    Returns
    -------
    set
        Set of all projected point coordinates.
    """
    return {pt for group in groups for pt in group['projected']}


def _build_indexed_projected(groups: List[Dict]) -> List[Tuple[int, int]]:
    """
    Build list of (group_index, node_index) tuples for all projected points.
    
    Parameters
    ----------
    groups : list of dict
        Projected point groups.
    
    Returns
    -------
    list of tuple
        List of (igroup, inode) index pairs.
    """
    return [
        (igroup, inode)
        for igroup, group in enumerate(groups)
        for inode in range(len(group['projected']))
    ]


def _find_wrapping_indices(
    groups: List[Dict],
    polygon: List,
    iprojected: List[Tuple[int, int]]
) -> Tuple[int, int]:
    """
    Find start and end indices in polygon boundary that match projected points.
    
    Parameters
    ----------
    groups : list of dict
        Projected point groups.
    polygon : list
        Polygon boundary coordinates.
    iprojected : list of tuple
        Indexed projected points.
    
    Returns
    -------
    tuple of int
        (first_index, last_index) in polygon boundary.
    """
    # Build a dict for O(1) lookup: coord -> polygon_index
    polygon_index_map = {coord: idx for idx, coord in enumerate(polygon)}
    
    start = groups[0]['projected'][0]
    end = groups[-1]['projected'][-1]
    
    # Try to find start point
    first = polygon_index_map.get(start)
    if first is None:
        # Search through indexed projections
        for igroup, inode in iprojected:
            coord = groups[igroup]['projected'][inode]
            first = polygon_index_map.get(coord)
            if first is not None:
                break
    
    # Try to find end point
    last = polygon_index_map.get(end)
    if last is None:
        # Search in reverse
        for igroup, inode in reversed(iprojected):
            coord = groups[igroup]['projected'][inode]
            last = polygon_index_map.get(coord)
            if last is not None:
                break
    
    return first, last


def _extract_polygon_part(polygon: List, first: int, last: int) -> List:
    """
    Extract a portion of polygon boundary between two indices.
    
    Handles wrapping around the polygon boundary.
    
    Parameters
    ----------
    polygon : list
        Polygon boundary coordinates.
    first : int
        Start index.
    last : int
        End index (inclusive).
    
    Returns
    -------
    list
        Extracted coordinates.
    """
    if first > last:
        return polygon[first:] + polygon[:last + 1]
    else:
        return polygon[first:last + 1]


def _handle_crossing(dilated1: List, dilated2: List) -> Tuple[List, List]:
    """
    Handle cases where the two dilated sides cross each other.
    
    Truncates both sides at their intersection point.
    
    Parameters
    ----------
    dilated1 : list
        First dilated line coordinates.
    dilated2 : list
        Second dilated line coordinates.
    
    Returns
    -------
    tuple of list
        (truncated_dilated1, truncated_dilated2)
    """
    line1 = LineString(dilated1)
    line2 = LineString(dilated2)
    
    if not shapely.crosses(line1, line2):
        return dilated1, dilated2
    
    # Find the crossing point
    for i1 in range(len(dilated1) - 1):
        seg1 = LineString([dilated1[i1], dilated1[i1 + 1]])
        
        if shapely.crosses(seg1, line2):
            for i2 in range(len(dilated2) - 1):
                seg2 = LineString([dilated2[i2], dilated2[i2 + 1]])
                
                if shapely.crosses(seg1, seg2):
                    intersection = shapely.intersection(seg1, seg2)
                    index1, index2 = i1 + 1, i2 + 1
                    
                    # Truncate both lines at intersection
                    dilated1 = [intersection.coords[0]] + dilated1[index1:]
                    dilated2 = [intersection.coords[0]] + dilated2[index2:]
                    
                    return dilated1, dilated2
    
    return dilated1, dilated2


def _process_dilated_line(
    dilated: List,
    groups: List[Dict],
    iprojected: List[Tuple[int, int]],
    projected_set: set,
    holes: List
) -> Tuple[List[LineString], List[Dict]]:
    """
    Process a dilated line to create final LineString parts and breakpoints.
    
    Parameters
    ----------
    dilated : list
        Dilated line coordinates.
    groups : list of dict
        Projected point groups.
    iprojected : list of tuple
        Indexed projected points.
    projected_set : set
        Set of all projected coordinates.
    holes : list
        List of hole coordinate lists.
    
    Returns
    -------
    tuple of (list, list)
        (parts, all_breaks) where parts are LineString objects and all_breaks are breakpoint dicts.
    """
    # Find breakpoints (points in dilated but not in projected)
    breakpoints = [i for i, p in enumerate(dilated) if p not in projected_set]
    
    # Process main part
    main_part, main_breaks = _treat_part(
        dilated, breakpoints, groups, iprojected, is_hole=False
    )
    
    parts = [LineString(main_part)]
    all_breaks = main_breaks
    
    # Process each hole
    for hole in holes:
        hole_coords = hole[:-1]  # Remove duplicate last point
        hole_breakpoints = [i for i, h in enumerate(hole_coords) if h not in projected_set]
        
        hole_part, hole_breaks = _treat_part(
            hole_coords, hole_breakpoints, groups, iprojected, is_hole=True
        )
        
        parts.append(LineString(hole_part))
        all_breaks.extend(hole_breaks)
    
    return parts, all_breaks


def _treat_part(
    coords: List,
    breaks: List[int],
    groups: List[Dict],
    iprojected: List[Tuple[int, int]],
    is_hole: bool
) -> Tuple[List, List[Dict]]:
    """
    Process a part (main or hole) to handle breakpoints.
    
    Parameters
    ----------
    coords : list
        Coordinates of the part.
    breaks : list of int
        Indices of breakpoints in coords.
    groups : list of dict
        Projected point groups.
    iprojected : list of tuple
        Indexed projected points.
    is_hole : bool
        Whether this is a hole (True) or main part (False).
    
    Returns
    -------
    tuple of (list, list)
        (processed_coords, breakpoint_dicts)
    """
    if len(breaks) == 0:
        return coords, []
    
    # Remove boundary breaks for non-holes
    if not is_hole:
        breaks = [b for b in breaks if b != 0 and b != len(coords) - 1]
        if len(breaks) == 0:
            return coords, []
    
    breakpoint_info = []
    
    # Build a map from coordinates to (igroup, inode) for faster lookup
    coord_to_index = {}
    for ip, (igroup, inode) in enumerate(iprojected):
        coord = groups[igroup]['projected'][inode]
        if coord not in coord_to_index:
            coord_to_index[coord] = (igroup, inode, ip)
    
    # Process each breakpoint
    min_ip = float('inf')
    min_coords_idx = 0
    
    for break_idx in breaks:
        break_point = Point(coords[break_idx])
        break_buffer = break_point.buffer(1e-8)
        
        # Find segments that intersect this breakpoint
        s1, s2 = None, None
        ip1 = None
        
        for ip in range(len(iprojected) - 1):
            igroup, inode = iprojected[ip]
            current = groups[igroup]['projected'][inode]
            
            # Determine next point
            if inode == len(groups[igroup]['projected']) - 1:
                igroup2, inode2 = igroup + 1, 0
                group_indices = [igroup, igroup + 1]
            else:
                igroup2, inode2 = igroup, inode + 1
                group_indices = [igroup]
            
            following = groups[igroup2]['projected'][inode2]
            segment_geom = shapely.LineString([current, following])
            
            if shapely.intersects(segment_geom, break_buffer):
                segment_info = {
                    'n1': {'igroup': igroup, 'inode': inode, 'geometry': current},
                    'n2': {'igroup': igroup2, 'inode': inode2, 'geometry': following},
                    'geometry': segment_geom
                }
                
                if s1 is None:
                    s1 = segment_info
                    o1 = group_indices
                    ip1 = ip
                elif s2 is None:
                    s2 = segment_info
                    o2 = group_indices
                    break  # Found both segments
        
        if s1 is not None and s2 is not None:
            # Track minimum for hole processing
            if ip1 < min_ip:
                min_ip = ip1
                min_coords_idx = break_idx
            
            # Calculate intersection point
            crosspoint = shapely.intersection(s1['geometry'], s2['geometry']).coords[0]
            
            breakpoint_info.append({
                's1': s1,
                's2': s2,
                'o1': o1,
                'o2': o2,
                'geometry': crosspoint,
                'coords_idx': break_idx
            })
    
    # Build final part
    if is_hole:
        # For holes, rotate to start at minimum breakpoint
        part = coords[min_coords_idx:] + coords[:min_coords_idx] + [coords[min_coords_idx]]
        
        # Update breakpoint coordinates in place
        for bp in breakpoint_info:
            coords[bp['coords_idx']] = bp['geometry']
        
        # Rotate breakpoint_info to match
        break_idx_in_list = next(
            i for i, bp in enumerate(breakpoint_info) if bp['coords_idx'] == min_coords_idx
        )
        breakpoint_info = (
            breakpoint_info[break_idx_in_list:] + breakpoint_info[:break_idx_in_list]
        )
    else:
        # For main parts, build incrementally
        part = []
        previous_node_idx = None
        
        # Build coord to index map for faster lookup
        coord_index_map = {coord: idx for idx, coord in enumerate(coords)}
        
        for bp in breakpoint_info:
            s1 = bp['s1']
            crosspoint = bp['geometry']
            
            # Add segment from previous to this breakpoint
            if s1['n1']['geometry'] in coord_index_map:
                start_idx = 0 if previous_node_idx is None else previous_node_idx
                end_idx = coord_index_map[s1['n1']['geometry']]
                part.extend(coords[start_idx:end_idx + 1])
            
            part.append(crosspoint)
            
            # Update for next iteration
            s2 = bp['s2']
            if s2['n2']['geometry'] in coord_index_map:
                previous_node_idx = coord_index_map[s2['n2']['geometry']]
        
        # Add remaining coordinates
        if previous_node_idx is not None:
            part.extend(coords[previous_node_idx:])
        else:
            part = coords
    
    return part, breakpoint_info