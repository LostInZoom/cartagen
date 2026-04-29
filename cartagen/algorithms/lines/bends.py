import numpy as np
import shapely
from shapely import LineString, Point
from shapely import affinity, is_ccw, intersection

from cartagen.utils.math.vector import Vector2D
from cartagen.utils.geometry.bends import BendSeries
from cartagen.utils.geometry.line import extend_line_by_length, get_line_middle_point, merge_linestrings
from cartagen.utils.geometry.segment import Segment

def accordion(line, width, exaggeration=1.0, sigma=15.0, sample=None):
    """
    Stretch a series of bends to enlarge each bend.

    This algorithm was proposed by Plazanet. :footcite:p:`plazanet:1996`
    It is dedicated to the caricature of sinuous bend series. Like the musical instrument,
    the Accordion algorithm stretches the road to enlarge each bend of the series.

    The algorithm is part of the toolbox to generalise mountain roads that contain sinuous bend series
    and is rather used when there is room to enlarge the bend series.
    When the diffusion of the enlargement to the connected roads causes more legibility problems than the
    ones solved by accordion, schematization should be preferred.

    Parameters
    ----------
    line : LineString
        The line to apply the accordion algorithm.
    width : float
        The width in meters of the casing of the symbol.
    exaggeration : float, optional
        A multiplicator for the width parameter.
    sigma : float, optional
        Gaussian smoothing strength.
    sample : int, optional
        Gaussian smoothing sample size.

    Returns
    -------
    LineString

    See Also
    --------
    schematization :
        Remove bends from series of bends to simplify it.
    smooth_gaussian :
        Smooth a line and attenuate its inflexion points.

    References
    ----------
    .. footbibliography::
    """
    # Gestion des cas limites
    if line is None or line.is_empty or line.length < (width * exaggeration):
        return line

    # Detect individual bends inside the smoothed line
    bs = BendSeries(line, sigma, sample)

    # from cartagen.utils.debug import plot_debug
    # plot_debug(*[b.get_geom() for b in bs.get_bends()])

    # Storage for the future distorted individual lines
    distorted = []
    
    for b in bs.bends:
        # Get the translation vector
        v = __get_vector(b, width * exaggeration)

        # If the vector could not be found, keep the original bend
        if v is None:
            distorted.append(b.bend)
            continue

        # Retrieve bend coordinates
        coords = list(b.bend.coords)
        
        # Vérifier qu'il y a assez de points
        if len(coords) < 2:
            distorted.append(b.bend)
            continue
        
        # Précalculer les Points
        points = [Point(c) for c in coords]
        start, end = points[0], points[-1]

        # Calculate the vector formed by the start and end point of the bend
        vab = Vector2D.from_points(start, end)

        # Calculate the scalar product between the translation vector
        # and the vector formed by the start and end point of the bend
        scalar = v.scalar_product(vab)

        inverted = False

        # If the scalar product is negative, inverse the translation vector
        if scalar < 0:
            v = v.opposite()

        # Calculate the summit
        summit = b.get_bend_summit()

        # Calculate the vector formed by the start and the summit of the bend
        vas = Vector2D.from_points(start, summit)

        # Calculate the vector product of vas and vab
        product = vas.product(vab)

        # Check that the closed bend is clockwise
        clockwise = not shapely.is_ccw(LineString(coords + [coords[0]]))

        # Inverse the vector by comparing product and clockwise flag
        if (clockwise and product > 0) or (not clockwise and product < 0):
            v = v.opposite()

        # Calculer la longueur totale et les distances entre points
        length = b.bend.length
        
        # Calculer les distances cumulées depuis la FIN (comme l'original)
        cumulative_from_end = [0.0]
        for i in range(len(points) - 1, 0, -1):
            cumulative_from_end.insert(0, cumulative_from_end[0] + points[i].distance(points[i - 1]))

        # Translation avec facteur progressif (décroissant depuis la fin)
        translated = []
        for i, p in enumerate(points):
            factor = (length - cumulative_from_end[i]) / length
            vdist = v.copy()
            vdist.scalar_multiplication(factor)
            translated.append(vdist.translate(p))
        
        # Add the list of translated points to the full list of the whole bend serie
        distorted.append(LineString(translated))

    # Fusion optimisée des lignes
    if not distorted:
        return line
    
    all_coords = list(distorted[0].coords)
    for dline in distorted[1:]:
        # Calculer l'offset une seule fois
        offset = (all_coords[-1][0] - dline.coords[0][0], all_coords[-1][1] - dline.coords[0][1])
        # Ajouter directement les coordonnées décalées (sans le premier point)
        all_coords.extend((x + offset[0], y + offset[1]) for x, y in dline.coords[1:])

    displaced = LineString(all_coords)

    # Get the middle of the provided bend serie
    bm = get_line_middle_point(line)
    # Get the middle of the distorted line
    dm = get_line_middle_point(displaced)

    # Translate the distorted line by the difference between both middle points
    return shapely.affinity.translate(displaced, bm.x - dm.x, bm.y - dm.y)


def __get_vector(bend, width):
    """
    Calculate the translation vector for a bend.
    
    Parameters
    ----------
    bend : Bend
        The bend object.
    width : float
        The width in meters of the casing of the symbol.
    
    Returns
    -------
    Vector2D or None
        The translation vector, or None if it cannot be calculated.
    """
    # Calculate the angle perpendicular from the orientation of the bend
    o = bend.get_orientation()
    angle = o + np.pi / 2 if o < np.pi else o - np.pi / 2
    
    # Create the vector from the angle and the norm
    w = bend.get_width()
    v = Vector2D.from_angle(angle, w * 2)
    
    # Try first the start point, then the end point
    coords = list(bend.bend.coords)
    for base_coord in [coords[0], coords[-1]]:
        base = Point(base_coord)
        
        # Translate the base point using the vector
        tpoint = v.translate(base)
        # Calculate the line between the base point and the translated point
        tline = LineString([base, tpoint])
        # Extend the line by twice the width of the bend on both sides
        eline = extend_line_by_length(tline, w * 2)

        # Create the bend without the base vertex to avoid intersection there
        crossline_coords = [c for c in coords if c != base_coord]
        
        # Need at least 2 points for a LineString
        if len(crossline_coords) < 2:
            continue
            
        crossline = LineString(crossline_coords)
        # Calculate the intersection
        inter = intersection(eline, crossline)

        # Check that intersection is a point
        if inter.geom_type == 'Point':
            return v.change_norm(width)
    
    return None

def schematization(line, sigma=30, sample=None):
    """
    Remove bends from series of bends to simplify it.

    This algorithm proposed by Lecordix *et al.* :footcite:p:`lecordix:1997`
    is a caricature algorithm that removes one (or more) bend of a
    series while preserving the global shape of the bend series.

    Bends in the series are identified using inflexion points, and when the middle
    bends are removed, the inflexion points are displaced along the axis of the series,
    the road is distorted to cushion the displacement of the inflexion points.

    Parameters
    ----------
    line : LineString
        The line to apply the accordion algorithm.
    sigma : float, optional
        Gaussian smoothing strength.
    sample : int, optional
        Gaussian smoothing sample size.

    Returns
    -------
    LineString

    See Also
    --------
    accordion :
        Stretch a series of bends to enlarge each bend.
    smooth_gaussian :
        Smooth a line and attenuate its inflexion points.
    inflexion_points :
        Extract inflexion points from a sinuous line.

    References
    ----------
    .. footbibliography::
    """
    
    def schematize_part(coords, wpoint, summits, length):
        """
        Apply schematization to a part of the line.
        
        OPTIMIZATIONS:
        - Pre-compute segment lengths
        - Use set for O(1) summit lookup
        - Reduce redundant Point() creations
        - Cache vector calculations
        """
        # Get the start vertex and create the translation vector
        start = Point(coords[0])
        vector = Vector2D.from_points(start, wpoint)

        # Get the orientation of the first segment of the part
        main = Segment(coords[0], coords[1]).orientation()
        # Get the perpendicular orientation
        perp = main + np.pi / 2
        # Normalize to [0, 2π)
        if perp > (2 * np.pi):
            perp -= 2 * np.pi

        # Project the vector on the main and perpendicular orientation
        vmain = vector.project(main)
        vperp = vector.project(perp)

        # Convert summits to set for O(1) lookup instead of O(n)
        summits_set = set(summits)

        # Pre-compute segment lengths to avoid redundant calculations
        n_coords = len(coords)
        segment_lengths = np.zeros(n_coords - 1)
        for i in range(1, n_coords):
            segment_lengths[i - 1] = LineString([coords[i], coords[i - 1]]).length

        # Isolate the first half of the first bend of part 1
        lasthalf_set = set()
        halflength = 0
        # Loop through points starting at the penultimate vertex towards the start
        for i in range(1, n_coords):
            point = coords[i]
            halflength += segment_lengths[i - 1]
            # Break when reaching the 1st summit
            if point in summits_set:
                break
            lasthalf_set.add(point)

        schematized = []
        total = 0
        
        # Pre-compute constants
        inv_halflength = 1.0 / halflength if halflength > 0 else 0
        inv_length = 1.0 / length if length > 0 else 0
        vmain_norm = vmain.get_norm()
        vperp_norm = vperp.get_norm()
        
        # Loop again starting at the penultimate vertex
        for i in range(1, n_coords):
            point = coords[i]
            # Increment the total length with the current segment length
            total += segment_lengths[i - 1]
            
            # Get the current point
            projected = Point(point)
            
            # If the point is in the first half of the first bend
            if point in lasthalf_set:
                # Cushion the point using the main vector
                factor_main = vmain_norm * (1 - total * inv_halflength)
                vect = vmain.change_norm(factor_main)
                projected = vect.translate(projected)
            
            # Cushion the translation using the perpendicular vector
            factor_perp = vperp_norm * (1 - total * inv_length)
            vect = vperp.change_norm(factor_perp)
            projected = vect.translate(projected)
            
            # Add the point to the list
            schematized.append(projected)
        
        return schematized

    # Detect individual bends inside the smoothed line
    bs = BendSeries(line, sigma, sample)
    bends = bs.bends
    nb = len(bends)

    # Early return if there are less than 3 bends
    if nb < 3:
        return line

    # Pre-compute all bend metrics in a single pass
    lengths = np.zeros(nb)
    shapes = np.zeros(nb)
    sizes = np.zeros(nb)
    
    for i, b in enumerate(bends):
        lengths[i] = b.bend.length
        height, width = b.get_height(), b.get_width()
        shapes[i] = (4 * height + width) / 5
        sizes[i] = (4 * height * width + height * width) / 5

    # Compute totals using numpy for efficiency
    length_total = np.sum(lengths)
    shape_total = np.sum(shapes)
    size_total = np.sum(sizes)
    
    # Avoid division by zero
    inv_length = 1.0 / length_total if length_total > 0 else 0
    inv_shape = 1.0 / shape_total if shape_total > 0 else 0
    inv_size = 1.0 / size_total if size_total > 0 else 0

    # Calculate indices for all candidate bends
    indiceremove = []
    sizeremove = []

    for i, b in enumerate(bends):
        # Skip first and last bend unless there are exactly 5 bends
        if i in [0, nb - 1] and nb != 5:
            continue

        # Calculate the indice depending on the number of bends
        if nb < 6:
            indice = (2 * shapes[i] * inv_shape + sizes[i] * inv_size + lengths[i] * inv_length) / 4
        else:
            symetry = b.get_symmetry()
            if nb < 11:
                indice = (shapes[i] * inv_shape + sizes[i] * inv_size + lengths[i] * inv_length + symetry) / 4
            else:
                indice = (shapes[i] * inv_shape + 2 * sizes[i] * inv_size + 0.5 * lengths[i] * inv_length + 0.5 * symetry) / 4
        
        # Store bend information
        indiceremove.append({
            'bend': b,
            'index': i,
            'indice': indice
        })
        sizeremove.append({
            'bend': b,
            'index': i,
            'size': sizes[i] * inv_size
        })
    
    # Early return if not enough bends to remove
    if len(indiceremove) < 2:
        return line
        
    # Sort both lists by indice and size
    indiceremove.sort(key=lambda d: d['indice'])
    sizeremove.sort(key=lambda d: d['size'])

    # Handle particular case: remove largest bend from indiceremove
    if nb > 4:
        largest_index = sizeremove[-1]['index']
        # Use list comprehension instead of loop + pop
        indiceremove = [idx for idx in indiceremove if idx['index'] != largest_index]

    # Early return if not enough bends remain
    if len(indiceremove) < 2:
        return line

    # Find two consecutive bends to remove
    bend1, bend2 = None, None

    if len(indiceremove) == 2:
        bend1, bend2 = indiceremove[0], indiceremove[1]
    else:
        b1 = indiceremove[0]
        b1_idx = b1['index']
        
        # Use a more efficient search for consecutive bends
        for candidate in indiceremove[1:]:
            cid = candidate['index']
            if cid == b1_idx + 1:
                bend1, bend2 = b1, candidate
                break
            elif cid == b1_idx - 1:
                bend1, bend2 = candidate, b1
                break
    
    # If no consecutive bends found, return original line
    if bend1 is None or bend2 is None:
        return line

    # Build coordinate lists for both parts
    coords1, coords2 = [], []
    summit1, summit2 = [], []

    bend1_idx = bend1['index']
    bend2_idx = bend2['index']

    for i, b in enumerate(bends):
        if i < bend1_idx:
            summit1.append(b.get_bend_summit())
            coords1.extend(b.bend.coords)
        elif i > bend2_idx:
            summit2.append(b.get_bend_summit())
            coords2.extend(b.bend.coords)
    
    # Calculate length of part 1 and 2
    length1 = LineString(coords1).length
    length2 = LineString(coords2).length
    
    # Avoid division by zero
    if length2 == 0:
        return line
        
    ratio = length1 / length2

    # Get start and end point of the joining line
    start = Point(coords1[-1])
    end = Point(coords2[0])
    linejoin = LineString([start, end])

    # Create the vector from those two points
    v = Vector2D.from_points(start, end)

    # Weight the vector according to the ratio
    if ratio == 1:
        wpoint = linejoin.interpolate(linejoin.length / 2)
    elif ratio < 1:
        v.scalar_multiplication(ratio / 2)
        wpoint = v.translate(start)
    else:
        v.scalar_multiplication(3 / (2 * ratio))
        wpoint = v.translate(start)

    # Reverse the first part to match the order
    coords1.reverse()

    # Schematize both parts
    part1 = schematize_part(coords1, wpoint, summit1, length1)
    part2 = schematize_part(coords2, wpoint, summit2, length2)
    part1.reverse()

    coords = list(line.coords)
    
    # Return the constructed linestring, keeping the original start and end point
    return LineString([coords[0]] + part1[1:] + [wpoint] + part2[:-1] + [coords[-1]])