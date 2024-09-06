import numpy as np
import shapely
from shapely import LineString, Point

from cartagen.utils.math.vector import Vector2D
from cartagen.utils.geometry.bends import BendSerie
from cartagen.utils.geometry.line import extend_line_by_length, get_line_middle_point, merge_linestrings
from cartagen.utils.geometry.segment import Segment

def accordion(line, width, sigma=30, sample=None):
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
    gaussian_smoothing :
        Smooth a line and attenuate its inflexion points.

    References
    ----------
    .. footbibliography::
    """
    # Detect individual bends inside the smoothed line
    bs = BendSerie(line, sigma, sample)

    # Storage for the future distorded individual lines
    distorded = []
    
    for b in bs.bends:
        # Get the translation vector
        v = __get_vector(b, width)

        # If the vector could not be found, keep the original bend
        if v is None:
            distorded.append(b.bend)
            continue

        # Retrieve start and end of bend coordinates
        coords = list(b.bend.coords)
        start, end = Point(coords[0]), Point(coords[-1])

        # Calculate the vector formed by the start and end point of the bend
        vab = Vector2D.from_points(start, end)

        # Calculate the scalar product between the translation vector
        # and the vector formed by the start and end point of the bend
        scalar = v.scalar_product(vab)

        # If the scalar product is negative, inverse the translation vector
        if scalar < 0:
            v = v.opposite()

        # Calculate the summit
        summit = b.get_summit()

        # Calculate the vector formed by the start and the summit of the bend
        vas = Vector2D.from_points(start, summit)

        # Calculate the vector product of vas and vab
        product = vas.product(vab)

        # Check that the closed bend is clockwise
        clockwise = not shapely.is_ccw(LineString(coords + [coords[0]]))

        # Inverse the vector by comparating product and clockwise flag
        if clockwise and product > 0:
            v = v.opposite()
        elif clockwise is False and product < 0:
            v = v.opposite()

        # Storage for the translated points
        translated = []

        total = 0
        length = b.bend.length
        # Loop through coordinates of the bend vertices starting from the end
        for i in range(len(coords) - 1, -1, -1):
            # Get the coordinates of the current point
            p = Point(coords[i])
            # If it's the end point of the bend
            if i == len(coords) - 1:
                # Displace the point by the full vector
                translated.append(v.translate(p))
            else:
                # Here, it's not the end point
                # Increment the current processed length of the line
                total += Point(p).distance(Point(coords[i + 1]))
                # Copy the vector
                vdist = v.copy()
                # Calculate the portion of the processed length
                factor = (length - total) / length
                # Reduce the vector
                vdist.scalar_multiplication(factor)
                # Displace the point by the reduced vector and add it at the start of the list
                translated.insert(0, vdist.translate(p))
        
        # Add the list of translated points to the full list of the whole bend serie
        distorded.append(LineString(translated))

    # Set the displaced line to be the first distorded line
    displaced = distorded[0]
    # Set the next start to be the end of the first line
    lstart = displaced.coords[-1]

    # Loop through remaining distorded lines
    for i in range(1, len(distorded)):
        # Get the line
        dline = distorded[i]
        # Get the start of the line
        dstart = dline.coords[0]
        # Translate the line with dx and dy which are the difference between line start and end of previous line
        tline = shapely.affinity.translate(dline, lstart[0] - dstart[0], lstart[1] - dstart[1])
        # Merge the displaced lines with the current displaced line
        displaced = merge_linestrings(displaced, tline)
        # Set the last start as the current displaced line end point 
        lstart = displaced.coords[-1]

    # Get the middle of the provided bend serie
    bm = get_line_middle_point(line)
    # Get the middle of the distorded line
    dm = get_line_middle_point(displaced)

    # Translate the distorded line by the difference between both middle points
    return shapely.affinity.translate(displaced, bm.x - dm.x, bm.y - dm.y)

def __get_vector(bend, width):
    # Check on both sides if the translation crosses the bend
    def check_translation(base, bend, width, angle, v):
        # Translate the base point using the vector
        tpoint = v.translate(base)
        # Calculate the line between the base point and the translated point
        tline = LineString([Point(base), tpoint])
        # Extend the line by twice the width of the bend on both sides
        eline = extend_line_by_length(tline, bend.width * 2)

        # Create the bend without the base vertex to avoid intersection there
        crossline = LineString([ x for x in list(bend.bend.coords) if x != base ])
        # Calculate the intersection
        intersection = shapely.intersection(eline, crossline)

        # Check that intersection is a point.
        if intersection.geom_type == 'Point':
            # Get the distance between the base point and the intersection
            dist = Point(base).distance(intersection)
            v2 = v.change_norm(width - dist)
            return v2
        else:
            # If the interseection is not a point
            # Check that the end point has not been tried already
            if base == list(bend.bend.coords)[-1]:
                # Here, no intersection where found
                return None
            else:
                # Restart the function using the end point as the base
                return check_translation(Point(bend.bend.coords[-1]), bend, width, angle, v)

    # Calculate the angle perpendicular from the orientation of the bend
    o = bend.get_orientation()
    if o < np.pi:
        angle = o + np.pi / 2
    else:
        angle = o - np.pi / 2
   
    # Create the vector from the angle and the norm
    w = bend.get_width()
    v = Vector2D.from_angle(angle, w*2)

    return check_translation(Point(bend.bend.coords[0]), bend, width, angle, v)

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
    gaussian_smoothing :
        Smooth a line and attenuate its inflexion points.
    inflexion_points :
        Extract inflexion points from a sinuous line.

    References
    ----------
    .. footbibliography::
    """
    # Treat a part of the schematization
    def schematize_part(coords, wpoint, summits, length):
        # Get the start vertex and create the translation vector
        start = Point(coords[0])
        vector = Vector2D.from_points(start, wpoint)

        # Get the orientation of the first segment of the part
        main = Segment(coords[0], coords[1]).orientation()
        # Get the perpendicular orientation
        perp = main + np.pi / 2
        # Get it between 0 & 2pi
        if perp > (2 * np.pi):
            perp -= 2 * np.pi

        # Project the vector on the main and perpendicular orientation
        vmain = vector.project(main)
        vperp = vector.project(perp)

        # Isolate the first half of the first bend of part 1
        lasthalf = []
        halflength = 0
        # Loop through points starting at the penultimate vertex towards the start
        for i in range(1, len(coords)):
            # Get the point coordinates
            point = coords[i]
            # Calculate the length of the segment using the previous point
            halflength += LineString([point, coords[i - 1]]).length
            # Break when reaching the 1st summit
            if point in summits:
                break
            # Append the point to the list
            lasthalf.append(point)

        schematized = []
        total = 0
        # Loop again starting at the penultimate vertex
        for i in range(1, len(coords)):
            point = coords[i]
            # Increment the total length with the current segment length
            total += LineString([point, coords[i - 1]]).length
            # Get the current point
            projected = Point(point)
            # If the point is in the first half of the first bend
            if point in lasthalf:
                # Cushion the point using the main vector, i.e. the orientation of the last semgent
                vect = vmain.change_norm(vmain.get_norm() * (1 - total / halflength))
                projected = vect.translate(projected)
            # Then, for every point, cushion the translation using also the perpendicular vector
            vect = vperp.change_norm(vperp.get_norm() * (1 - total / length))
            projected = vect.translate(projected)
            # Add the point to the start of the list
            schematized.append(projected)
        
        return schematized

    schematized = []

    # Detect individual bends inside the smoothed line
    bs = BendSerie(line, sigma, sample)
    bends = bs.bends

    # If there are less than 3 bends, do not schematize
    if len(bends) < 3:
        return line
    
    # If there are more than 3 bends, keep the first and the last
    if len(bends) > 3:
        schematized = [ bends[0].bend, bends[-1].bend ]

    length = shape = size = 0
    # Loop through bends
    lengths, shapes, sizes = [], [], []
    for b in bends:
        # Increment by the length of the bend
        le = b.bend.length
        length += le
        lengths.append(le)
        # Increment the shape and size measures
        height, width = b.get_height(), b.get_width()
        sh = (4 * height + width) / 5
        shape += sh
        shapes.append(sh)
        si = (4 * height * width + height * width) / 5
        size += si
        sizes.append(si)

    # List for index to remove by indices and size
    indiceremove, sizeremove = [], []

    # Loop through bends
    nb = len(bends)
    for i, b in enumerate(bends):
        # Make sure the bend is to be schematized (not the first, not the last)
        if i in [0, nb - 1] and len(bends) != 5:
            continue

        indice = 0
        # Calculate the indice depending on the number of bends inside the bend serie
        if nb < 6:
            indice = (2 * shapes[i] / shape + sizes[i] / size + lengths[i] / length) / 4
        else:
            symetry = b.get_symetry()
            if nb < 11:
                indice = (shapes[i] / shape + sizes[i] / size + lengths[i] / length + symetry) / 4
            else:
                indice = (shapes[i] / shape + 2 * sizes[i] / size + 0.5 * lengths[i] / length + 0.5 * symetry) / 4
        
        # Add the bend with index, indice and size
        indiceremove.append({
            'bend': b,
            'index': i,
            'indice': indice
        })
        sizeremove.append({
            'bend': b,
            'index': i,
            'size': sizes[i] / size
        })
        
    # Order both list by indice and size
    indiceremove = sorted(indiceremove, key=lambda d: d['indice'])
    sizeremove = sorted(sizeremove, key=lambda d: d['size'])

    # Handle particular case
    if len(bends) > 4:
        i = sizeremove[-1]['index']
        topop = None
        for j, idx in enumerate(indiceremove):
            if idx['index'] == i:
                topop = j
        
        if topop is not None:
            indiceremove.pop(topop)

    bend1, bend2 = None, None

    # If the number of bends to remove is only one or zero,
    # schematization should not be used, returning the original line.
    if len(indiceremove) < 2:
        return line

    # If only 2 bends are to be remove, set them to be the removed ones
    if len(indiceremove) == 2:
        bend1 = indiceremove[0]
        bend2 = indiceremove[1]
    else:
        b1 = indiceremove[0]
        # Loop through bends to remove avoiding the first one
        for i in range(1, len(indiceremove)):
            # Get the original bend id
            cid = indiceremove[i]['index']
            # If this bend id is directly after the first one
            if cid == b1['index'] + 1:
                # Set both of them as bends to remove
                bend1 = b1
                bend2 = indiceremove[i]
                # Two bends to remove were found, breaking the loop
                break
            # If this bend id is directly before the first one, invert the two
            elif cid == b1['index'] - 1:
                bend1 = indiceremove[i]
                bend2 = b1
                break

    coords1, coords2 = [], []
    summit1, summit2 = [], []

    # Populate part 1 and part 2 of the discontinuous linestring
    for i, b in enumerate(bends):
        if i < bend1['index']:
            summit1.append(b.get_summit())
            coords1.extend(b.bend.coords)
            continue
        if i > bend2['index']:
            summit2.append(b.get_summit())
            coords2.extend(b.bend.coords)
    
    # Calculate length of part 1 and 2
    length1, length2 = LineString(coords1).length, LineString(coords2).length
    # Calculate the ratio of those two lengths
    ratio = length1 / length2

    # Get start and end point of the line formed by the last vertex 
    # of part 1 and the first vertex of part 2
    start, end = Point(coords1[-1]), Point(coords2[0])
    # Create the line
    linejoin = LineString([start, end])

    # Create the vector from those two points
    v = Vector2D.from_points(start, end)

    # Depending on the ratio of both lengths, weight the vector accordingly
    if ratio == 1:
        wpoint = linejoin.interpolate(linejoin.length / 2)
    elif ratio < 1:
        v.scalar_multiplication(ratio/2)
    else:
        v.scalar_multiplication(3/(2*ratio))

    # Create the anchor point
    wpoint = v.translate(start)

    # Reverse the first part to match the order
    coords1.reverse()

    # Schematize the first part
    part1 = schematize_part(coords1, wpoint, summit1, length1)
    part2 = schematize_part(coords2, wpoint, summit2, length2)
    part1.reverse()

    coords = list(line.coords)
    
    # Return the constructed linestring, keeping the original start and end point
    return LineString([coords[0]] + part1[1:] + [wpoint] + part2[:-1] + [coords[-1]])