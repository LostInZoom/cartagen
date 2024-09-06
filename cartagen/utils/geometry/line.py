import numpy as np
import shapely
import geopandas as gpd
from shapely.ops import split, nearest_points, snap, transform
from shapely.geometry import LineString, Point, Polygon, MultiPoint

from cartagen.utils.geometry.angle import angle_3_pts
from cartagen.utils.partitioning.tessellation import HexagonalTessellation
from cartagen.utils.partitioning import partition_grid

def douglas_peucker(line, threshold, preserve_topology=True):
    """
    Distance-based line simplification.

    This algorithm was proposed by Ramer :footcite:p:`ramer:1972` and by Douglas and Peucker.
    :footcite:p:`douglas:1973` It is a line filtering algorithm, which means that it
    filters the vertices of the line (or polygon) to only retain the most important ones
    to preserve the shape of the line. The algorithm iteratively searches the most
    characteristics vertices of portions of the line and decides to retain
    or remove them given a distance threshold.

    The algorithm tends to unsmooth geographic lines, and is rarely used to simplify geographic features.
    But it can be very useful to quickly filter the vertices of a line inside another algorithm.

    This is a simple wrapper around :func:`shapely.simplify() <shapely.simplify()>`.

    Parameters
    ----------
    line : LineString
        The line to simplify.
    threshold : float
        The distance thresholdto remove the vertex from the line.
    preserve_topology : bool, optional
        If set to True, the algorithm will prevent invalid geometries
        from being created (checking for collapses, ring-intersections, etc).
        The trade-off is computational expensivity.

    Returns
    -------
    LineString

    See Also
    --------
    visvalingam_whyatt :
        Area-based line simplification.
    raposo :
        Hexagon-based line simplification.
    li_openshaw :
        Square grid-based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> douglas_peucker(line, 1.0)
    <LINESTRING (0 0, 2 0, 5 3)>
    """
    return shapely.simplify(line, threshold, preserve_topology=preserve_topology)

def visvalingam_whyatt(line, area_tolerance):
    """
    Area-based line simplification.

    This algorithm proposed by Visvalingam and Whyatt :footcite:p:`visvalingam:1993` performs a
    line simplification that produces less angular results than the filtering algorithm of Ramer-Douglas-Peucker.
    The principle of the algorithm is to select the vertices to delete (the less characteristic ones)
    rather than choosing the vertices to keep (in the Douglas and Peucker algorithm).
    To select the vertices to delete, there is an iterative process,
    and at each iteration, the triangles formed by three consecutive vertices are computed. If the area of the smallest
    triangle is smaller than a threshold (“area_tolerance” parameter), the middle vertex is deleted, and another iteration starts.

    The algorithm is relevant for the simplification of natural line or polygon features such as rivers, forests, or coastlines.

    Parameters
    ----------
    line : LineString
        The line to simplify.
    area_tolerance : float
        The minimum triangle area to keep a vertex in the line.

    Returns
    -------
    LineString

    See Also
    --------
    douglas_peucker :
        Distance-based line simplification.
    raposo :
        Hexagon-based line simplification.
    li_openshaw :
        Square grid-based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> visvalingam_whyatt(line, 5.0)
    <LINESTRING (0 0, 2 0, 5 3)>
    """
    def contains_another_point(line, pt, index):
        first = line.coords[index-1]
        last = line.coords[index+1]
        triangle_coords = [first, pt, last, first]
        triangle = Polygon(triangle_coords)
        for vertex in line.coords:
            if (vertex == first):
                continue
            if (vertex == pt):
                continue
            if (vertex == last):
                continue
            if (triangle.contains(Point(vertex))):
                return True
        return False
    def compute_area_point(line, pt, index):
        first = line.coords[index-1]
        last = line.coords[index+1]
        triangle = [first, pt, last, first]
        polygon = Polygon(triangle)
        return polygon.area

    final_coords = list(line.coords)
    while (True):
        pt_area_min = None
        area_min = line.envelope.area
        current_line = LineString(tuple(final_coords))
        for i in range(1,len(final_coords)-2):
            pt = final_coords[i]
            if(contains_another_point(current_line, pt, i)):
                continue
            area = compute_area_point(current_line, pt, i)
            if (area < area_min):
                area_min = area
                if(area < area_tolerance):
                    pt_area_min = pt
        if(pt_area_min is None):
            break
        # remove pt from the line
        final_coords.remove(pt_area_min)

    return LineString(tuple(final_coords))

def raposo(line, initial_scale, final_scale, centroid=True, tobler=False):
    """
    Hexagon-based line simplification.
    
    This algorithm proposed by Raposo :footcite:p:`raposo:2013` simplifies lines based on a
    hexagonal tessellation. The algorithm also works for the simplification of the border of a polygon object.
    The idea of the algorithm is to put a hexagonal tessallation on top of the line to simplify,
    the size of the cells depending on the targeted granularity of the line.
    Similarly to the Li-Openshaw algorithm, only one vertex is kept inside each cell.
    This point can be the centroid of the removed vertices, or a projection on the initial line of this centroid.
    The shapes obtained with this algorithm are less sharp than the ones obtained with other algorithms such as Douglas-Peucker.

    The algorithm is dedicated to the smooth simplification of natural features such as rivers, forests, coastlines, lakes.

    Parameters
    ----------
    line : LineString
        The line to simplify.
    initial_scale : float
        Initial scale of the provided line (25000.0 for 1:25000 scale).
    final_scale : float
        Final scale of the simplified line.
    centroid : bool, optional
        If true, uses the center of the hexagonal cells as the new vertex,
        if false, the center is projected on the nearest point in the initial line.
    tobler : bool, optional
        If True, compute cell resolution based on Tobler’s formula, else uses Raposo's formula

    Returns
    -------
    LineString   

    See Also
    --------
    douglas_peucker :
        Distance-based line simplification.
    visvalingam_whyatt :
        Area-based line simplification.
    li_openshaw :
        Square grid-based line simplification.

    Notes
    -----
    The Tobler based formula to compute hexagonal cell size is :math:`cellsize=5·l·s`
    where :math:`l` is the width of the line in the map in meters
    (*e.g.* 0.0005 for 0.5 mm), and :math:`s` is the target scale denominator.

    Raposo’s formula to compute hexagonal cell size is :math:`cellsize={l/n}·{t/d}`
    where :math:`l` is the length of the line, :math:`n` the number of vertices of the line,
    :math:`t` the denominator of the target scale, and :math:`d` the denominator of the initial scale.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.raposo(line, 5000, 10000)
    <LINESTRING (0 0, 1 0.3333333333333333, 5 3)>
    """
    width = 0
    if tobler:
        width = final_scale * 5 / 4000 
    else:
        firstFactor = line.length / (len(line.coords)-1)
        secondFactor = final_scale / initial_scale
        width = firstFactor * secondFactor

    # compute hexagon tessellation
    tessellation = HexagonalTessellation(line.envelope, width)

    current_index = 0
    final_coords = []
    # append the first point of the line, but without the z coordinate
    twod_point = transform(lambda *args: args[:2], Point(line.coords[0]))
    final_coords.append(twod_point.coords[0])
    previous_cell = None
    while (current_index < len(line.coords)-1):
        current_cell = None
        # now loop on the vertices from current index
        # builds a point cloud as a multi-point geometry with all line vertices
        # contained in the current cell.
        point_cloud = []
        for i in range(current_index,len(line.coords)):
            # get the cells containing the point
            point = line.coords[i]
            containing_cells = tessellation.get_containing_cells(point)
            if(current_cell is None):
                if(previous_cell in containing_cells):
                    containing_cells.remove(previous_cell)
                current_cell = containing_cells[0]
                point_cloud.append(point)
                continue

            if (current_cell in containing_cells):
                point_cloud.append(point)
                current_index = i+1
            else:
                current_index = i
                previous_cell = current_cell
                break
        multipoint = MultiPoint(point_cloud)
        if (centroid):
            # replace the points by the centroid of the vertices in the cell
            final_coords.append(multipoint.centroid.coords[0])
        else:
            # find the nearest vertex to the centroid
            nearest = nearest_points(multipoint,multipoint.centroid)
            final_coords.append(nearest[0].coords[0])
        if (current_index == len(line.coords) - 1):
            twod_final_point = transform(lambda *args: args[:2], Point(line.coords[current_index]))
            final_coords.append(twod_final_point.coords[0])
    # add the final point if it is not in the line
    twod_final_point = transform(lambda *args: args[:2], Point(line.coords[len(line.coords) - 1]))
    if(twod_final_point.coords[0] not in final_coords):
        final_coords.append(twod_final_point.coords[0])
    
    # print(final_coords)
    return LineString(final_coords)

def li_openshaw(line, cell_size):
    """
    Regular grid-based line simplification.

    This algorithm proposed by Li & Openshaw :footcite:p:`li:1993` simplifies lines based on a
    regular square grid. It first divide the line vertexes into groups partionned by a regular
    grid, then each group of vertexes is replaced by their centroid.

    Parameters
    ----------
    line : LineString
        The line to simplify.
    cell_size : float
        The size of the regular grid used to divide the line.

    Returns
    -------
    LineString   

    See Also
    --------
    douglas_peucker :
        Distance-based line simplification.
    visvalingam_whyatt :
        Area-based line simplification.
    raposo :
        Hexagon-based line simplification.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.li_openshaw(line, 1)
    <LINESTRING (0 0, 0.5 0.5, 2 0, 5 3)>
    """
    vertexes = [ Point(x) for x in list(line.coords) ]

    gdf = gpd.GeoDataFrame(geometry=vertexes)
    groups, squares = partition_grid(gdf, cell_size)

    simplified = []
    previous = None
    for i in range(0, len(vertexes)):
        index = None
        for gi, j in enumerate(groups):
            if i in j:
                index = gi
        
        if index is not None and index != previous:
            group = groups[index]
            centroid = shapely.centroid(MultiPoint([ vertexes[v] for v in group ]))
            simplified.append(centroid)
            previous = index

    if vertexes[0] != simplified[0]:
        simplified.insert(0, vertexes[0])
    if vertexes[-1] != simplified[-1]:
        simplified.append(vertexes[-1])

    return LineString(simplified)

def gaussian_smoothing(geometry, sigma=30, sample=None, densify=True):
    """
    Smooth a line or a polygon and attenuate its inflexion points.

    The gaussian smoothing has been studied by Babaud *et al.* :footcite:p:`babaud:1986`
    for image processing, and by Plazanet :footcite:p:`plazanet:1996`
    for the generalisation of cartographic features.

    Parameters
    ----------
    geometry : LineString or Polygon
        The line or polygon to smooth.
        If a line is provided, the first and last vertexes are kept.
        If a polygon is provided, every vertex is smoothed.
    sigma : float, optional
        Gaussian filter strength. Default value to 30, which is a high value.
    sample : float, optional
        The length in meter between each nodes after resampling the geometry.
        If not provided, the sample is derived from the geometry and is the average distance between
        each consecutive vertex.
    densify : bool, optional
        Whether the resulting geometry should keep the new vertex density. Default to True.

    Returns
    -------
    LineString or Polygon

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.gaussian_smoothing(line, 1)
    <LINESTRING (0 0, 1.666666666666667 0.6051115971014416, 3.333333333333334 1.6051115971014418, 5 3)>

    >>> polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    >>> c4.gaussian_smoothing(polygon, 1)
    <POLYGON ((0.1168459780814714 0.3005282653219513, ... 0.1168459780814714 0.3005282653219513))>
    """
    # Extend the given set of points at its first and last points of k points using central inversion.
    def extend(line, interval):
        # Compute the central inversion of a position. origin is the center of symmetry and p is the point to inverse.
        def central_inversion(origin, p):
            x = 2 * origin[0] - p[0]
            y = 2 * origin[1] - p[1]
            return (x,y)
        
        # Get the coordinates of the vertices
        coords = list(line.coords)
        # Get the first and last vertex
        first, last = coords[0], coords[-1]

        # Get the index of the penultimate vertex
        # -2 is to avoid taking the last vertex
        pen = len(coords) - 2

        # Set the start of the line as the central inversion of n first vertices (n = interval)
        result = [central_inversion(first, coords[i]) for i in range(interval, 0, -1)]

        # Add the full line as the middle part of the line
        result.extend(coords)

        # Add the end of the line as the central inversion of n last vertices (n = interval)
        result.extend([central_inversion(last, coords[i]) for i in range(pen, pen - interval, -1)])

        return LineString(result)

    polygon = None
    ring = None
    geomtype = geometry.geom_type
    if geomtype == 'LineString':
        polygon = False
        ring = geometry
    elif geomtype == 'Polygon':
        polygon = True
        ring = geometry.exterior
    else:
        raise Exception("{0} geometry cannot be smoothed.".format(geomtype))

    coords = list(ring.coords)

    if sample is None:
        distances = []
        for i in range(0, len(coords) - 1):
            v1, v2 = coords[i], coords[i + 1]
            distances.append(shapely.Point(v1).distance(shapely.Point(v2)))
        avg = (sum(distances) / len(distances))
        sample = avg

    # First resample the line, making sure there is a maximum distance between two consecutive vertices
    resampled = resample_line(ring, sample)

    # Calculate the interval (number of vertex to take into consideration when smoothing)
    interval = round(4 * sigma / sample)
    # If the interval is longer than the input line, we change the interval and recalculate the sigma
    if interval >= len(resampled.coords):
        interval = len(resampled.coords) - 1
        sigma = interval * sample / 4
    
    # Compute gaussian coefficients
    c2 = -1.0 / (2.0 * sigma * sigma)
    c1 = 1.0 / (sigma * np.sqrt(2.0 * np.pi))

    # Compute the gaussian weights and their sum
    weights = []
    total = 0
    for k in range (0, interval + 1):
        weight = c1 * np.exp(c2 * k * k)
        weights.append(weight)
        total += weight
        if k > 0:
            total += weight
    
    rline = list(resampled.coords)
    length = len(rline)

    if polygon:
        extended = LineString(rline[-interval:] + rline + rline[:interval])
    else:
        # Extend the line at its first and last points with central inversion
        extended = extend(resampled, interval)
   
    smoothed_coords = []
    for i in range(0, length):
        x, y = 0, 0
        for k in range(-interval , interval + 1):
            p1 = extended.coords[i - k + interval]
            x += weights[abs(k)] * p1[0] / total
            y += weights[abs(k)] * p1[1] / total
        smoothed_coords.append((x,y))

    if densify:
        final_coords = smoothed_coords
    else:
        # Only return the points matching the input points in the resulting filtered line
        final_coords = []
        # Stores for index of already treated vertices
        done = []
        # Loop through initial vertices
        for point in coords:
            # Set the distance to infinite
            distance = float("inf")
            nearest = None
            # Loop through smoothed coordinates
            for i in range(len(smoothed_coords)):
                # Check that the index has not been already added
                if i not in done:
                    # Calculate distance from the point
                    d = Point(smoothed_coords[i]).distance(Point(point))
                    if d < distance:
                        # Update distance and nearest index if below existing
                        distance, nearest = d, i

            # If a nearest point has been found, add it to the new line
            if nearest is not None:
                final_coords.append(smoothed_coords[nearest])
                # Add the index as treated already
                done.append(nearest)
            else:
                final_coords.append(point)
    
    result = None
    if polygon:
        final_coords.append(final_coords[0])
        result = Polygon(final_coords)
    else:
        # Replace first and last vertex by the line's original ones
        final_coords[0] = Point(coords[0])
        final_coords[-1] = Point(coords[-1])
        result = LineString(final_coords)
    
    return result

def get_bend_side(line):
    """
    Return the side of the interior of the bend.

    Parameters
    ----------
    line : LineString
        The line to get the bend side from.

    Returns
    -------
    side : str
        'right' or 'left'

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0)])
    >>> get_bend_side(line)
    right
    """
    # Get the list of nodes
    coords = list(line.coords)

    # The total angle of the bend
    total = 0

    # Get the start point of the bend
    start = shapely.Point(coords[0])

    # Loop through the nodes of the linestring, starting on the seconde node
    for i in range(1, len(coords) - 1):
        # Get the current and the next node coordinates
        c1 = coords[i]
        c2 = coords[i + 1]

        # Add to the total the angle between the starting point and those two nodes
        angle = angle_3_pts(start, shapely.Point(c2), shapely.Point(c1))

        # Set the angle to be between 0 and 2*pi
        if angle % (2 * np.pi) >= 0:
            angle = abs(angle % (2 * np.pi))
        else:
            angle = angle % (2 * np.pi) + 2 * np.pi

        # Set the angle to be between -pi and pi
        if angle > np.pi:
            angle = angle - 2 * np.pi

        total += angle

    # If the total is above 0, the bend is left sided, otherwise it right sided
    if total < 0:
        return 'right'
    else:
        return 'left'

def get_shortest_edge_length(geom):
    min_length = float('inf')
    segments = get_linestring_segments(geom)
    for segment in segments:
        print(segment)
        length = Point(segment[0]).distance(Point(segment[1]))
        if(length < min_length):
            min_length = length

    return min_length

def get_linestring_segments(line):
    segments = []
    prev_coord = line.coords[0]
    for coord in line.coords[1:]:
        segments.append((prev_coord,coord))
        prev_coord = coord

    return segments


# for a 3d sequence of coordinates, returns the 2d sequence. Useful to convert 3d points to 2d points
def to_2d(x, y, z):
    return tuple(filter(None, [x, y]))


# for an array of 3d sequences of coordinates, returns an array containing the 2d sequences. Useful to convert 3d geometries to 2d geometries
def array_to_2d(array):
    array_2d = []
    for tuple in array:
       array_2d.append(to_2d(tuple[0],tuple[1],tuple[2]))
    return array_2d

def polygons_3d_to_2d(polygons):
    polygons_2d = []
    for polygon in polygons:
       # outer ring
       coords_2d = array_to_2d(polygon.exterior.coords)
       # inner rings
       interiors = []
       for interior in polygon.interiors:
          interiors.append(array_to_2d(interior.coords))
       polygons_2d.append(Polygon(coords_2d,interiors))
    return polygons_2d

def resample_line(line, step, keep_vertices=False):
    """
    Densify a line by adding vertices.

    This function densifies a line by adding a vertex every 'step' along the line.
    It preserves the first and last vertex of the line.
    
    Parameters
    ----------
    line : LineString
        The line to densify.
    step : float
        The step (in meters) to resample the geometry.
    keep_vertices : bool, optional
        If set to true, original vertices of the line are kept.
        This is useful to keep the exact same geographical shape.

    Returns
    -------
    LineString

    Examples
    --------
    >>> line = LineString([(1, 1), (5, 1)])
    >>> resample_line(line, 1)
    <LINESTRING (1 1, 2 1, 3 1, 4 1, 5 1)>
    """

    coords = list(line.coords)
    original = []

    if keep_vertices:
        for i in range(1, len(coords) - 1):
            c = coords[i]
            d = shapely.line_locate_point(line, shapely.Point(c))
            original.append((c, d))
    
    # Get the length of the line
    length = line.length 
    
    # Storage for final vertices, starting with the start of the line
    xy = [(line.coords[0])]
    
    for distance in np.arange(step, int(length), step):
        if keep_vertices:
            remove = 0
            for i, o in enumerate(original):
                if o[1] < distance:
                    xy.append(o[0])
                    remove += 1
            for r in range(0, remove):
                original.pop(0)

        # Interpolate a point every step along the old line
        point = line.interpolate(distance)
        # Add the tuple of coordinates
        xy.append((point.x, point.y))
    
    # Add the last point of the line if it doesn't already exist
    if xy[-1] != line.coords[-1]:
        xy.append(line.coords[-1])
        
    # Here, we return a new line with densified points.
    return LineString(xy)

# returns the index of a vertex in a line. Returns -1 if the point given is not a vertex of the line
def get_index_of_vertex(line, vertex, tolerance = 0.01):
    for i in range(0, len(line.coords)):
        point = Point(line.coords[i])
        if point.equals_exact(vertex, tolerance):
            return i
    return -1

# returns the index of the nearest vertex in a line. Returns -1 if the point given is not a vertex of the line
def get_index_of_nearest_vertex(line, vertex):
    mindist = float("inf")
    nearest_index = -1
    for i in range(0, len(line.coords)):
        point = Point(line.coords[i])
        dist = point.distance(vertex)
        if dist < mindist:
            nearest_index = i
            mindist = dist
    return nearest_index

def geometry_flatten(geom):
  if hasattr(geom, 'geoms'):  # Multi<Type> / GeometryCollection
    for g in geom.geoms:
      yield from geometry_flatten(g)
  elif hasattr(geom, 'interiors'):  # Polygon
    yield geom.exterior
    yield from geom.interiors
  else:  # Point / LineString
    yield geom

def geometry_len(geom):
  return sum(len(g.coords) for g in geometry_flatten(geom))

# Gets the nearest point of the geometry to a point. If the point is a vertex, it is not chosen.
def get_nearest_vertex(geometry, point):
    min_dist = float('inf')
    nearest = None
   
    for vertex in geometry.coords:
        vertex_pt = Point(vertex)
        distance = vertex_pt.distance(point)
        if distance == 0:
            continue
        if distance < min_dist:
            min_dist = distance
            nearest = vertex_pt
    
    return nearest

def extend_line_with_point(line, point, position='start'):
    """
    Extend the line with a given point, depending on the position, adds it at the start or the end
    """
    new_line = []
    p = [point.x, point.y]

    if position == 'start':
        new_line.append(p)

    for n in line.coords:
        new_line.append(n)

    if position == 'end':
        new_line.append(p)

    return LineString(new_line)

def extend_line_by_length(line, length, position='both'):
    """
    Extend the line by a given length.
    This function extend the line by adding a vertice at the start and/or the end of the line and extend it by a given length
    following the direction of the first and/or last segment of the line.
    Position is by default set to 'both' which will add the length to both
    the start and the end of the line. Can be set to 'start' or 'end'.
    """

    def __calculate_coords(a, b, length):
        xa, ya, xb, yb = a[0], a[1], b[0], b[1]
        ab = np.sqrt((yb - ya)**2 + (xb -xa)**2)
        xc = xb - (((xb - xa) / ab) * (ab + length))
        yc = yb - (((yb - ya) / ab) * (ab + length))
        return (xc, yc)

    # Get the list of vertex
    vertex = list(line.coords)

    new_line = []

    if position == 'both' or position == 'start':
        new_line.append(__calculate_coords(vertex[0], vertex[1], length))

    new_line.extend(vertex)

    if position == 'both' or position == 'end':
        new_line.append(__calculate_coords(vertex[-1], vertex[-2], length))

    return LineString(new_line)

def get_line_middle_point(line):
    """
    Return the point located along the line at half its length.
    """
    return line.interpolate(line.length / 2)

        
def split_line_at_point(line, point):
    """
    Split a line at a given point along this line.
    Return the two new linestrings.
    Return None if the line and the point doesn't intersect.
    """
    projected = nearest_points(point, line)[1]
    splitted = split(snap(line, projected, 0.0001), projected)

    lines = []
    for s in splitted.geoms:
        lines.append(s)

    return lines[0], lines[1]

def get_segment_center(segment):
    """
    Return the center of the given shapely segment (a linestring with only two vertices) as a shapely Point.
    """
    c = segment.coords

    if len(c) > 2:
        raise Exception("The provided segment has more than two vertices.")

    x1, y1, x2, y2 = c[0][0], c[0][1], c[1][0], c[1][1]
    return shapely.Point([(x1 + x2) / 2, (y1 + y2) / 2])

def merge_linestrings(line1, line2):
    """
    Merge two linestring that are connected by their start or end points.
    Return the new linestring. The first line defines the direction of the merged linestring.
    """
    coords1, coords2 = line1.coords, line2.coords

    start1, end1 = coords1[0], coords1[-1]
    start2, end2 = coords2[0], coords2[-1]
    
    # Lines connected by...
    # ...start and start...
    if start1 == start2:
        # Reverse the second line, and concatenate the first (drop the start)
        return shapely.LineString(list(coords2[::-1]) + list(coords1[1:]))
    # ...start and end...
    elif start1 == end2:
        # Concatenate the second and the first line (drop the start)
        return shapely.LineString(list(coords2) + list(coords1[1:]))
    # ...end and start...
    elif end1 == start2:
        # Concatenate the first and the second line (drop the start)
        return shapely.LineString(list(coords1) + list(coords2[1:]))
    # ...end and end
    elif end1 == end2:
        # Reverse the second line (drop the start) and add it to the first line
        return shapely.LineString(list(coords1) + list(coords2[::-1][1:]))
    # Here, linstrings are not connected
    else:
        raise Exception("Provided LineStrings are not connected by their start or end vertex.")

def inflexion_points(line, min_dir=120.0):
    """
    Detect inflexion points inside a curved line.

    This algorithm extract inflexion points from a line using
    angles calculation. Micro inflexions are removed.

    Parameters
    ----------
    line : LineString
        The line to extract the inflexion points from.
    min_dir : float, optional
        The minimum direction change (in degrees) between two consecutive inflexion points.
        This parameter allows to remove micro inflexions from the results. If set to 0,
        every micro-inflexions will be kept.

    Returns
    -------
    list of int
        A list of index of the line vertices considered to be inflexion points
    
    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3), (8, 5)])
    >>> c4.inflexion_points(line)
    [2, 3]
    """

    coords = list(line.coords)

    # Storage for inflexion points
    inflexion = []

    # Storage for parts to check when micro inflexion points are detected
    part = []

    # Stores for the previous angle and the previous direction
    prevangle = None
    prevdir = None

    # Loop through vertices
    for i in range(1, len(coords) - 1):
        # Get previous, current and next point
        p1, p2, p3 = coords[i - 1], coords[i], coords[i + 1]

        # Calculate angles formed by p1 and p2, and p1 and p3
        a1 = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])
        a2 = np.arctan2(p3[1] - p1[1], p3[0] - p1[0])

        # Calculate the direction of p1 p2
        direction = np.degrees(a1) % 360
        # Calculate the angle formed by the three points
        alpha = a2 - a1

        # Store if angle is positive or negative
        angle = 1 if alpha > 0 else -1

        # Check that previous angle has been calculated
        if prevangle is not None:
            # Check if the angle is not the same as the previous -> inflexion point
            if angle != prevangle:
                # Check that the previous direction has been calculated
                if prevdir is not None:
                    # If the absolute difference between the previous direction and
                    # the current one is above the direction threshold
                    # It means that this is not a micro inflexion
                    if abs(prevdir - direction) > min_dir:
                        # Adding the middle of the part list, i.e. the middle of the micro inflexions
                        inflexion.append(part[len(part) // 2])
                        # Restart the part list with the current inflexion point
                        part = [i]
                    else:
                        # Add the current point to the part list
                        part.append(i)
                # If not, append the point to the part
                else:
                    part.append(i)

                # Set the previous angle as the direction of the current 
                prevdir = direction
        
        # Set previous angle as current
        prevangle = angle

    # Append the last part's middle point
    inflexion.append(part[len(part) // 2])

    return inflexion
