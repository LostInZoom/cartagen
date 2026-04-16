import numpy as np
from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

from cartagen.utils.geometry.line import resample_line

def smooth_gaussian(geometry, sigma=30, sample=None, densify=True):
    """
    Smooth a line or a polygon and attenuate its inflexion points.

    The gaussian smoothing has been studied by Babaud *et al.* :footcite:p:`babaud:1986`
    for image processing, and by Plazanet :footcite:p:`plazanet:1996`
    for the generalisation of cartographic features.

    Accept Multi geometries. If a polygon is provided,
    it also apply the smoothing to its holes using the same parameters.

    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to smooth.
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
    LineString, Polygon, MultiLineString, MultiPolygon

    See Also
    --------
    smooth_catmull_rom :
        Smooth a line or polygon and preserve vertexes.
    smooth_chaikin :
        Smooth a line or polygon by cutting corners.
    smooth_platre :
        Smooth a line and preserve the integrity of sharp turns.
    smooth_taubin :
        Smooth a line or polygon and prevent shrinkage.
    smooth_wma :
        Smooth a line or polygon using a low-pass filter.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.smooth_gaussian(line, 1)
    <LINESTRING (0 0, 1.666666666666667 0.6051115971014416, 3.333333333333334 1.6051115971014418, 5 3)>

    >>> polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    >>> c4.smooth_gaussian(polygon, 1)
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

    # --- 1. Recursive handling for Multi-geometries ---
    if geometry.geom_type == 'MultiLineString':
        geoms = [smooth_gaussian(geometry, sigma, sample, densify) for g in geometry.geoms]
        return MultiLineString(geoms)

    if geometry.geom_type == 'MultiPolygon':
        geoms = [smooth_gaussian(geometry, sigma, sample, densify) for g in geometry.geoms]
        return MultiPolygon(geoms)

    # --- 2. Handling Polygons ---
    if geometry.geom_type == 'Polygon':
        # Simplify the exterior ring
        simplified_exterior = smooth_gaussian(geometry.exterior, sigma, sample, densify)
        
        # VALIDITY CHECK: A LinearRing must have at least 4 coordinates
        if len(simplified_exterior.coords) < 4:
            # Option A: Return an empty geometry (it will be filtered out in the main loop)
            return Polygon() 
            # Option B: return geom (if you want to keep the original instead of deleting it)
        
        # Ensure it's a LinearRing (closed)
        exterior_ring = LinearRing(simplified_exterior.coords)

        # Simplify interior rings (holes)
        simplified_interiors = []
        for interior in geometry.interiors:
            s_int = smooth_gaussian(geometry, sigma, sample, densify)
            # Only keep holes that are still large enough to be rings
            if len(s_int.coords) >= 4:
                simplified_interiors.append(LinearRing(s_int.coords))
        
        poly = Polygon(exterior_ring, simplified_interiors)
        
        # Final topological repair
        if not poly.is_valid:
            poly = poly.buffer(0)
        return poly

    # --- 3. Core Linear Logic (for LineString, LinearRing, etc.) ---
    if geometry.geom_type not in ['LineString', 'LinearRing']:
        raise ValueError(f'{geometry.geom_type} geometry type cannot be simplified.')

    coords = list(geometry.coords)

    if sample is None:
        distances = []
        for i in range(0, len(coords) - 1):
            v1, v2 = coords[i], coords[i + 1]
            distances.append(Point(v1).distance(Point(v2)))
        avg = (sum(distances) / len(distances))
        sample = avg

    # First resample the line, making sure there is a maximum distance between two consecutive vertices
    resampled = resample_line(geometry, sample)

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

    if geometry.geom_type == 'LinearRing':
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
    if geometry.geom_type == 'LinearRing':
        final_coords.append(final_coords[0])
    else:
        # Replace first and last vertex by the line's original ones
        final_coords[0] = Point(coords[0])
        final_coords[-1] = Point(coords[-1])

    return LineString(final_coords)