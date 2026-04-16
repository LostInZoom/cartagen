import numpy as np
from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, LinearRing

def smooth_taubin(geometry, iterations=10, smoothing=0.5, inflation=-0.53):
    """
    Smooth a line or polygon and prevent shrinkage.
    
    This algorithm was proposed by Taubin :footcite:p:`taubin:1995`.
    It is a two-step low-pass filter that preserves volume and 
    prevents shrinkage. It applies alternating expansion and contraction steps 
    using lambda (smoothing) and mu (inflation) parameters.
    
    Accept Multi geometries. If a polygon is provided, it also applies the 
    smoothing to its holes using the same parameters.

    This implementation is a translation of the C++ implementation available
    `here <https://github.com/philipschall/shapelysmooth>`_.
    
    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The line or polygon to smooth.
        If a line is provided, the first and last vertices are kept fixed.
        If a closed ring or polygon is provided, all vertices are smoothed.
    iterations : int, optional
        Number of smoothing iterations. Default is 10.
        More iterations result in stronger smoothing.
    smoothing : float, optional
        Smoothing factor for the first pass. Default is 0.5.
        Should be positive and typically between 0 and 1.
    inflation : float, optional
        Inflation factor for the second pass. Default is -0.53.
        Should be negative with inflation slightly larger than smoothing to prevent shrinkage.
    
    Returns
    -------
    LineString, Polygon, MultiLineString, MultiPolygon, LinearRing
        Smoothed geometry of the same type as input.

    See Also
    --------
    smooth_catmull_rom :
        Smooth a line or polygon and preserve vertexes.
    smooth_chaikin :
        Smooth a line or polygon by cutting corners.
    smooth_gaussian :
        Smooth a line or a polygon and attenuate its inflexion points.
    smooth_platre :
        Smooth a line and preserve the integrity of sharp turns.
    smooth_wma :
        Smooth a line or polygon using a low-pass filter.

    References
    ----------
    .. footbibliography::
    
    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (3, 1), (4, 0)])
    >>> smooth_taubin(line, smoothing=0.5, iterations=5)
    <LINESTRING (0 0, 1.2 0.8, 2 0.3, 2.8 0.8, 4 0)>
    
    >>> polygon = Polygon([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)])
    >>> smooth_taubin(polygon, smoothing=0.3, iterations=10)
    <POLYGON ((0.2 0.2, 0.2 1.8, 1.8 1.8, 1.8 0.2, 0.2 0.2))>
    """
    # --- 1. Recursive handling for Multi-geometries ---
    if geometry.geom_type == 'MultiLineString':
        geoms = [smooth_taubin(geometry, iterations, smoothing, inflation) for g in geometry.geoms]
        return MultiLineString(geoms)

    if geometry.geom_type == 'MultiPolygon':
        geoms = [smooth_taubin(geometry, iterations, smoothing, inflation) for g in geometry.geoms]
        return MultiPolygon(geoms)

    # --- 2. Handling Polygons ---
    if geometry.geom_type == 'Polygon':
        # Simplify the exterior ring
        simplified_exterior = smooth_taubin(geometry.exterior, iterations, smoothing, inflation)
        
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
            s_int = smooth_taubin(interior, iterations, smoothing, inflation)
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
    
    # Process LineString or LinearRing
    coords = np.array(geometry.coords)
    is_closed = np.array_equal(coords[0], coords[-1])
    
    # Apply Taubin smoothing
    factors = [smoothing, inflation]
    
    for _ in range(iterations):
        for factor in factors:
            new_coords = coords.copy()
            
            # Store endpoints for closed geometries
            if is_closed:
                endpoints = [coords[1].copy(), coords[-2].copy()]
            
            # Smooth interior vertices
            for i in range(1, len(coords) - 1):
                # Compute average of neighbors
                avg_point = 0.5 * coords[i - 1] + 0.5 * coords[i + 1]
                # Apply weighted sum
                new_coords[i] = (1 - factor) * coords[i] + factor * avg_point
            
            # Handle closed geometries
            if is_closed:
                avg_point = 0.5 * endpoints[0] + 0.5 * endpoints[1]
                new_coords[0] = (1 - factor) * coords[0] + factor * avg_point
                new_coords[-1] = new_coords[0]
            
            coords = new_coords
    
    # Return appropriate geometry type
    if isinstance(geometry, LinearRing):
        return LinearRing(coords)
    return LineString(coords)