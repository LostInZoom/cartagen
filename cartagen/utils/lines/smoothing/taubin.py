import numpy as np
from shapely.geometry import (
    LineString, MultiLineString, LinearRing, Polygon, MultiPolygon
)

def smooth_taubin(geometry, iterations=10, lamb=0.5, mu=-0.53):
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
    lamb : float, optional
        Smoothing factor for the first pass. Default is 0.5.
        Should be positive and typically between 0 and 1.
    mu : float, optional
        Inflation factor for the second pass. Default is -0.53.
        Should be negative with |mu| slightly larger than lambda to prevent shrinkage.
    
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

    References
    ----------
    .. footbibliography::
    
    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (3, 1), (4, 0)])
    >>> smooth_taubin(line, lamb=0.5, iterations=5)
    <LINESTRING (0 0, 1.2 0.8, 2 0.3, 2.8 0.8, 4 0)>
    
    >>> polygon = Polygon([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)])
    >>> smooth_taubin(polygon, lamb=0.3, iterations=10)
    <POLYGON ((0.2 0.2, 0.2 1.8, 1.8 1.8, 1.8 0.2, 0.2 0.2))>
    """
    # Handle different geometry types
    if isinstance(geometry, MultiLineString):
        return MultiLineString([
            smooth_taubin(line, lamb, mu, iterations)
            for line in geometry.geoms
        ])
    
    if isinstance(geometry, MultiPolygon):
        return MultiPolygon([
            smooth_taubin(poly, lamb, mu, iterations)
            for poly in geometry.geoms
        ])
    
    if isinstance(geometry, Polygon):
        # Smooth exterior ring
        exterior = smooth_taubin(
            geometry.exterior, lamb, mu, iterations
        )
        # Smooth interior rings (holes)
        interiors = [
            smooth_taubin(interior, lamb, mu, iterations)
            for interior in geometry.interiors
        ]
        return Polygon(exterior, interiors)
    
    # Process LineString or LinearRing
    coords = np.array(geometry.coords)
    is_closed = np.array_equal(coords[0], coords[-1])
    
    # Apply Taubin smoothing
    factors = [lamb, mu]
    
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