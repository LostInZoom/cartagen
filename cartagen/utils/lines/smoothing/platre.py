import numpy as np
from shapely import LineString, MultiLineString

def smooth_platre(line, sigma=2, curvature=0.05):
    """
    Smooth a line and preserve the integrity of sharp turns.

    The PLATRE algorithm was created by Emmanuel Fritsch :footcite:p:`fritsch:1998` to
    attenuate minor bends in a polyline while preserving the position and integrity of the sharpest turns.
    Unlike simple coordinate-based filters, it operates on the line's intrinsic 
    geometry (curvature/angle) to maintain structural consistency.
    
    Parameters
    ----------
    line : LineString or MultiLineString
        The line to smooth.
    sigma : float, optional
        Strength of the smoothing.
    curvature : float, optional
        The curvature threshold to identify sharpest turns.
        Low values is closer to the original line.

    Returns
    -------
    LineString, MultiLineString

    See Also
    --------
    smooth_catmull_rom :
        Smooth a line or polygon and preserve vertexes.
    smooth_chaikin :
        Smooth a line or polygon by cutting corners.
    smooth_gaussian :
        Smooth a line or a polygon and attenuate its inflexion points.
    smooth_taubin :
        Smooth a line or polygon and prevent shrinkage.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (5, 3)])
    >>> c4.platre(line, sigma=2, curvature=5)
    <LINESTRING (0 0, 0.9950928016739138 0.6091961570208742, 2.0147215949782598 1.172411528937378, 5 3)>
    """

    if line.geom_type not in ['LineString', 'MultiLineString', 'LinearRing']:
        raise ValueError(f'{line.geom_type} geometry type cannot be simplified.')
    
    if line.geom_type == 'MultiLineString':
        geoms = [platre(geom, sigma, curvature) for geom in line.geoms]
        return MultiLineString(geoms)
    
    # --- Préparation des données ---
    coords = np.array(line.coords)
    x, y = coords[:, 0], coords[:, 1]
    
    # Calcul de l'abscisse curviligne (distance cumulée)
    dx = np.diff(x)
    dy = np.diff(y)
    ds = np.sqrt(dx**2 + dy**2)
    s = np.concatenate(([0], np.cumsum(ds)))
    
    # --- 1. Changement de représentation ---
    # Calcul des angles des segments
    angles = np.arctan2(dy, dx)
    # Déroulage des angles pour éviter les sauts de 2pi
    angles_unwrapped = np.unwrap(angles)
    
    # --- 2. Lissage de la courbure ---
    # Note : Lisser l'angle revient à intégrer le lissage de la courbure.
    # On crée un noyau Gaussien manuel
    size = int(6 * sigma + 1)
    if size % 2 == 0: size += 1
    t = np.arange(size) - (size - 1) / 2
    kernel = np.exp(-0.5 * (t / sigma)**2)
    kernel /= kernel.sum()
    
    # Application du lissage (padding pour garder la même longueur)
    padded_angles = np.pad(angles_unwrapped, (size//2, size//2), mode='edge')
    smooth_angles = np.convolve(padded_angles, kernel, mode='valid')
    
    # --- 3. Retour à la représentation en polyligne ---
    # On reconstruit les segments avec les nouveaux angles mais les distances originales
    new_dx = ds * np.cos(smooth_angles)
    new_dy = ds * np.sin(smooth_angles)
    
    x_smooth = np.concatenate(([x[0]], x[0] + np.cumsum(new_dx)))
    y_smooth = np.concatenate(([y[0]], y[0] + np.cumsum(new_dy)))
    coords_smooth = np.column_stack((x_smooth, y_smooth))
    
    # --- 4. Repositionnement des virages serrés ---
    # Calcul de la courbure locale (différence d'angle / distance moyenne)
    # On utilise les angles originaux pour détecter les virages
    curv = np.zeros_like(s)
    curv[1:-1] = np.abs(np.diff(angles_unwrapped)) / (0.5 * (ds[:-1] + ds[1:]))
    
    # Identification des points à fixer (extrémités + courbure > seuil)
    fixed_indices = np.where(curv > curvature)[0]
    fixed_indices = np.unique(np.concatenate(([0], fixed_indices, [len(s)-1])))
    
    # Calcul des vecteurs de déplacement aux points fixes
    # (Position Initiale - Position Lissée)
    displacements = coords[fixed_indices] - coords_smooth[fixed_indices]
    
    # Interpolation linéaire des déplacements pour tous les points
    # On interpole chaque composante (x, y) en fonction de l'abscisse curviligne
    interp_dx = np.interp(s, s[fixed_indices], displacements[:, 0])
    interp_dy = np.interp(s, s[fixed_indices], displacements[:, 1])
    
    final_coords = coords_smooth + np.column_stack((interp_dx, interp_dy))
    
    return LineString(final_coords)