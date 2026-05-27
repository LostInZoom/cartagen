import warnings
import numpy as np
from shapely.geometry import (
    LineString, MultiLineString, Polygon, MultiPolygon,
    LinearRing, GeometryCollection,
)

def smooth_snake(geometry, iterations=10, alpha=0.1, beta=0.1, gamma=1.0, kappa=0.1):
    """
    Smooth a line or polygon using the snake method.

    This algorithm proposed by Burghardt :footcite:p:`burghardt:2005`
    treats the geometry as an elastic snake optimizing an energy
    functional. Internal energy terms (``alpha`` and ``beta``) control the tension
    and rigidity of the line to smooth out sharp noise, while an external spring
    force (``kappa``) anchors the vertices to their original positions to prevent
    geometric shrinkage and uncontrolled displacement (Burghardt's shape-preservation
    term).

    The smoothed position is obtained by solving the linear system

    .. math::

        (A + (\\gamma + \\kappa)\\,I)\\,\\mathbf{x}^* =
            \\gamma\\,\\mathbf{x}_0 + \\kappa\\,\\mathbf{x}_0

    Because both the system matrix and the right-hand side are constant across
    iterations, the snake converges to the same fixed point regardless of the
    number of iterations.  The ``iterations`` parameter therefore controls a
    **pseudo-temporal relaxation**: starting from ``x_0``, the snake is advanced
    ``iterations`` steps of implicit Euler integration, which can be useful to
    observe intermediate states or to limit over-smoothing on short geometries.
    With ``iterations=1`` the result is the direct closed-form solution.

    Multi-geometries are handled recursively.  For polygons, the smoothing is
    applied independently to the exterior ring and to each interior ring (hole).

    Parameters
    ----------
    geometry : LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
        The geometry to smooth.
        If an open line is provided, the endpoints are preserved.
        If a closed ring or polygon is provided, the smoothing wraps around.
    iterations : int, optional
        Number of implicit Euler steps for the snake evolution.  With
        ``iterations=1`` the direct closed-form solution is returned.
        Default is 10.
    alpha : float, optional
        Tension parameter (continuity).  Higher values penalise line stretching
        and pull vertices toward their neighbours.  Default is 0.1.
    beta : float, optional
        Rigidity parameter (curvature).  Higher values penalise sharp turns and
        produce smoother, more rounded curves.  Default is 0.1.
    gamma : float, optional
        Step-size / weighting of the current position in the RHS.  Acts as a
        regularisation weight on the current iterate rather than a physical
        viscosity coefficient.  Higher values slow down the evolution per step.
        Default is 1.0.
    kappa : float, optional
        Shape-preservation parameter (Burghardt's external spring force).  Higher
        values pull the snake back to its original shape, preventing shrinkage
        and structural collapse.  Setting ``kappa=0`` disables this term and
        closed loops will shrink toward their centroid over many iterations.
        Default is 0.1.

    Returns
    -------
    LineString, Polygon, MultiLineString, MultiPolygon, LinearRing
        Smoothed geometry of the same type as the input.

    Notes
    -----
    The stiffness matrix ``A`` is pentadiagonal: the tension band uses
    coefficients ``[-α, 2α, -α]`` and the rigidity band uses
    ``[β, -4β, 6β, -4β, β]``.  For closed rings the matrix is circulant
    (wrap-around indices); for open lines the rows are clamped at the
    boundaries and the first/last rows are replaced by identity rows to
    enforce exact endpoint preservation.

    The evolution at each step solves::

        M · x_new = γ · x_current + κ · x_original

    where ``M = A + (γ + κ) · I``.  This is equivalent to one step of
    implicit gradient descent on the snake energy.

    See Also
    --------
    smooth_gaussian :
        Smooth a line or a polygon and attenuate its inflexion points.
    smooth_platre :
        Smooth a line and preserve the integrity of sharp turns.
    smooth_taubin :
        Smooth a line or polygon and prevent shrinkage.
    smooth_topographic :
        Smooth a line or polygon and mimic hand-made cartographic generalization.
    smooth_wma :
        Smooth a line or polygon using a low-pass filter.

    Examples
    --------
    >>> line = LineString([(0, 0), (1, 1), (2, 0), (3, 1), (4, 0)])
    >>> smooth_snake(line)
    <LINESTRING (0 0, 0.875 0.256, 2 0.298, 3.125 0.256, 4 0)>
    """
    if geometry.is_empty or geometry.geom_type in ("Point", "MultiPoint"):
        return geometry

    # ------------------------------------------------------------------
    # Recursive handling of multi-part and polygon geometries
    # ------------------------------------------------------------------
    if geometry.geom_type in ("MultiLineString", "MultiPolygon", "GeometryCollection"):
        parts = [
            smooth_snake(part, iterations, alpha, beta, gamma, kappa)
            for part in geometry.geoms
        ]
        if geometry.geom_type == "MultiLineString":
            return MultiLineString(parts)
        if geometry.geom_type == "MultiPolygon":
            return MultiPolygon(parts)
        return GeometryCollection(parts)

    if geometry.geom_type == "Polygon":
        ext = smooth_snake(geometry.exterior, iterations, alpha, beta, gamma, kappa)
        if len(ext.coords) < 4:
            warnings.warn(
                "Exterior ring degenerated after smoothing; "
                "returning original geometry.",
                UserWarning,
                stacklevel=2,
            )
            return geometry
        holes = []
        for ring in geometry.interiors:
            smoothed = smooth_snake(ring, iterations, alpha, beta, gamma, kappa)
            if len(smoothed.coords) >= 4:
                holes.append(smoothed)
            else:
                warnings.warn(
                    "An interior ring degenerated after smoothing and was dropped.",
                    UserWarning,
                    stacklevel=2,
                )
        poly = Polygon(ext, holes)
        return poly if poly.is_valid else poly.buffer(0)

    # ------------------------------------------------------------------
    # Core logic: LineString or LinearRing
    # ------------------------------------------------------------------
    coords = np.array(geometry.coords)

    is_closed = geometry.geom_type == "LinearRing" or (
        len(coords) > 1 and np.allclose(coords[0], coords[-1])
    )

    if len(coords) < 3:
        return geometry

    # Work with the open vertex list for closed rings (strip the repeated vertex)
    if is_closed:
        coords = coords[:-1]

    n = len(coords)
    if n < 3:
        # Cannot build a meaningful stiffness matrix; reconstruct and return
        if is_closed:
            closed = np.vstack([coords, coords[0]])
            return LinearRing(closed) if geometry.geom_type == "LinearRing" else LineString(closed)
        return LineString(coords)

    # Preserve endpoints for open lines (used both in the matrix and the RHS)
    start_pt = coords[0].copy()
    end_pt   = coords[-1].copy()

    # ------------------------------------------------------------------
    # Build the stiffness matrix A (pentadiagonal)
    # ------------------------------------------------------------------
    A = np.zeros((n, n))

    if is_closed:
        # Circulant matrix: wrap-around indices via modular arithmetic
        for i in range(n):
            # Tension  (1st-derivative penalty)
            A[i, (i - 1) % n] += -alpha
            A[i,  i          ] +=  2 * alpha
            A[i, (i + 1) % n] += -alpha

            # Rigidity (2nd-derivative penalty)
            A[i, (i - 2) % n] +=  beta
            A[i, (i - 1) % n] += -4 * beta
            A[i,  i          ] +=  6 * beta
            A[i, (i + 1) % n] += -4 * beta
            A[i, (i + 2) % n] +=  beta
    else:
        # Open line: boundary rows will be replaced by identity rows in M,
        # so only interior rows matter here.  We still build a full matrix
        # with clamped indices to keep the code uniform, but rows 0 and n-1
        # are overwritten below.
        for i in range(1, n - 1):
            # Tension
            A[i, i - 1] += -alpha
            A[i, i    ] +=  2 * alpha
            A[i, i + 1] += -alpha

            # Rigidity (clamped at boundaries)
            A[i, max(0,     i - 2)] +=  beta
            A[i,            i - 1 ] += -4 * beta
            A[i,            i     ] +=  6 * beta
            A[i, min(n - 1, i + 1)] += -4 * beta
            A[i, min(n - 1, i + 2)] +=  beta

    # ------------------------------------------------------------------
    # Build the system matrix M = A + (γ + κ) I
    # For open lines, replace boundary rows with identity to clamp endpoints.
    # ------------------------------------------------------------------
    M = A + (gamma + kappa) * np.eye(n)

    if not is_closed:
        M[0,  :] = 0.0;  M[0,  0 ] = 1.0
        M[-1, :] = 0.0;  M[-1, -1] = 1.0

    # Pre-factorise the (constant) system matrix once — the RHS changes
    # between iterations but M does not.
    M_lu = np.linalg.inv(M)   # For small n this is fine; for large n use scipy LU

    # Original positions, kept fixed throughout (shape-preservation anchor)
    original = coords.copy()
    current  = coords.copy()

    # ------------------------------------------------------------------
    # Iterative pseudo-temporal evolution
    # Each step: current ← M⁻¹ · (γ · current + κ · original)
    # The fixed point of this iteration is the closed-form solution.
    # ------------------------------------------------------------------
    for _ in range(iterations):
        rhs = gamma * current + kappa * original

        if not is_closed:
            rhs[0 ] = start_pt
            rhs[-1] = end_pt

        current = M_lu @ rhs   # shape (n, 2) — X and Y solved simultaneously

    # ------------------------------------------------------------------
    # Reconstruct the geometry
    # ------------------------------------------------------------------
    if is_closed:
        result = np.vstack([current, current[0]])
        if geometry.geom_type == "LinearRing":
            return LinearRing(result) if len(result) >= 4 else geometry
        return LineString(result)

    return LineString(current)