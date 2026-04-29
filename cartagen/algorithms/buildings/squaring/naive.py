import numpy as np
from shapely.geometry import Polygon, Point

from cartagen.utils.geometry.polygon import orientation
from cartagen.utils.math.vector import Vector2D

def square_polygon_naive(polygon, orient='primary', angle_tolerance=8.0, correct_tolerance=0.6, remove_flat=True):
    """
    Squares a polygon according to its orientation.

    This method, described in Touya, :footcite:p:`touya:2016` first
    calculates the orientation of the polygon. Sides are then
    corrected depending on the angles formed at the vertexes and
    on their alignment regarding the calculated orientation.

    Parameters
    ----------
    polygon : Polygon
        The polygon to square.
    orient : str, optional
        The method to calculate the orientation. Be aware that the orientation
        of the polygon defines how the sides are corrected.

        - **'primary'** calculates the orientation of the
          longest side of the provided polygon.
        - **'mbr'** calculates the orientation
          of the long side of the minimum rotated bounding rectangle.
        - **'mbtr'** calculates the orientation
          of the long side of the minimum rotated bounding touching rectangle.
          It is the same as the mbr but the rectangle and the polygon
          must have at least one side in common.
        - **'swo'** or statistical weighted orientation described in
          Duchêne, :footcite:p:`duchene:2003` calculates
          the orientation of a polygon using the statistical weighted orientation.
          This method relies on the length and orientation of the longest and
          second longest segment between two vertexes of the polygon.

    angle_tolerance : float, optional
        Tolerance in degrees to square the considered angle.
    correct_tolerance : float, optional
        Tolerance in degrees to consider the angle to be already flat or right.
    remove_flat : bool, optional
        If set to True, vertexes with an angle detected or corrected as flat
        are removed. Thus, the resulting polygon can have less vertexes than
        the provided one.

    Returns
    -------
    Polygon

    See Also
    --------
    square_polygon_ls :
        Squares a polygon using the method of least squares.
    orientation :
        Calculate the orientation of a polygon.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 1), (1.1, 1), (1, 0)])
    >>> square_polygon_naive(polygon)
    """
    angle_tolerance = angle_tolerance * np.pi / 180
    correct_tolerance = correct_tolerance * np.pi / 180

    # Calculate the orientation using SWO
    o = orientation(polygon, orient)
    vo = Vector2D.from_angle(o, 1)

    points = list(polygon.exterior.coords)[:-1]
    nb_edges = len(points) - 1

    vecs, angles, align = [None] * nb_edges, [None] * nb_edges, [None] * nb_edges

    primary_axe = 0
    norm_max = 0

    for i in range(0, nb_edges):
        vector = Vector2D.from_points(Point(points[i]), Point(points[(i + 1) % nb_edges]))
        vecs[i] = vector

        norm = vector.get_norm()
        if norm > norm_max:
            norm_max = norm
            primary_axe = i

    def __get_vecs_around(i):
        v1 = (len(vecs) - 1) if i == 0 else (i - 1)
        v2 = i
        return v1, v2
    
    def __get_point_to_move(vertice, vect_to_move):
        if vect_to_move == vertice:
            return (vertice + 1) % nb_edges
        if vertice == 0:
            return nb_edges - 1
        return vertice - 1
    
    def __update():
        for i in range(0, nb_edges):
            vecs[i] = Vector2D.from_points(Point(points[i]), Point(points[(i + 1) % nb_edges]))
        angles[0] = np.pi - vecs[-1].angle(vecs[0])
        for i in range(0, len(vecs) - 1):
            angles[i + 1] = np.pi - vecs[i].angle(vecs[i + 1])
            align[i] = abs(vecs[i].product(vo))
        align[len(align) - 1] = abs(vecs[len(align) - 1].product(vo))
    
    def __signed_angle(v1, v2):
        return np.arctan2(v2.y, v2.x) - np.arctan2(v1.y, v1.x)

    def __rotate_vec(deg, vec_to_move, vec_fixed, i):
        vsquared = None
        if deg == 90:
            vsquared = Vector2D.from_point(Point(-vecs[vec_fixed].y, vecs[vec_fixed].x))
            vsquared.normalize()
            vsquared = vsquared.const_product(vecs[vec_to_move].get_norm())
            if vec_to_move == i:
                if vsquared.scalar_product(vecs[vec_to_move]) < 0:
                    vsquared = vsquared.const_product(-1)
            else:
                if vsquared.scalar_product(vecs[vec_to_move]) > 0:
                    vsquared = vsquared.const_product(-1)
            return vsquared
        elif deg == 45:
            vp = vecs[vec_fixed]
            v1, v2 = vec_to_move, vec_fixed
            if vec_to_move == i:
                vp = vecs[vec_fixed].const_product(-1)
                v1, v2 = vec_fixed, vec_to_move
                if __signed_angle(vecs[v1], vecs[v2]) < 0:
                    x = 0.7071067811865476 * (vp.x - vp.y)
                    y = 0.7071067811865476 * (vp.x + vp.y)
                    vsquared = Vector2D.from_point(Point(x, y))
                else:
                    x = 0.7071067811865476 * (vp.x + vp.y)
                    y = 0.7071067811865476 * (-1 * vp.x + vp.y)
                    vsquared = Vector2D.from_point(Point(x, y))
            else:
                if __signed_angle(vecs[v1], vecs[v2]) > 0:
                    x = 0.7071067811865476 * (vp.x - vp.y)
                    y = 0.7071067811865476 * (vp.x + vp.y)
                    vsquared = Vector2D.from_point(Point(x, y))
                else:
                    x = 0.7071067811865476 * (vp.x + vp.y)
                    y = 0.7071067811865476 * (-1 * vp.x + vp.y)
                    vsquared = Vector2D.from_point(Point(x, y))
            return vsquared
        elif deg == 0:
            vp = vecs[vec_fixed]
            if vec_to_move != i:
                vpn = vecs[vec_fixed].const_product(-1)
                vp = Vector2D.from_point(vpn.x, vpn.y).normalize()

            vpn = vp.const_product(-1)
            vp = Vector2D.from_point(vpn.x, vpn.y).normalize().const_product(vecs[vec_to_move].norm())
            return vp
        return vsquared
    
    angles[0] = np.pi - vecs[-1].angle(vecs[0])
    for i in range(0, len(vecs) - 1):
        angles[i + 1] = np.pi - vecs[i].angle(vecs[i + 1])
        align[i] = abs(vecs[i].product(vo))
    align[len(align) - 1] = abs(vecs[len(align) - 1].product(vo))

    for i in range(0, len(angles)):
        if abs(np.pi - angles[i]) <= angle_tolerance and abs(np.pi - angles[i]) > correct_tolerance:
            v = __get_vecs_around(i)
            vv = vecs[v[0]].add(vecs[v[1]])
            vv.normalize()
            norm = vecs[v[0]].get_norm() * np.cos(vecs[v[0]].angle(vv))
            vsquared = vv.const_product(norm)
            points[i] = (points[v[0]][0] + vsquared.x, points[v[0]][1] + vsquared.y)
            __update()

    for i in range(0, len(angles)):
        if abs(np.pi/2 - angles[i]) <= angle_tolerance and abs(np.pi/2 - angles[i]) > correct_tolerance:
            v = __get_vecs_around(i)
            vec_to_move, vec_fixed = v[0], v[1]
            if align[vec_to_move] < align[vec_fixed]:
                vec_to_move, vec_fixed = v[1], v[0]
            
            point_to_move = __get_point_to_move(i, vec_to_move)
            angle_before = v[0]

            angle_prec_right = True if abs(angles[angle_before] - np.pi / 2) < correct_tolerance else False
            if angle_prec_right and vec_to_move == point_to_move:
                continue

            vsquared = __rotate_vec(90, vec_to_move, vec_fixed, i)
            points[point_to_move] = (points[i][0] + vsquared.x, points[i][1] + vsquared.y)
            __update()
        
        elif abs(np.pi/4 - angles[i]) <= angle_tolerance and abs(np.pi/4 - angles[i]) > 0.01:
            v = __get_vecs_around(i)
            vec_to_move, vec_fixed = v[0], v[1]
            if align[vec_to_move] < align[vec_fixed]:
                vec_to_move, vec_fixed = v[1], v[0]
            
            point_to_move = __get_point_to_move(i, vec_to_move)
            angle_before = v[0]
            angle_after = (i + 1) % nb_edges

            angle_before_isright = abs(angles[angle_before] - angles[angle_before]) < correct_tolerance
            angle_after_isright = abs(angles[angle_after] - angles[angle_after]) < correct_tolerance

            if angle_before_isright and angle_after_isright:
                continue

            vsquared = __rotate_vec(45, vec_to_move, vec_fixed, i)
            points[point_to_move] = (points[i][0] + vsquared.x, points[i][1] + vsquared.y)
            __update()
    
    result = []
    removed = []
    if remove_flat:
        for i in range(0, len(points)):
            angle = angles[i % nb_edges]
            if abs(np.pi - angle) > correct_tolerance:
                result.append(points[i])
            else:
                removed.append(points[i])
    else:
        result = points

    return Polygon(result)