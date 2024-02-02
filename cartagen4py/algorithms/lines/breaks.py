import geopandas as gpd
import shapely, numpy

from cartagen4py.utils.geometry.angle import *
from cartagen4py.utils.geometry.dilation import *

def max_break(line, width):
    """
    
    """

    def __get_bend_side(coords):
        total = 0
        start = shapely.Point(coords[0])

        for i in range(1, len(coords) - 1):
            c1 = coords[i]
            c2 = coords[i + 1]

            p1 = shapely.Point(c1)
            p2 = shapely.Point(c2)

            angle = angle_3_pts(start, p2, p1)

            total += angle

        if total > 0:
            return 'left'
        else:
            return 'right'

    side = __get_bend_side(list(line.coords))

    offset = width / 2

    if side == 'left':
        offset = - offset

    dilated = offset_curve(line, offset, cap_style='flat', quad_segs=8)
    
    return shapely.LineString(dilated[0])

    