# this file contains several line simplification algorithms, but not the Douglas & Peucker algorithm because it is already implemented in shapely

from shapely.geometry import LineString, Point, Polygon, MultiPoint
from cartagen4py.utils.tessellation.hexagonal import *
from shapely.ops import nearest_points, transform

# Visvalingam-Whyatt algorithm (1993)
def visvalingam_whyatt(line, area_tolerance):
    final_coords = list(line.coords)
    while (True):
        pt_area_min = None
        area_min = line.envelope.area
        current_line = LineString(tuple(final_coords))
        for i in range(1,len(final_coords)-2):
            pt = final_coords[i]
            if(__contains_another_point(current_line, pt, i)):
                continue
            area = __compute_area_point(current_line, pt, i)
            if (area < area_min):
                area_min = area
                if(area < area_tolerance):
                    pt_area_min = pt
        if(pt_area_min is None):
            break
        # remove pt from the line
        final_coords.remove(pt_area_min)

    return LineString(tuple(final_coords))


def __compute_area_point(line, pt, index):
    first = line.coords[index-1]
    last = line.coords[index+1]
    triangle = [first, pt, last, first]
    polygon = Polygon(triangle)
    return polygon.area

def __contains_another_point(line, pt, index):
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

# Raposo simplification algorithm (2010): uses an hexagonal tessallation, with a size related to the final scale, 
# and the algorithm only retains one vertex per hexagonal cell. Be careful, it uses the scale as parameter. If the centroid parameter 
# is True, the vertices inside an hexagon cell are replaced by the centroid of the cell; if it is False, they are replaced by the nearest
# vertex to the centroid of the cell.
def raposo_simplification(line, initial_scale, final_scale, centroid=True, tobler=False):
    width = 0
    if tobler:
        width = __compute_tobler_width(final_scale)
    else:
        width = __compute_raposo_width(line, initial_scale, final_scale)
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

def __compute_raposo_width(line,initial_scale,final_scale):
    firstFactor = line.length / (len(line.coords)-1)
    secondFactor = final_scale / initial_scale
    return firstFactor * secondFactor

# Based on Tobler's rule-of-thumb (1987)
def __compute_tobler_width(final_scale):
    return final_scale * 5 / 4000 

# Li-Openshaw simplification algorithm (1993). The simplification factor is the size of the regular grid applied on the line
def __li_openshaw_simplification(line, cell_size):
    # TODO
    return line
