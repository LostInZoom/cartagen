# This file contains several algorithms to amalgamate or aggregate two or more buildings

import math

from shapely.geometry import Polygon,MultiPolygon,Point

from cartagen4py.utils.math.morphology import closing_multi_polygon, opening
from cartagen4py.utils.math.vector import Vector2D
from cartagen4py.utils.geometry.segment import get_segment_list

# This the amalgamation algorithm from Damen et al. 2008 (https://www.semanticscholar.org/paper/High-Quality-Building-Generalization-by-Extending-Damen-Kreveld/b64618584b3ae3725da7eeb5a545d1580e5f2113)
def morphological_amalgamation(buildings, buffer_size, edge_length):
    output_collection = []
    clusters = []
    multipolygon = MultiPolygon(buildings)

    # make a morphological closing on the multipolygon
    closed = closing_multi_polygon(multipolygon, buffer_size, cap_style=2)
    print(closed)
    merged = opening(closed, buffer_size, cap_style=2)

    if(merged.geom_type == 'Polygon'):
        clusters.append(merged)
    elif (merged.geom_type == 'MultiPolygon'):
        for simple in merged.geoms:
            clusters.append(simple)
    
    for newbuilding in clusters:
        simplified = __edge_removal(newbuilding, edge_length)
        output_collection.append(simplified)
    
    return output_collection

def __edge_removal(polygon, edge_length):
    # first filter unnecessary vertices with a Douglas & Peucker filter with a very small threshold
    filtered = polygon.simplify(0.2)

    # then initialise the final list of vertices and add the first vertex
    final_coords = []
    final_coords.append(filtered.exterior.coords[0])

    # get the list of segments of the polygon
    segment_list = get_segment_list(polygon)

    # then, loop on the segments to remove the ones that are too short
    previousEdge = segment_list[len(segment_list)-1]
    for i in range(0,len(segment_list)-1):
        segment = segment_list[i]
        final_coords.append(segment.point2)

        # check the segment length
        if(segment.length() >= edge_length):
            continue

        # arrived here, the edge is too short
        nextEdge = None
        if(i == len(segment_list)-1):
            nextEdge = segment_list[0]
        else:
            nextEdge = segment_list[i+1]
        
        # we compute the angle between previousEdge and nextEdge
        angle = abs(previousEdge.orientation() - nextEdge.orientation())
        if(angle < math.pi/4):
            # offset case
            # create a vector from the edge
            vector = Vector2D(nextEdge)
            # keep last vertex of final_coords in memory
            last_vertex = final_coords[len(final_coords)-1]
            # remove the last two vertices
            l_element = len(tuple)-2
            final_coords = final_coords[:l_element]
            # get the antepenultimate vertex of final_coords (which is now the last)
            antepenultimate = final_coords[len(final_coords)-1]
            # translate this vertex with the vector
            translated = vector.translate(Point(antepenultimate))
            # add translated and then lastVertex to the list of vertices
            final_coords.append(translated.coords)
            final_coords.append(last_vertex)
        elif(angle < 3*math.pi/4):
            # it is a corner case
            # get the intersection point of previousEdge and nextEdge considered as straight lines
            intersection = previousEdge.straight_line_intersection(nextEdge)
            # keep last vertex of final_coords in memory
            last_vertex = final_coords[len(final_coords)-1]
            # remove last vertex
            # TODO
            # add intersection and then lastVertex to the list of vertices
            final_coords.append(intersection.coords)
            final_coords.append(last_vertex)
        else:
            # intrusion or protrusion case
            # create a vector from the edge
            vector = Vector2D(nextEdge)
            # remove the last two vertices
            l_element = len(tuple)-2
            final_coords = final_coords[:l_element]
            # get the antepenultimate vertex of final_coords (which is now the last)
            antepenultimate = final_coords[len(final_coords)-1]
            # translate this vertex with the vector
            translated = vector.translate(Point(antepenultimate))
            final_coords.append(intersection.coords)
            final_coords.append(nextEdge.point2)
    
    return Polygon(final_coords)

