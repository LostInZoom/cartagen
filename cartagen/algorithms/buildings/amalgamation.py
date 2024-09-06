# This file contains several algorithms to amalgamate or aggregate two or more buildings
import numpy as np
from shapely.geometry import Polygon,MultiPolygon,Point

from cartagen.utils.math.morphology import close_multipolygon, open_polygon
from cartagen.utils.math.vector import Vector2D
from cartagen.utils.geometry.segment import get_segment_list
 
def morphological_amalgamation(buildings, buffer, edge_length, threshold=0.2):
    """
    Amalgamate buildings using dilation and erosion.
    
    The amalgamation algorithm proposed by Damen *et al.*
    :footcite:p:`damen:2008` is based on morphological dilations and erosions with a square cap.
    It is particularly useful to keep the overall shape of building blocks.

    Parameters
    ----------
    buildings : list of Polygon
        Buildings to amalgamate.
    buffer : float
        Size of the buffer used for dilation (in meters).
        Buildings closer than 2 times the buffer size are amalgamated.
    edge_length : float
        Minimum length of edges in the amalgamated geometries
        (a simplification process is carried out).
    threshold : float, optional
        The threshold used for the Douglas-Peucker simplification of the
        resulting polygon before the edge removal. Do not change.

    Returns
    -------
    list of Polygon

    See Also
    --------
    boffet_areas :
        Calculate urban areas from buildings. Useful for smaller scale maps.

    References
    ----------
    .. footbibliography::
    
    Examples
    --------
    >>> buildings = [Polygon([(1, 0), (9, 0), (9, 6), (1, 6), (1, 0)]), Polygon([(10, 0), (17, 0), (17, 6), (10, 6), (10, 0)])]
    >>> morphological_amalgamation(buildings, 1.0, 1.0)
    <POLYGON ((1.207 1.983, 2.547 5.885, 16.768 4.282, 15.42 0.148, 1.207 1.983))>
    """
    output_collection = []
    clusters = []
    multipolygon = MultiPolygon(buildings)

    # make a morphological closing on the multipolygon
    closed = close_multipolygon(multipolygon, buffer, quad_segs=2)
    merged = open_polygon(closed, buffer, quad_segs=2)

    if(merged.geom_type == 'Polygon'):
        clusters.append(merged)
    elif (merged.geom_type == 'MultiPolygon'):
        for simple in merged.geoms:
            clusters.append(simple)
    
    for newbuilding in clusters:
        simplified = __edge_removal(newbuilding, edge_length, threshold)
        output_collection.append(simplified)
    
    return output_collection

def __edge_removal(polygon, edge_length, threshold):
    # first filter unnecessary vertices with a Douglas & Peucker filter with a very small threshold
    filtered = polygon.simplify(threshold)

    # then initialise the final list of vertices
    final_coords = []
    final_coords.append(filtered.exterior.coords[0])

    # get the list of segments of the polygon
    segment_list = get_segment_list(filtered)

    # then, loop on the segments to remove the ones that are too short
    previousEdge = segment_list[-1]
    for i in range(0, len(segment_list)):
        segment = segment_list[i]
        final_coords.append(segment.point2)

        # check the segment length
        if segment.length() < edge_length:
            # arrived here, the edge is too short
            nextEdge = None
            if(i == len(segment_list) - 1):
                nextEdge = segment_list[0]
            else:
                nextEdge = segment_list[i + 1]
            
            # we compute the angle between previousEdge and nextEdge
            angle = abs(previousEdge.orientation() - nextEdge.orientation())

            if(angle < np.pi/4):
                if len(final_coords) > 2:
                    # offset case
                    # create a vector from the edge
                    vector = Vector2D.from_segment(segment)
                    # keep last vertex of final_coords in memory
                    last_vertex = final_coords[-1]
                    # remove the last two vertices
                    del final_coords[-2:]
                    # # get the antepenultimate vertex of final_coords (which is now the last)
                    antepenultimate = final_coords[len(final_coords)-1]
                    # translate this vertex with the vector
                    translated = vector.translate(Point(antepenultimate))
                    # add translated and then lastVertex to the list of vertices
                    final_coords.append(translated.coords[0])
                    final_coords.append(last_vertex)

            elif(angle < 3*np.pi/4):
                # it is a corner case
                # get the intersection point of previousEdge and nextEdge considered as straight lines
                intersection = previousEdge.straight_line_intersection(nextEdge)
                # keep last vertex of final_coords in memory
                last_vertex = final_coords[-1]
                # remove last vertex
                del final_coords[-1:]
                # add intersection and then lastVertex to the list of vertices
                final_coords.append(intersection.coords[0])
                final_coords.append(last_vertex)
            else:
                if len(final_coords) > 2:
                    # intrusion or protrusion case
                    # create a vector from the edge
                    vector = Vector2D.from_segment(segment)
                    # remove the last two vertices
                    del final_coords[-2:]
                    # get the antepenultimate vertex of final_coords (which is now the last)
                    antepenultimate = final_coords[-1]
                    # translate this vertex with the vector
                    translated = vector.translate(Point(antepenultimate))
                    final_coords.append(translated.coords[0])
                    final_coords.append(nextEdge.point2)

        previousEdge = segment

    return Polygon(final_coords)

