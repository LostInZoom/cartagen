# This file contains several algorithms to amalgamate or aggregate two or more buildings
import numpy as np
from shapely.geometry import Polygon,MultiPolygon,Point,LineString
from shapely.ops import unary_union
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

class BuildingAggregator:
    def __init__(self, poly_a, poly_b):
        self.poly_a = poly_a
        self.poly_b = poly_b

    def _get_adjacent_edges(self, poly, point):
        """Find the two edges of a polygon that share the given vertex."""
        coords = list(poly.exterior.coords)
        for i, p in enumerate(coords[:-1]):
            if Point(p).equals(point):
                # Arêtes précédente et suivante
                # Previous and next edges
                e1 = LineString([coords[i-1], p])
                e2 = LineString([p, coords[i+1]])
                return e1, e2
        return None, None

    def _extend_line(self, line, length=100):
        """Extend a segment to enable the intersections."""
        coords = list(line.coords)
        p1, p2 = np.array(coords[0]), np.array(coords[1])
        vec = p2 - p1
        unit_vec = vec / np.linalg.norm(vec)
        # On prolonge dans les deux sens pour être sûr de croiser l'autre segment
        # We extend in both directions to ensure intersection with the other segment
        new_p1 = p1 - unit_vec * length
        new_p2 = p2 + unit_vec * length
        return LineString([new_p1, new_p2])

    def aggregate(self):
        # --- ÉTAPE 1 : Convex Hull ---
        combined = unary_union([self.poly_a, self.poly_b])
        hull = combined.convex_hull
        
        # --- ÉTAPE 2 : Search for bridges between polygons (a1b1, a2b2) ---
        # On cherche les segments du hull qui ne sont pas sur les polygones d'origine
        # We look for hull segments that are not part of the original polygons
        hull_coords = list(hull.exterior.coords)
        bridge_segments = []
        for i in range(len(hull_coords)-1):
            seg = LineString([hull_coords[i], hull_coords[i+1]])
            if not (seg.within(self.poly_a.buffer(1e-7)) or seg.within(self.poly_b.buffer(1e-7))):
                bridge_segments.append(seg)
        
        if len(bridge_segments) < 2:
            return combined # no clear bridges found, return the convex hull

        # --- ÉTAPE 3 : Identification of the Internal/External links ---
        solutions = []
        # On traite les deux ponts principaux (a1-b1 et a2-b2)
        # For the two main bridges (a1-b1 and a2-b2)
        for bridge in bridge_segments[:2]:
            p_a = Point(bridge.coords[0]) if Point(bridge.coords[0]).distance(self.poly_a) < 1e-7 else Point(bridge.coords[1])
            p_b = Point(bridge.coords[0]) if p_a.equals(Point(bridge.coords[1])) else Point(bridge.coords[1])
            
            # Pour chaque point, identifier arête interne (dans le hull) et externe (sur le hull)
            # For each point, identify internal edge (inside the hull) and external edge (on the hull)
            edges_a = self._get_adjacent_edges(self.poly_a, p_a)
            edges_b = self._get_adjacent_edges(self.poly_b, p_b)
            
            # Une arête est 'externe' si elle est colinéaire à un segment du hull
            # the edge is 'external' if it lies on the hull boundary
            def classify(edge):
                return "ext" if any(edge.within(hull.exterior.buffer(1e-7)) for _ in [0]) else "int"

            # Store the four possible extensions
            potential_lines = {
                "a_int": self._extend_line(edges_a[0] if classify(edges_a[0]) == "int" else edges_a[1]),
                "a_ext": self._extend_line(edges_a[1] if classify(edges_a[0]) == "int" else edges_a[0]),
                "b_int": self._extend_line(edges_b[0] if classify(edges_b[0]) == "int" else edges_b[1]),
                "b_ext": self._extend_line(edges_b[1] if classify(edges_b[0]) == "int" else edges_b[0])
            }
            solutions.append(potential_lines)

        # --- ÉTAPE 4 : Choice of the best solution to join the polygons  ---
        # Nicolas Regnauld suggests that we should try first the internal links
        # then try the other combinations if the first one is not valid 
        try:
            # Intersection of the internal links for the first bridge
            inter1 = solutions[0]["a_int"].intersection(solutions[0]["b_int"])
            # Intersection of the internal links for the second bridge
            inter2 = solutions[1]["a_int"].intersection(solutions[1]["b_int"])
            
            # Create the filling polygon
            filling_coords = [
                bridge_segments[0].coords[0], inter1.coords[0], bridge_segments[0].coords[1],
                bridge_segments[1].coords[0], inter2.coords[0], bridge_segments[1].coords[1]
            ]
            filling_poly = Polygon(filling_coords)
            final_result = unary_union([self.poly_a, self.poly_b, filling_poly])
            
            good =True
            if final_result.is_valid == False:
                good = False
            # Control the area ratio (1/2 threshold) 
            ratio = (final_result.area - (self.poly_a.area + self.poly_b.area)) / (self.poly_a.area + self.poly_b.area)
            if ratio > 0.5:
                print(f"Ratio d'augmentation ({ratio:.2f}) dépasse le seuil de 0.5 ")
                good = False
            
            # Case when the default solution is not valid; we have to try the three other solutions and select the best one
            if good == False:
                best_polygon = unary_union([self.poly_a, self.poly_b, hull])

                # Intersection des prolongements externes pour le premier pont
                inter3 = solutions[0]["a_ext"].intersection(solutions[0]["b_ext"])
                # Intersection des prolongements externes pour le second pont
                inter4 = solutions[1]["a_ext"].intersection(solutions[1]["b_ext"])

                #********************  Ext-Ext  ********************#
                # create filling polygon for ext-ext
                filling_coords_ext_ext = [
                    bridge_segments[0].coords[0], inter3.coords[0], bridge_segments[0].coords[1],
                    bridge_segments[1].coords[0], inter4.coords[0], bridge_segments[1].coords[1]
                ]
                filling_poly = Polygon(filling_coords_ext_ext)
                poly_ext_ext = unary_union([self.poly_a, self.poly_b, filling_poly])

                # Control the area and topology for ext-ext
                if poly_ext_ext.is_valid:
                    if poly_ext_ext.area < best_polygon.area:
                        best_polygon = poly_ext_ext

                #********************  Int-Ext  ********************#
                # create filling polygon for int-ext
                filling_coords_int_ext = [
                    bridge_segments[0].coords[0], inter1.coords[0], bridge_segments[0].coords[1],
                    bridge_segments[1].coords[0], inter4.coords[0], bridge_segments[1].coords[1]
                ]
                filling_poly = Polygon(filling_coords_int_ext)
                poly_int_ext = unary_union([self.poly_a, self.poly_b, filling_poly])
                # Control the area and topology for int-ext
                if poly_int_ext.is_valid:
                    if poly_int_ext.area < best_polygon.area:
                        best_polygon = poly_int_ext

                #********************  Ext-Int  ********************#
                # create filling polygon for ext-int
                filling_coords_ext_int = [
                    bridge_segments[0].coords[0], inter3.coords[0], bridge_segments[0].coords[1],
                    bridge_segments[1].coords[0], inter2.coords[0], bridge_segments[1].coords[1]
                ]
                filling_poly = Polygon(filling_coords_ext_int)
                poly_ext_int = unary_union([self.poly_a, self.poly_b, filling_poly])
                # Control the area and topology for int-ext
                if poly_ext_int.is_valid:
                    if poly_ext_int.area < best_polygon.area:
                        best_polygon = poly_ext_int

                final_result = best_polygon
            return final_result
        except Exception as e:
            print("Error in the computation of intersections, return the convex hull.")
            return combined

def building_amalgamation(building1, building2):
    """
    Amalgamate two buildings by computing the most compact junction between them.
    
    The amalgamation algorithm was proposed by Nicolas Regnauld 
    :footcite:p:`regnauld_generalisation_1998` (pp. 150-154). The algorithm is based on the edges 
    of the convex hull that join both buildings. The algorithm chooses the links that create the smallest joined polygons.

    Parameters
    ----------
    building1 : the geometry of the first building to amalgamate.
    building2 : the geometry of the second building to amalgamate.
    
    Returns
    -------
    a Polygon

    See Also
    --------
    morphological_amalgamation :
        Amalgamates a list of close buildings.

    References
    ----------
    .. footbibliography::
    
    Examples
    --------
    >>> buildings = [Polygon([(1, 0), (9, 0), (9, 6), (1, 6), (1, 0)]), Polygon([(10, 0), (17, 0), (17, 6), (10, 6), (10, 0)])]
    >>> morphological_amalgamation(buildings, 1.0, 1.0)
    <POLYGON ((1.207 1.983, 2.547 5.885, 16.768 4.282, 15.42 0.148, 1.207 1.983))>
    """
    agg = BuildingAggregator(building1, building2)
    return agg.aggregate()