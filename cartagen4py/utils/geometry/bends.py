import shapely

from cartagen4py.utils.geometry.angle import *
from cartagen4py.utils.geometry.line import *

class BendSerie:
    """
    Create a serie of bends from a given linestring object.
    """
    def __init__(self, line):
        self.inflexion = inflexion_points(line)
        self.bends = self.__detect_bends(line)

    def get_bend_sides(self):
        return [get_bend_side(x) for x in self.bends]

    def __detect_bends(self, line):
        """
        Detect bends from a given LineString.
        """
        # Get the list of input vertices
        vertices = list(line.coords)

        bends = []
        bend = []

        # Loop through vertices
        for i, vertex in enumerate(vertices):
            # Add the vertex to to the current bend
            bend.append(vertex)
            # If the vertex is an inflexion point
            if i in self.inflexion:
                # Add the bend to the list
                bends.append(shapely.LineString(bend))
                # Restart a new bend
                bend = [vertex]

        # Add the last bend to the list
        bends.append(shapely.LineString(bend))
        
        return bends