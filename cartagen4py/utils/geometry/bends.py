import shapely

from cartagen4py.utils.geometry.angle import *
from cartagen4py.utils.geometry.line import *

class Bend:
    """
    Create a bend object from a given linestring.
    """
    def __init__(self, line):
        self.side = get_bend_side(line)

        # Get the bend line as a vertex list
        coords = list(line.coords)
        # Calculate the base segment
        base = shapely.LineString([coords[0], coords[-1]])

        # Calculate the middle of the base segment
        middle = base.interpolate(base.length / 2)

        summit = None
        distmax = 0
        for i in range(1, len(coords) - 2):
            point = shapely.Point(coords[i])
            dist = point.distance(middle)
            if dist > distmax:
                summit = point
                distmax = dist

        self.summit = summit

class BendSerie:
    """
    Create a serie of bends from a given linestring object.
    """
    def __init__(self, line, sigma, sample):
        # Smooth the line to avoid unnecessary micro inflexion points
        smoothed = gaussian_smoothing(line, sigma, sample, densify=False)

        # Get inflexion points without first and last
        inflexion = get_inflexion_points(smoothed)[1:-1]

        # Get the list of input vertices
        vertices = list(line.coords)

        self.bends = []
        
        # Loop through vertices
        bend = []
        for i, vertex in enumerate(vertices):
            # Add the vertex to to the current bend
            bend.append(vertex)
            # If the vertex is an inflexion point
            if i in inflexion:
                # Add the bend to the list
                self.bends.append(Bend(shapely.LineString(bend)))
                # Restart a new bend
                bend = [vertex]

        # Add the last bend to the list
        self.bends.append(Bend(shapely.LineString(bend)))