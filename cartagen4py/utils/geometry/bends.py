import numpy as np

from shapely import Point, LineString
from shapely.ops import nearest_points

from cartagen4py.utils.geometry.angle import *
from cartagen4py.utils.geometry.line import *

class Bend:
    """
    Create a bend object from a given linestring.
    """
    def __init__(self, line):
        self.bend = line
        self.side = get_bend_side(line)

        # Get the bend line as a vertex list
        coords = list(self.bend.coords)

        # Calculate the base segment
        base = LineString([coords[0], coords[-1]])

        self.width = base.length

        # Calculate the base point of the base segment
        self.middle = base.interpolate(base.length / 2)

        # Smooth the bend
        smoothed = gaussian_smoothing(self.bend, densify=False)
        # Get the vertices list
        scoords = list(smoothed.coords)

        curvature_max = 0
        smooth_summit = None
        # Loop through vertices without first and last
        for i in range(1, len(scoords) - 2):
            # Get the previous, current and next vertex
            p1, p2, p3 = Point(scoords[i - 1]), Point(scoords[i]), Point(scoords[i + 1])
            # Calculate the curvature
            curvature = get_curvature(p1, p2, p3)
            # Update the max curvature if above
            if curvature > curvature_max:
                curvature_max = curvature
                smooth_summit = p2

        # Set the distance to infinite
        distance = float("inf")
        for c in coords:
            d = Point(c).distance(smooth_summit)
            if d < distance:
                # Update distance and nearest index if below existing
                distance = d
                self.summit = Point(c)

        # Calculate the orientation of the bend
        orientation = np.arctan2(self.summit.y - self.middle.y, self.summit.x - self.middle.x)

        # Set the orientation between 0 and 2pi
        modulo = orientation % (2 * np.pi)
        if modulo >= 0:
            self.orientation = abs(modulo)
        else:
            self.orientation = modulo + 2 * np.pi

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
                self.bends.append(Bend(LineString(bend)))
                # Restart a new bend
                bend = [vertex]

        # Add the last bend to the list
        self.bends.append(Bend(LineString(bend)))