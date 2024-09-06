import numpy as np

from shapely import Point, LineString

from cartagen.utils.geometry.angle import get_curvature
from cartagen.utils.geometry.line import gaussian_smoothing, get_bend_side, inflexion_points

class Bend:
    """
    Create a bend object from a given linestring.
    """
    def __init__(self, line):
        self.bend = line
        self.side = get_bend_side(line)
        self.base = None
        self.middle = None
        self.height = None
        self.width = None
        self.summit = None
        self.isummit = None
        self.orientation = None

    def get_base(self):
        """
        Get the base of the bend. i.e. the line connecting the start and end point
        of the bend.
        """
        coords = list(self.bend.coords)
        self.base = LineString([coords[0], coords[-1]])
        return self.base

    def get_middle(self):
        """
        Get the middle of the base of the bend.
        """
        base = self.base
        if base is None:
            base = self.get_base()
        self.middle = base.interpolate(base.length / 2)
        return self.middle

    def get_width(self):
        """
        Get the width of the bend. i.e. the length of the line connecting the start
        and end point of the bend.
        """
        base = self.base
        if base is None:
            base = self.get_base()
        self.width = base.length
        return self.width

    def get_height(self):
        """
        Get the height of the bend.
        """
        summit = self.summit
        if summit is None:
            summit = self.get_summit()
        middle = self.middle
        if middle is None:
            middle = self.get_middle()
        self.height = LineString([middle, summit]).length
        return self.height

    def get_summit(self):
        """
        Calculate the summit of the bend.
        """
        if self.summit is not None:
            return self.summit
        
        # Smooth the bend
        smoothed = gaussian_smoothing(self.bend, densify=False)
        # Get the vertices list
        scoords = list(smoothed.coords)

        curvature_max = 0
        smooth_summit = None
        summit = None
        isummit = None
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
        for i, c in enumerate(list(self.bend.coords)):
            d = Point(c).distance(smooth_summit)
            if d < distance:
                # Update distance and nearest index if below existing
                distance = d
                summit = Point(c)
                isummit = i

        self.isummit = isummit
        self.summit = summit
        return summit

    def get_orientation(self):
        """
        Calculate the orientation of the bend.
        """
        summit = self.summit
        if summit is None:
            summit = self.get_summit()
        middle = self.middle
        if middle is None:
            middle = self.get_middle()

        # Calculate the orientation of the bend
        orientation = np.arctan2(summit.y - middle.y, summit.x - middle.x)
        # Set the orientation between 0 and 2pi
        modulo = orientation % (2 * np.pi)
        if modulo >= 0:
            self.orientation = abs(modulo)
        else:
            self.orientation = modulo + 2 * np.pi
        return self.orientation

    def get_symetry(self):
        """
        Retrieve the measure of the bend symetry. Value between 0 and 1.
        """
        # Get the summit
        summit = self.summit
        if summit is None:
            summit = self.get_summit()

        # Get vertices coordinates
        coords = list(self.bend.coords)
        line = []
        l1, l2 = None, None

        # Loop through coordinates
        for i, c in enumerate(coords):
            # Append the coordinate to the list
            line.append(c)
            # If it's the summit
            if i == self.isummit:
                # Set l1 to be the length of the first half of the bend
                l1 = LineString(line).length
                # Restart the line
                line = [c]
        # Set l2 to be the length of the second part of the bend
        l2 = LineString(line).length
        # Calculate the symetry measure
        self.symetry = min(l1, l2) / max(l1, l2)
        return self.symetry
        
class BendSerie:
    """
    Create a serie of bends from a given linestring object.
    """
    def __init__(self, line, sigma=30, sample=None):
        # Smooth the line to avoid unnecessary micro inflexion points
        smoothed = gaussian_smoothing(line, sigma, sample, densify=False)

        # Get inflexion points without first and last
        inflexion = inflexion_points(smoothed)[1:-1]

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