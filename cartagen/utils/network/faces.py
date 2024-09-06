import numpy as np
import shapely

class NetworkFace:
    """
    An object representing a network face along with specific geometric properties.
    """
    def __init__(self, polygon):
        ring = polygon.exterior

        # Area and perimeter
        self.area = polygon.area
        self.perimeter = ring.length

        # Compactness
        self.compactness = 4 * np.pi * self.area / (self.perimeter * self.perimeter)

        # Concavity, which is the factor between polygon surface and its convex hull surface
        hull = shapely.convex_hull(polygon)
        self.concavity = self.area / hull.area

        # Calculate length, width and elongation
        points = ring.coords
        ssr = shapely.MultiPoint(points).minimum_rotated_rectangle.exterior.coords

        x0, x1, x2 = ssr[0][0], ssr[1][0], ssr[2][0]
        y0, y1, y2 = ssr[0][1], ssr[1][1], ssr[2][1]

        length = np.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0))
        width = np.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))

        if length < width:
            length, width = width, length

        self.length = length
        self.width = width
        self.elongation = length / width

        # Calculate huber's width
        insqrt = (self.perimeter * self.perimeter) - (16.0 * self.area)
        if insqrt < 0:
            self.huber = 100
        else:
            self.huber = (self.perimeter - np.sqrt(insqrt)) / 4
