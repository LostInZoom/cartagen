# this file contains a data structure useful for several point set operations: the PointSetQuadTree, proposed by Bereuter & Weibel (2012)
import numpy as np
from shapely.geometry import Polygon, Point

class PointSetQuadTree():

    """A class implementing a point set quadtree."""

    def __init__(self, envelope, max_points=1, depth=0):
        """Initialize this node of the quadtree.

        envelope is a shapely Polygon object defining the region from which points are
        placed into this node; max_points is the maximum number of points the
        node can hold before it must divide (branch into four more nodes);
        depth keeps track of how deep into the quadtree this node lies.

        """

        self.envelope = envelope
        self.max_points = max_points
        self.points = []
        self.depth = depth
        # A flag to indicate whether this node has divided (branched) or not.
        self.divided = False
    
    def __str__(self):
        """Return a string representation of this node of the QuadTree, suitably formatted."""
        sp = ' ' * self.depth * 2
        s = str(self.boundary) + '\n'
        s += sp + ', '.join(str(point) for point in self.points)
        if not self.divided:
            return s
        return s + '\n' + '\n'.join([
                sp + 'nw: ' + str(self.nw), sp + 'ne: ' + str(self.ne),
                sp + 'se: ' + str(self.se), sp + 'sw: ' + str(self.sw)])
    
    def divide(self):
        """Divide (branch) this node by spawning four children nodes."""

        cx, cy = self.envelope.centroid.x, self.envelope.centroid.y
        length = abs(self.envelope.bounds[0]) - (self.envelope.bounds[2])
        # The boundaries of the four children nodes are "northwest",
        # "northeast", "southeast" and "southwest" quadrants within the
        # boundary of the current node.
        self.sw = PointSetQuadTree(Polygon([(cx - length/2, cy - length/2), (cx, cy - length/2), (cx, cy), (cx - length/2, cy), 
                                            (cx - length/2, cy - length/2)]),
                                    self.max_points, self.depth + 1)
        self.se = PointSetQuadTree(Polygon([(cx, cy - length/2), (cx + length/2, cy - length/2), (cx + length/2, cy), (cx, cy), 
                                            (cx, cy - length/2)]),
                                    self.max_points, self.depth + 1)
        self.nw = PointSetQuadTree(Polygon([(cx - length/2, cy), (cx, cy), (cx, cy + length/2), (cx - length/2, cy + length/2), 
                                            (cx - length/2, cy)]),
                                    self.max_points, self.depth + 1)
        self.ne = PointSetQuadTree(Polygon([(cx, cy), (cx + length/2, cy), (cx + length/2, cy + length/2), (cx, cy + length/2), 
                                            (cx, cy)]),
                                    self.max_points, self.depth + 1)
        self.divided = True
    
    def insert(self, point):
        """Try to insert a point from a GeoDataFrame into this QuadTree."""

        if not self.envelope.contains(point['geometry']):
            # The point does not lie inside boundary: bail.
            return False
        if len(self.points) < self.max_points:
            # There's room for our point without dividing the QuadTree.
            self.points.append(point)
            return True

        # No room: divide if necessary, then try the sub-quads.
        if not self.divided:
            self.divide()
            for prev_point in self.points:
              (self.ne.insert(prev_point) or
                self.nw.insert(prev_point) or
                self.se.insert(prev_point) or
                self.sw.insert(prev_point))  

        return (self.ne.insert(point) or
                self.nw.insert(point) or
                self.se.insert(point) or
                self.sw.insert(point))
    
    def __len__(self):
        """Return the number of points in the quadtree."""

        npoints = len(self.points)
        if self.divided:
            npoints += len(self.nw)+len(self.ne)+len(self.se)+len(self.sw)
        return npoints

    def draw(self, ax):
        """Draw a representation of the quadtree on Matplotlib Axes ax."""

        x1, y1, x2, y2 = self.envelope.bounds
        ax.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1], c='k', lw=1)
        if self.divided:
            self.nw.draw(ax)
            self.ne.draw(ax)
            self.se.draw(ax)
            self.sw.draw(ax)
    
    def populate(self, geodataframe):
        """Populate the quadtree with the points contained in a GeoDataFrame, with a Point geometry"""

        for index, feature in geodataframe.iterrows():
            self.insert(feature)

    def get_all_points(self):
        all_points = []
        for point in self.points:
            all_points.append([point, self.depth])
        if self.divided:
            for point, depth in self.nw.get_all_points():
                all_points.append([point, depth])
            for point, depth in self.sw.get_all_points():
                all_points.append([point, depth])
            for point, depth in self.ne.get_all_points():
                all_points.append([point, depth])
            for point, depth in self.se.get_all_points():
                all_points.append([point, depth])
        return all_points
