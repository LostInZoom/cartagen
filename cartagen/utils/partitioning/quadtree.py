from shapely.geometry import Polygon, LineString

def partition_quadtree(objects, max_points=1, max_depth=None):
    """
    Partition objects using a quadtree.

    This algorithm inserts the centroids (or surface points) of the provided
    objects into a quadtree built over the extent of the dataset, then collects
    the leaf cells and assigns each object to the leaf it was placed in.

    Parameters
    ----------
    objects : GeoDataFrame
        The objects to partition. Supports Point, LineString, and Polygon
        geometries.
    max_points : int, optional
        Maximum number of points a quadtree node can hold before it
        subdivides. Lower values produce finer partitions.
        Defaults to 1.
    max_depth : int, optional
        Maximum depth the quadtree is allowed to reach. Nodes at this depth
        will not subdivide further, even if they exceed ``max_points``.
        If set to None, depth is unlimited.
        Defaults to None.

    Returns
    -------
    partition : tuple
        A tuple containing two elements:

        #. A list of lists of index ordered by the quadtree leaf cells
        #. A list of the geometry of the quadtree leaf cells

    See Also
    --------
    partition_grid :
        Partition objects using a grid of a given shape.
    partition_networks :
        Partition objects using one or multiple networks.

    Examples
    --------
    >>> points = gpd.GeoDataFrame(geometry=[Point(1, 1), Point(3, 3)])
    >>> partition_quadtree(points, max_points=1)
    ([[0], [1]], [<POLYGON (...)>, <POLYGON (...)>])
    """
    xmin, ymin, xmax, ymax = objects.total_bounds

    # The quadtree envelope must be square for correct subdivision
    side = max(xmax - xmin, ymax - ymin)
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    half = side / 2
    envelope = Polygon([
        (cx - half, cy - half),
        (cx + half, cy - half),
        (cx + half, cy + half),
        (cx - half, cy + half),
    ])

    tree = PointSetQuadTree(envelope, max_points=max_points)

    # Build centroid proxies — store original index alongside geometry
    centroids = []
    for i, row in objects.iterrows():
        geom = row["geometry"]
        if geom.geom_type == "Polygon":
            centroids.append((i, geom.point_on_surface()))
        else:
            centroids.append((i, geom.centroid))

    # Insert synthetic GeoDataFrame-like dicts so the tree accepts them
    for orig_idx, point in centroids:
        tree.insert({"geometry": point, "_orig_idx": orig_idx})

    # Collect every leaf node envelope, respecting max_depth
    def _collect_leaves(node):
        """Recursively yield leaf PointSetQuadTree nodes."""
        if not node.divided or (max_depth is not None and node.depth >= max_depth):
            yield node
        else:
            for child in (node.nw, node.ne, node.se, node.sw):
                yield from _collect_leaves(child)

    leaves = list(_collect_leaves(tree))

    # Build the partition: one entry per leaf, preserving leaf order
    # When max_depth is set, a capped node may contain multiple points,
    # so all of them are gathered recursively from its subtree.
    def _collect_points(node):
        """Recursively collect all points in a node's subtree."""
        pts = list(node.points)
        if node.divided:
            for child in (node.nw, node.ne, node.se, node.sw):
                pts.extend(_collect_points(child))
        return pts

    index_lists = []
    geometries = []
    for leaf in leaves:
        points = _collect_points(leaf)
        indices = [pt["_orig_idx"] for pt in points]
        if indices:  # skip empty leaves
            index_lists.append(indices)
            geometries.append(leaf.envelope)

    return (index_lists, geometries)

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
        length = abs((self.envelope.bounds[0]) - (self.envelope.bounds[2]))
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
        if not self.envelope.intersects(point['geometry']):
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

        return (
            self.ne.insert(point) or
            self.nw.insert(point) or
            self.se.insert(point) or
            self.sw.insert(point)
        )
    
    def __len__(self):
        """Return the number of points in the quadtree."""

        npoints = len(self.points)
        if self.divided:
            npoints += len(self.nw)+len(self.ne)+len(self.se)+len(self.sw)
        return npoints

    def draw(self, ax, max_depth=None):
        """Draw a representation of the quadtree on Matplotlib Axes ax."""
        linewidth, color = 0.5, 'gray'

        if max_depth is not None:
            if self.depth <= max_depth:
                linewidth, color = 1, 'blue'
        
        x1, y1, x2, y2 = self.envelope.bounds
        ax.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1], color=color, linewidth=linewidth)
        if self.divided:
            self.nw.draw(ax, max_depth)
            self.ne.draw(ax, max_depth)
            self.se.draw(ax, max_depth)
            self.sw.draw(ax, max_depth)

    def geometry(self, max_depth=None, index=None):
        """Get the geometry of the quadtree as individual lines for drawing"""
        lines = []
        border = list(self.envelope.boundary.coords)
        borderline = []
        for i in range(0, len(border) - 1):
            borderline.append(LineString([border[i], border[i + 1]]))

        if index is None:
            lines += borderline
        else:
            lines.append(borderline[index])

        if max_depth is not None and self.depth >= max_depth:
            return lines
        else:
            if self.divided:
                lines += self.nw.geometry(max_depth, 0)
                lines += self.ne.geometry(max_depth, 3)
                lines += self.se.geometry(max_depth, 2)
                lines += self.sw.geometry(max_depth, 1)
            return lines
    
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
