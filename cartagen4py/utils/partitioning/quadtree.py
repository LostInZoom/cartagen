import shapely

class QuadTree:
    """
        Create a quad tree object
    """
    def __init__(self, objects, max_objects, max_depth):
        self.OBJECTS = objects
        self.MAX_OBJECTS = max_objects
        self.MAX_DEPTH = max_depth

        self.EXTENT = shapely.bounds(objects).tolist()

    def subdivide():
        pass