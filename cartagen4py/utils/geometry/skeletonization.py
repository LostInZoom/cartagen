import shapely
from shapely import LineString

class SkeletonTIN:
    def __init__(self, polygon):
        """
        Create a polygon skeleton from a Delaunay triangulation.
        """
        self.triangles = []
        self.edges = []
        self.nodes = []

        self.bones = []
        self.joints = []

        self.edges_object = []
        self.edges_virtual = []
        self.triangles_ear = []
        self.triangles_middle = []
        self.triangles_interior = []

        self.__delaunay_triangulation(polygon)
        self.__skeletonize()
        
        # self.__compute_edges(polygon)

        # Check if multiple triangles have been calculated
        # If only one triangle is present, it obviously means
        # the provided polygon is only a simple triangle
        # if len(self.triangles) > 1:
        #     # Derive the triangles from the edges
        #     self.__compute_triangles(polygon)
        #     # Compute the TIN skeleton
        #     self.skeleton = self.__skeletonize()

    def __skeletonize(self):
        """
        Recursively create the TIN skeleton.
        """
        start = None
        virtual = None

        # Loop through the triangles and break the loop when we find the first ear triangle,
        # i.e. the first triangle to have only one virtual edge
        for index, triangle in enumerate(self.triangles):
            nb_virtual = 0
            for edge in triangle:
                for t in self.triangles:
                    for e in t:
                        if t != triangle and e == edge:
                            # Retrieve the virtual edge index
                            virtual = edge
                            nb_virtual += 1
            if nb_virtual == 1:
                start = index
                break

        # If an ear triangle was found, starting the skeletonization
        if start is not None and virtual is not None:
            nextt = None
            for triangle in self.triangles:
                if virtual in triangle:
                    nextt = triangle
            print(nextt)

    def __delaunay_triangulation(self, polygon):
        """
        Create the delaunay triangulation and populate triangles, edges and nodes.
        The triangles created are a list of edge index.
        The edges created are a list of two nodes index.
        The nodes created are a list of two coordinates x, y.
        """
        index_node = 0
        index_edge = 0

        # Retrieve the index of the concerned node or add it to the list
        def handle_node(coords, index_node):
            # If the node already exists
            if coords in self.nodes:
                # Get the index of the node
                node = self.nodes.index(coords)
            # If it doesn't exists
            else:
                # Append the coordinates to the nodes list
                self.nodes.append(coords)
                # Set the node as the current index and increment it by one
                node = index_node
                index_node += 1
            # Return the node index and the current index
            return node, index_node

        # Retrieve the index of the concerned edge or add it to the list
        def handle_edge(start, end, index_edge):
            # If an edge already exists, get the index of the edge
            if [start, end] in self.edges:
                edge = self.edges.index([start, end])
            elif [end, start] in self.edges:
                edge = self.edges.index([end, start])
            else:
                # Add the edge to the list
                self.edges.append([start, end])
                # Set the edge as the current index and increment it by one
                edge = index_edge
                index_edge += 1
            # Return the edge index and the current index
            return edge, index_edge

        # Retrieve the boundary of the polygon
        boundary = polygon.boundary

        # Compute de delaunay triangulation
        delaunay = shapely.delaunay_triangles(polygon).geoms

        # Looping through all triangles
        for face in delaunay:
            # Check if the triangle is contained within the polygon
            if shapely.contains(polygon, face):
                # Retrieve the boundary of the triangle
                boundary = face.boundary.coords
                # Retrieve edges of the triangles
                edges = [LineString(boundary[k:k+2]) for k in range(len(boundary) - 1)]
                triangle = []
                # Loop through the triangle edges
                for e in edges:
                    # Get the start and end nodes index or create them
                    start, index_node = handle_node(e.coords[0], index_node)
                    end, index_node = handle_node(e.coords[1], index_node)
                    # Get the edge index or create it
                    edge, index_edge = handle_edge(start, end, index_edge)
                    # Append the edge to the triangle
                    triangle.append(edge)
                # Append the triangle to the list 
                self.triangles.append(triangle)