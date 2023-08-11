import shapely, itertools
from shapely import LineString

class SkeletonTIN:
    def __init__(self, polygon):
        """
        Create a polygon skeleton from a Delaunay triangulation.
        """
        # Stores future nodes and edges of the triangulation
        self.nodes = []
        self.edges = []

        # Stores future triangles composed of a list of three edges index
        # This list should be empty at the end of the skeletonization as
        # triangles are recursively removed from the list
        self.__triangles = []
        
        # Stores future bones and joints, i.e. links and nodes created by the skeletonization
        self.bones = []
        self.joints = []

        # Launching Delaunay's triangulation to populate nodes, edges and triangles
        self.__delaunay_triangulation(polygon)

        # Actual skeletonization.
        self.__skeletonize()

    def __skeletonize(self):
        """
        Recursively create the TIN skeleton.
        """

        def calculate_joint(edge, previous_joint):
            """
            Calculate the position of the joint on the edge.
            """
            # Retrieve the coordinates of the edge end and start point
            v = self.edges[edge]
            v1, v2 = self.nodes[v[0]], self.nodes[v[1]]
            x1, y1, x2, y2 = v1[0], v1[1], v2[0], v2[1]

            joint = [(x1 + x2) / 2, (y1 + y2) / 2]

            # Adding it to the list
            self.joints.append(shapely.Point(joint))

            # Create the current bone from the current and the next joint
            self.bones.append(shapely.LineString([previous_joint, joint]))

            return joint

        def recursive_bone_growth(edge, joint):
            """
            Recursively create bones and joints inside the skeleton.
            """
            # Retrieve the only left triangle having that edge, this is the current triangle
            current_triangle = [(i, t) for i, t in enumerate(self.__triangles) if edge in t][0]
            index, triangle = current_triangle[0], current_triangle[1]

            # Retrieve the triangles that share one edge with the current triangle
            followings = []
            virtuals = []
            for e in triangle:
                for t in self.__triangles:
                    if e in t and t != triangle:
                        followings.append(t)
                        virtuals.append(e)

            # Depending on the number of triangles sharing an edge with the current one,
            # we can determinate the kind of triangle it is
            typet = len(followings)

            # Here, without triangles left sharing this edge, it's an ear triangle
            if typet == 0:
                # Retrieve the two object edges of the triangle
                tedge = [self.edges[e] for e in self.__triangles[index] if e != edge]

                # Retrieve the common node of those two edges
                node = self.nodes[set.intersection(set(tedge[0]), *itertools.islice(tedge, 1, None)).pop()]
                
                # Creating end point
                end = shapely.Point(node)

                # Add the end node to the list
                self.joints.append(end)

                # Add the bone to the list
                self.bones.append(shapely.LineString([joint, end]))

                # Removes the triangle from the list
                self.__triangles.pop(index)

            # Here, with one triangles sharing this edge, it's a middle triangle
            elif typet == 1:
                virtual = virtuals[0]

                # Create the following joint in the skeleton, in the middle of the virtual edge
                nextjoint = calculate_joint(virtual, joint)

                # Remove the current triangle from the list
                self.__triangles.pop(index)

                # Continue the bone creation
                recursive_bone_growth(virtual, nextjoint)

            # Here, with two triangles sharing this edge, it's an interior triangle
            elif typet == 2:
                # Get all the nodes making the polygon
                tnodes = list({x for l in [self.edges[e] for e in triangle] for x in l})
                # Create the polygon and calculate the centroid
                centroid = shapely.Polygon([self.nodes[n] for n in tnodes]).centroid

                # Add the centroid to the joints list
                self.joints.append(centroid)
                # Create the line from the current joint to the centroid
                line = shapely.LineString([joint, centroid])
                # Add the line to the bones list
                self.bones.append(line)

                # Removes the current triangle from the list
                self.__triangles.pop(index)

                # Loop through virtual edges
                for virtual in virtuals:
                    # Creating the joint and bone
                    vjoint = calculate_joint(virtual, centroid)

                    # Continue the bone creation from this joint
                    recursive_bone_growth(virtual, vjoint)

            # Here, we have a problem
            else:
                raise Exception('There was a problem during the skeletonization.')

        start = None
        virtual = None

        # Loop through the triangles and break the loop when we find the first ear triangle,
        # i.e. the first triangle to have only one virtual edge
        for index, triangle in enumerate(self.__triangles):
            nb_virtual = 0
            for edge in triangle:
                for t in self.__triangles:
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
            # Retrieve the two object edges of the triangle
            tedge = [self.edges[e] for e in self.__triangles[start] if e != virtual]

            # Retrieve the common node of those two edges and add it to the list of joints
            node = self.nodes[set.intersection(set(tedge[0]), *itertools.islice(tedge, 1, None)).pop()]
            self.joints.append(shapely.Point(node))
            
            # Creating the joint, i.e. the point located in the middle of the virtual edge
            joint = calculate_joint(virtual, node)

            # Removes the triangle from the list
            self.__triangles.pop(start)

            # Starts the recursive bone creation
            recursive_bone_growth(virtual, joint)


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
                self.__triangles.append(triangle)