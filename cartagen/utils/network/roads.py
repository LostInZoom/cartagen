import shapely
from shapely.geometry.polygon import orient
from shapely.ops import linemerge, unary_union

from cartagen.utils.geometry.angle import angle_3_pts, angle_to_zero_pi

class Crossroad:
    """
    An object representing a crossroad from one or multiple face(s) of a network.
    """
    def __init__(self, roads, tree, *faces):
        self.face = unary_union(*faces)

        # Retrieve objects that intersects the considered network face using strtree
        self.original = tree.query(self.face, predicate='intersects').tolist()

        self.original_geoms = []
        for i in self.original:
            self.original_geoms.append(roads[i])

        # Fully dissolve and node the subnetwork
        unioned = unary_union(self.original_geoms)
        # Merge all contiguous lines
        merged = linemerge(unioned)

        self.network = []
        if merged.geom_type == 'LineString':
            self.network.append(merged)
        elif merged.geom_type == 'MultiLineString':
            for line in merged.geoms:
                self.network.append(line)

        rtree = shapely.STRtree(self.network)
        indexes = rtree.query(self.face)

        self.nodes, self.links = self.__create_graph_from_face(self.network, self.face)
        self.internals = self.__get_internal_roads(indexes, self.network, self.face)
        self.externals = self.__get_external_roads(indexes, self.network, self.face)

    def get_unchanged_roads(self, roads=None):
        """
        Return a list of index of unchanged roads, i.e. roads that have the same gemetry as the original.
        If the roads argument is None, returns all unchanged roads. Possible values are: 'externals' and 'internals'.
        Default is set to None.
        """
        unchanged = []
        # Loop through original geometry
        for oid, o in enumerate(self.original_geoms):
            # If externals roads are calculated
            if roads is None or roads == 'externals':
                # Loop through external roads
                for e in self.externals:
                    # Check if the external road geometry is the same as the original one
                    if shapely.equals(o, self.network[e]):
                        if self.original[oid] not in unchanged:
                            # If the id is not already added, append the id of the road to the result
                            unchanged.append(self.original[oid])
            # Same process for internal roads 
            if roads is None or roads == 'internals':
                for i in self.internals:
                    if shapely.equals(o, self.network[i]):
                        if self.original[oid] not in unchanged:
                            unchanged.append(self.original[oid])

        # Return the unchanged list of index
        return unchanged


    def __create_graph_from_face(self, network, face):
        """
        Return a list of nodes with their degree and a list a link between the index of the nodes.
        """
        coordinates = []
        nodes = []
        index = 0

        def update_nodes(node, index, coordinates):
            """
            Add a node to the list of nodes if it's not already present and if it's and internal node.
            i.e. it intersects with the network face considered.
            """
            if node in coordinates:
                # If the node is already in the coordinates list, increment its degree only
                i = coordinates.index(node)
                nodes[i][0] += 1
            else:
                # If it's not already there
                if shapely.intersects(face, shapely.Point(node)):
                    # Add to the node and coordinates list and increment the index
                    nodes.append([1, node])
                    coordinates.append(node)
                    index += 1

        # Loop though the lines in the provided network
        for link in network:
            # Retrieve all nodes of the line and keep the first and last node
            linenodes = link.coords
            firstnode = linenodes[0]
            lastnode = linenodes[-1]
            # Add the first and last to the list if not present
            update_nodes(firstnode, index, coordinates)
            update_nodes(lastnode, index, coordinates)

        links = []
        previous = None
        # Retrieve the oriented (so the start and end node of the exterior ring is always counterclockwise) linear ring of the network face
        linear_ring = orient(face).exterior.coords
        # Loop through the linear ring points
        for i, cc in enumerate(linear_ring):
            # If the point is present in the coordinates, it's a network node
            if cc in coordinates:
                # Retrieve the index of the node
                node = coordinates.index(cc)
                # Check if the 
                if len(nodes[node]) < 3:
                    # Get the previous and following point of the ring (if else is in case of the first and last point of the linear ring)
                    cp = linear_ring[-2] if i == 0 else linear_ring[i - 1]
                    cn = linear_ring[1] if i == len(linear_ring) - 1 else linear_ring[i + 1]
                    # Calculate the rounded interior angle of the polygon at the network node using the previous and following point
                    alpha = angle_3_pts(shapely.Point(cp), shapely.Point(cc), shapely.Point(cn))
                    # Add the angle to the node attributes
                    nodes[node].append(angle_to_zero_pi(alpha))

                # If it's not the first iteration, create the link with the previous node
                if previous is not None:
                    links.append([previous, node])

                # Set the previous node as the current one
                previous = node

        return nodes, links

    def __get_internal_roads(self, indexes, network, face):
        """
        Return the internal roads of a network. i.e. the roads that are contains by the face exterior ring.
        """
        internal = []
        for road in indexes:
            if shapely.contains(face.exterior, network[road]):
                internal.append(road)
        return internal

    def __get_external_roads(self, indexes, network, face):
        """
        Return the external roads of a network. i.e. the roads that only touches the face.
        """
        external = []
        for road in indexes:
            if shapely.touches(face.exterior, network[road]):
                external.append(road)
        return external