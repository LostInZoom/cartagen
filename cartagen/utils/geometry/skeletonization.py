import shapely, itertools
import numpy as np
from shapely.geometry import LineString
from copy import deepcopy
from cartagen.utils.geometry.line import gaussian_smoothing, get_segment_center, merge_linestrings

def skeletonize_natural(polygon, threshold_range=(0.7, 1.4)):
    """
    Create the natural TIN skeleton of a polygon.
    
    The natural skeleton means the entry points of the skeleton
    are naturally calculated during the skeleton creation, it depends
    on the Delaunay triangulation.

    This version is inspired by the enhanced TIN skeleton proposed by Wang :footcite:p:`wang:2010`

    Parameters
    ----------
    polygon : Polygon
        The polygon to compute the skeleton from.
    threshold_range : tuple, optional
        This value is only used for inner triangles calculated by the Delaunay triangulation.
        If two of the length ratio between each pair of edges (of an inner triangle) is outside the given range,
        the skeleton uses the middle of the line connecting the two center of the longest lines of the triangle,
        else it uses the centroid of the triangle

    Returns
    -------
    list of LineString

    See Also
    --------
    skeletonize_artificial : 
        Create a TIN skeleton by enforcing entry points.
    skeletonize_network :
        Create a TIN skeleton and blend it with the provided network.

    References
    ----------
    .. footbibliography::
    """
    # Generate the skeleton
    skeleton = SkeletonTIN(polygon, threshold_range)
    # Calculate the network from joinst and bones
    n = skeleton.create_network()

    return n

def skeletonize_artificial(polygon, entries, connection='joint', threshold_range=(0.7, 1.4)):
    """
    Create a TIN skeleton from a polygon by enforcing entry points.
    
    The artificial skeleton means the entry points of the skeleton
    are manually set and the natural entry points of the skeleton
    calculated during the skeleton creation are removed.

    This version is inspired by the enhanced TIN skeleton
    proposed by Wang :footcite:p:`wang:2010`

    Parameters
    ----------
    polygon : Polygon
        The polygon to compute the skeleton from.
    entries : list of Point
        The entry points of the skeleton.
    connection : str, optional
        Connection can be 'joint', which forces the entry point to connect to the closest joint,
        or can be 'interior', which first look if there is an interior triangle inside the triangulation
        and connect directly the entry point to its skeleton joint. If no interior triangle is found,
        apply 'joint' connection instead.
    threshold_range : tuple, optional
        This value is only used for inner triangles calculated by the Delaunay triangulation.
        If two of the length ratio between each pair of edges (of an inner triangle) is outside the given range,
        the skeleton uses the middle of the line connecting the two center of the longest lines of the triangle,
        else it uses the centroid of the triangle.

    Returns
    -------
    list of LineString

    See Also
    --------
    skeletonize_natural : 
        Create the natural TIN skeleton.
    skeletonize_network :
        Create a TIN skeleton and blend it with the provided network.

    References
    ----------
    .. footbibliography::
    """
    # Generate the skeleton
    skeleton = SkeletonTIN(polygon, threshold_range)
    # Storing the natural entries
    natural = skeleton.entries.copy()
    # Add the artificial entries
    skeleton.add_incoming_points(entries, connection)

    # Keep only natural entries different from provided entries
    remove = [ e for e in natural if e not in entries ]

    # Remove the natural entries
    skeleton.remove_entries(remove)

    # Calculate the network
    n = skeleton.create_network()

    return n

def skeletonize_network(polygon, network, sigma=None, attributes=None, threshold_range=(0.7, 1.4)):
    """
    Create an artificial TIN skeleton and blend it inside a network.
    
    This function first creates an artificial TIN skeleton by enforcing
    entry points derived from the provided network (the entry points are the extremities
    of the provided network touching the polygon ring). Then the network is blended
    with the skeleton and an optional gaussian smoothing is apply.

    This version is inspired by the enhanced TIN skeleton proposed by Wang :footcite:p:`wang:2010`

    Parameters
    ----------
    polygon : Polygon
        The polygon to compute the skeleton from.
    network : GeoDataFrame of LineString
        The network touching the ring of the polygon.
    sigma : float, optional
        Gaussian filter strength. By default, the skeleton is not smoothed.
    attributes : list of dict, optional
        A list of dict where the index matches the provided network. This dict represents the attributes
        of the lines to be propagated to the new network.
        If None is provided, the resulting network will lack the original attributes.
    threshold_range : tuple, optional
        This value is only used for inner triangles calculated by the Delaunay triangulation.
        If two of the length ratio between each pair of edges (of an inner triangle) is outside the given range,
        the skeleton uses the middle of the line connecting the two center of the longest lines of the triangle,
        else it uses the centroid of the triangle.

    Returns
    -------
    GeoDataFrame of LineString

    See Also
    --------
    skeletonize_natural : 
        Create the natural TIN skeleton.
    skeletonize_artificial :
        Create a TIN skeleton by enforcing entry points.

    References
    ----------
    .. footbibliography::
    """
    # Create the natural skeleton
    skeleton = SkeletonTIN(polygon, threshold_range)
    # Add incoming lines to the skeleton
    skeleton.add_incoming_lines(network)
    # Create the network
    skeleton.create_network()
    # Blend the network with the provided one
    skeleton.blend(attributes, sigma=sigma)

    return skeleton.blended

class SkeletonTIN:
    """
    Create a polygon skeleton from a Delaunay triangulation.
    
    The skeleton is represented by the lines (called bones) linking each middle point (called joints) of the edges
    of the triangles computed by the Delaunay triangulation.

    Parameters
    ----------
    polygon : shapely.Polygon
        The polygon to compute the skeleton from.
    threshold_range : tuple, Default=(0.7, 1.4)
        This value is only used for inner triangles calculated by the Delaunay triangulation.
        If two of the length ratio between each pair of edges (of an inner triangle) is outside the given range,
        the skeleton uses the middle of the line connecting the two center of the longest lines of the triangle,
        else it uses the centroid of the triangle (default value used by Wang, 2010).
    """
    def __init__(self, polygon, threshold_range=(0.7, 1.4)):
        # Stores polygon for future use
        self.polygon = polygon

        # Stores future nodes and edges of the triangulation
        self.nodes = []
        self.edges = []

        # Stores future triangles composed of a list of three edges index
        # This list should be empty at the end of the skeletonization as
        # triangles are recursively removed from the list
        self.__triangles = []
        
        # Stores future bones and joints,
        # i.e. links and nodes created by the skeletonization
        self.bones = []
        self.joints = []

        # Interiors joints are joints that are connected to more than two bones
        # i.e. joints created inside an interiors triangle during skeletonization
        # They are updated when adding and removing entry points
        self.interiors = []

        # Stores the skeleton as a network
        self.network = []

        # Storage for entry points
        self.entries = []

        self.incoming = None
        self.blended = []

        # Launching Delaunay's triangulation to populate nodes, edges and triangles
        self.__delaunay_triangulation(polygon)

        # Create the skeleton
        self.__skeletonize(polygon, threshold_range)

    def add_incoming_points(self, entries, connection='joint'):
        """
        Add new entry points that are not 'natural' skeleton entries,
        i.e. entries derived from the Delaunay triangulation.

        Parameters
        ----------
        entries : shapely Polygon
            The polygon to compute the skeleton from.
        connection : str optional
            Connection can be 'joint', which forces the entry point to connect to the closest joint,
            or can be 'interior', which first look if there is an interior triangle inside the triangulation
            and connect directly the entry point to its skeleton joint. If no interior triangle is found,
            apply 'joint' connection instead.
        """
        # Add the point
        def add_point(entry, connection):
            # If it doesn't exist, add it
            if entry not in self.joints:
                closest = None
                distance = None

                if connection == 'interior' and len(self.interiors) > 0:
                    joints = self.interiors.copy()
                else:
                    joints = self.joints.copy()

                # Looping through joints
                for joint in joints:
                    # Verify the joint is not an entry
                    if joint not in self.entries:
                        # Calculate distance between entry and joint
                        d = shapely.distance(entry, joint)
                        # If the previous distance has not been calculated yet or it's higher than the new one
                        if distance is None or d < distance:
                            # Updating it
                            distance = d
                            closest = joint

                # If a closest joint is found
                if closest is not None:
                    # Add the joint and entry to the lists
                    self.joints.append(entry)
                    self.entries.append(entry)

                    # Add the bone for connectivity
                    self.bones.append(shapely.LineString([entry, closest]))

                    # Add the closest point as a new interior joint
                    if closest not in self.interiors:
                        self.interiors.append(closest)


        if isinstance(entries, list):
            for entry in entries:
                add_point(entry, connection)
        else:
            if isinstance(entries, shapely.geometry.base.BaseGeometry):
                if entries.geom_type == 'Point':
                    add_point(entries, connection)
                else:
                    raise Exception('Cannot add entry point, wrong geometry type.')

    def add_incoming_lines(self, lines):
        """
        Add incoming lines to the skeleton by adding entry points to the 'natural' skeleton entries,
        i.e. entries derived from the Delaunay triangulation. Those lines must have their first or
        last node be the same as a node from the provided polygon to calculate the skeleton.

        Parameters
        ----------
        lines : Geopandas GeoDataFrame of LineStrings.
            The incoming lines to add to the skeleton.
        """
        entries = []
        boundary = list(self.polygon.boundary.coords)
        self.incoming = lines
        externals = [ i['geometry'] for i in lines ]

        # Loop through provided lines
        for line in externals:
            # Get start and end node of the line
            start, end = line.coords[0], line.coords[-1]
            # Add the start and/or end point of the line if it doesn't already exists
            if start in boundary:
                if start not in entries:
                    entries.append(shapely.Point(start))
            if end in boundary:
                if end not in entries:
                    entries.append(shapely.Point(end))

        # Add the list of new entry points
        self.add_incoming_points(entries)

        # Remove natural entries not in the provided ones
        self.remove_entries([ e for e in self.entries if e not in entries ])

    def remove_entries(self, entries):
        """
        Remove the provided points if they are entries of the skeleton.

        The skeleton is reconstructed by removing the provided skeleton entries as well as all the joints and bones
        that connects the entry to the rest of the skeleton. Be aware though, you can't remove an entry from
        a skeleton having only two entries.

        Parameters
        ----------
        entries : list of shapely Point.
            The entry points to remove from the skeleton.
        """
        # Recursively remove the bones and joints connecting the provided entry
        def recursive_removal(point):
            # Remove specific joint from the list
            def remove_joint(joint):
                if joint in self.joints:
                    self.joints.pop(self.joints.index(joint))
                    if joint in self.entries:
                        self.entries.pop(self.entries.index(joint))
            # Remove specific bone from the list
            def remove_bone(bone):
                bones = []
                if bone in self.bones:
                    self.bones.pop(self.bones.index(bone))

            # List all connected bones
            bones = []
            # Loop through all the bones
            for b in self.bones:
                # Add it to the list if it intersects the point
                if shapely.intersects(point, b):
                    bones.append(b)

            # If no bones were found it means the whole skeleton was removed
            if len(bones) == 0:
                # Removing the point from the joints
                remove_joint(point)
            # If one bone was found, continuing the recursive removal
            elif len(bones) == 1:
                bone = bones[0]
                start, end = shapely.Point(bone.coords[0]), shapely.Point(bone.coords[1])

                next_point = None
                # Get the next joint
                if shapely.equals(start, point):
                    next_point = end
                elif shapely.equals(end, point):
                    next_point = start

                if next_point is not None:
                    remove_joint(point)
                    remove_bone(bone)
                    recursive_removal(next_point)
            elif len(bones) == 2:
                # Here, the node intersects two bones, it's an interior joint
                # Recursion stops here, the interior joint is removed
                if point in self.interiors:
                    self.interiors.pop(self.interiors.index(point))

        # Start the recursion if the skeleton has more than two entries
        if len(self.entries) > 2:
            for e in entries:
                recursive_removal(e)

    def create_network(self):
        """
        Create the inside network from the skeleton bones and joints. Bones are merged together if only one joint links them.
        No smoothing or simplification are apply at this stage.
        Returns the network calculated.
        """
        bones = self.bones.copy()
        network = []

        def recursive_network_creation(point, polyline):
            # Add the current point to the polyline
            polyline.append(point)

            connected = []
            nextpoints = []
            # Loop through bones
            for bone in bones:
                # Retrieve start and end point of bone
                start, end = bone.coords[0], bone.coords[-1]
                # Add the start point or end point to the nextpoints list
                if start == point:
                    nextpoints.append(end)
                elif end == point:
                    nextpoints.append(start)
                else:
                    # If start or end doesn't match the given point, continue
                    continue

                # If start or end matches the given point, add the bone as connected
                connected.append(bone)

            # Stores the number of connexions
            connexions = len(connected)

            # If 1 connexions, continuing the polyline creation
            if connexions == 1:
                # Remove the connected bone from the bone list
                bones.pop(bones.index(connected[0]))
                # Continue the recursion with the same polyline
                recursive_network_creation(nextpoints[0], polyline)
            # If more than one connexion is found, it's a junction
            elif connexions > 1:
                # Add the current polyline to the network
                network.append(shapely.LineString(polyline))
                # Loop through all connected bones
                for index, b in enumerate(connected):
                    # Removes the considered connected bone
                    bones.pop(bones.index(b))
                    # Start the recursion again with a new polyline starting at the current point
                    recursive_network_creation(nextpoints[index], [point])
            # If no connexions were found
            else:
                # Add the current polyline to the network
                network.append(shapely.LineString(polyline))          

        # Starting the recursion at the first entry point
        if len(self.entries) > 0:
            coords = self.entries[0].coords[0]
            recursive_network_creation(coords, [])

        self.network = network
        return network

    def blend(self, attributes=None, sigma=None):
        """
        This method can only be used if incoming lines were provided using the add_incoming_lines() method.
        Blends the incoming lines with the created skeleton network. It relies on entry points.
        Parameters
        ----------
        attributes : list of dict optional.
            A list of dict where the index matches the procided incoming lines. This dict represents the attributes
            of the lines to be propagated to the new network. If None is provided, the resulting network will lack
            the original attributes. You can use any attributes propagation you see fit prior to calculate the skeleton.
        sigma : float optional.
            The strength of the gaussian smoothing applied to the resulting blended network. If None is provided, no smoothing
            is applied to the lines. Be careful that this method is destructive in the way that the vertices of incoming lines blended
            with the skeleton are moved, thus, the entry points are lost.
        """

        if self.incoming is None:
            raise Exception("Incoming lines were not provided during the skeletonization. Cannot blend with the network.")

        blended = deepcopy(self.incoming)
        network = self.network.copy()

        remove = []
        # Loop through entry points
        for entry in self.entries:
            # Counter for the degree of the entry point
            degree = 1
            # Will store the incoming line informations
            outside, outindex = None, None
            # Loop through incoming lines
            for iline, line in enumerate(self.incoming):
                # Here, the incoming line intersects the entry point
                if shapely.intersects(line['geometry'], entry):
                    # Increment the entry point degree and set the outside as the incoming line
                    degree += 1
                    outside = line['geometry']
                    outindex = iline

            # If the degree is 2, merge the incoming line with the skeleton line
            # The outside variable here should have been assigned only once
            if degree == 2:
                inside = None
                # Here, retrieve the skeleton line connected to the entry point
                for skline in network:
                    if shapely.intersects(skline, entry):
                        inside = skline

                # Update the geometry of the incoming line by merging it with the intersecting skeleton line
                blended[outindex]['geometry'] = merge_linestrings(outside, inside)

        # Add the unmodified skeleton lines to the results
        if attributes is not None:
            blended += [{ **attributes, **{"geometry": n} } for n in network if n not in remove]
        else:
            blended += [{ "geometry": n } for n in network if n not in remove]

        remove = []
        # Here, handle lines contained by an other
        for i1, b1 in enumerate(blended):
            for i2, b2 in enumerate(blended):
                geom1, geom2 = b1['geometry'], b2['geometry']
                if geom1 != geom2:
                    if shapely.contains(geom2, geom1):
                        if i1 not in remove:
                            remove.append(i1)

        # Here, handle specific case of two entry of degree 2
        if len(self.entries) == 2 and (len(blended) - len(remove)) == 2:
            final = blended[0]
            # Merge all the lines into one
            geometries = []
            for b in blended:
                geometries.append(b['geometry'])
            final['geometry'] = shapely.ops.linemerge(geometries)
            blended = [final]
        else:
            # Else, retrieve the blended network without contained lines
            blended = [b for i, b in enumerate(blended) if i not in remove]

        # Apply a gaussian smoothing to the blended network if selected
        if sigma is not None:
            for i, n in enumerate(blended):
                geom = n['geometry']
                blended[i]['geometry'] = gaussian_smoothing(geom, sigma=sigma)

        self.blended = blended
        return blended

    def __skeletonize(self, polygon, threshold_range):
        """
        Recursively create the TIN skeleton.
        Populate joints property, i.e. all the vertex of the polylines of the skeleton
        Populate bones property, i.e. all the segment composing the polylines of the skeleton
        """
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
                self.entries.append(end)

                # Add the bone to the list
                self.bones.append(shapely.LineString([joint, end]))

                # Removes the triangle from the list
                self.__triangles.pop(index)

            # Here, with one triangles sharing this edge, it's a middle triangle
            elif typet == 1:
                virtual = virtuals[0]

                # Create the following joint in the skeleton, in the middle of the virtual edge
                nextjoint = self.__calculate_joint(virtual, joint)

                # Remove the current triangle from the list
                self.__triangles.pop(index)

                # Continue the bone creation
                recursive_bone_growth(virtual, nextjoint)

            # Here, with two triangles sharing this edge, it's an interior triangle
            elif typet == 2:
                # Get all the nodes making the polygon
                tnodes = list({x for l in [self.edges[e] for e in triangle] for x in l})

                # Calculate edge ratios to reduce zigzags or not
                center = self.__retrieve_center(tnodes, threshold_range)

                # Add the center to the joints list and the interior list
                self.joints.append(center)
                self.interiors.append(center)

                # Create the line from the current joint to the center
                line = shapely.LineString([joint, center])
                # Add the line to the bones list
                self.bones.append(line)

                # Removes the current triangle from the list
                self.__triangles.pop(index)

                # Loop through virtual edges
                for virtual in virtuals:
                    # Creating the joint and bone
                    vjoint = self.__calculate_joint(virtual, center)

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
            self.entries.append(shapely.Point(node))
            
            # Creating the joint, i.e. the point located in the middle of the virtual edge
            joint = self.__calculate_joint(virtual, node)

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

    def __retrieve_center(self, nodes, threshold_range):
        """
        Return the center of the inner triangle. It depends on the given range.
        If two of the length ratio between each pair of edges is outside the given range,
        returns the middle of the line connecting the two center of the longest lines of the triangle,
        else it returns the centroid of the triangle.
        """
        # Get the two longest lines of an array of linestrings
        def get_two_longest_lines(lines):
            # Create a list of line lengths
            lengths = []
            for line in lines:
                lengths.append(line.length)

            # Convert the list to a numpy array and get the sorted indexes
            l = np.array(lengths)
            indexes = l.argsort()

            # Retrieve the two longest lines
            return lines[indexes[-1]], lines[indexes[-2]]

        count = 0
        lines = []
        # Loop through triangle nodes
        for i, n in enumerate(nodes):
            # Retrieve the node coordinates and the two following nodes coordinates
            n, n1, n2 = self.nodes[n], None, None
            # Following nodes depends on the index of the current node
            if i == 0:
                n1, n2 = self.nodes[nodes[1]], self.nodes[nodes[2]]
            elif i == 1:
                n1, n2 = self.nodes[nodes[2]], self.nodes[nodes[0]]
            else:
                n1, n2 = self.nodes[nodes[0]], self.nodes[nodes[1]]
            
            # Create the lines from current node to next, and next to last
            l1 = shapely.LineString([n, n1])
            l2 = shapely.LineString([n1, n2])

            # Add the two lines if they don't already exists
            if l1 not in lines:
                lines.append(l1)
            if l2 not in lines:
                lines.append(l2)

            # Calculate the length ratio between those two lines
            ratio = l1.length / l2.length

            # If the ratio is outside the given range, increment the counter by 1
            if ratio < threshold_range[0] or ratio > threshold_range[1]:
                count += 1
        
        # If more than 2 ratio are outside the given range...
        if count < 2:
                # ...return the centroid of the triangle
            return shapely.Polygon([self.nodes[n] for n in nodes]).centroid
        else:
            # Else, retrieve the two longest lines
            line1, line2 = get_two_longest_lines(lines)
            # Get both of those lines' center
            p1 = get_segment_center(line1)
            p2 = get_segment_center(line2)
            # Return the center of the line formed by the two centers
            return get_segment_center(shapely.LineString([p1, p2]))

    def __calculate_joint(self, edge, previous_joint):
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