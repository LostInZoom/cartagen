import shapely
import geopandas as gpd
import networkx as nx
from shapely.ops import unary_union, linemerge

from cartagen.utils.geometry.distances import group_intersecting
from cartagen.utils.geometry.line import (
    get_line_middle_point, merge_linestrings,
    resample_line, gaussian_smoothing
)
from cartagen.utils.network.graph import create_graph

def spinalize_polygon(polygon, densify=None, sigma=None, entries=None, structural=None):
    """
    Collapse a single polygon into one or multiple lines.

    This algorithm creates the spine of a polygon to generate a line.
    It uses a Voronoi diagram from the list
    of vertices of the polygon. Then, depending on the provided
    entries, it calculates the closest path (using Djikstra) between
    those entries or between all vertexes.

    More information about the usage of this algorithm in
    Touya & Girres :footcite:p:`touya:2013` :footcite:p:`touya:2014`

    Parameters
    ----------
    polygon : Polygon
        The polygon to spinalize.
    densify : float, optional
        The densification step used to resample the polygon ring.
        If not provided, it can sometimes render the spinalization
        impossible because the generated Voronoi diagram doesn't have
        enough nodes.
        If set to a too high value and no entries are set, computation
        time can skyrocket.
    sigma : float, optional
        Gaussian filter strength.
        If set to None, no gaussian smoothing is applied.
    entries : list of Point, optional
        A list of entry points to the polygon.

        - **If left to None**, the spine is the longest path between all the nodes
          of the Voronoi diagram (this can be really expensive in terms of resources).
        - **If one is set**, the spine is the longest path between this point and all the
          nodes of the Voronoi diagram.
        - **If two are set**, the spine is the path between those two nodes on the Voronoi
          diagram.
        - **If more than two are set**, multiple spines are calculated between each nodes on
          the Voronoi diagram.

        The last two are the most efficient in terms of computation time because it doesn't
        calculate shortest paths to each node.

    structural : list of Point, optional
        A list of point inside the polygon that will be kept by the algorithm.
        Each point is snapped to the closest vertex on the spine before the
        gaussian smoothing is applied in order to keep them.
        This is useful to preserve right angles inside the polygon for example.

    Returns
    -------
    list of LineString

    Warning
    -------
    The ``densify`` and ``entries`` parameters, if not set correctly,
    can lead to extreme computation time. Read carefully the documentation
    to properly use this algorithm.

    See Also
    --------
    spinalize_polygons:
        Collapse multiple polgygons into one or multiple lines.
    resample_line :
        Densify a line by adding vertices.
    gaussian_smoothing:
        Smooth a line and attenuate its inflexion points.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> polygon = Polygon([(0, 0), (0, 10), (20, 10), (20, 0), (0, 0)])
    >>> spinalize_polygon(polygon)
    [<LINESTRING (0 5, 10 5, 20 5)>]
    """
    def __get_closest_node(nodes, point):
        distance = None
        closest = None
        for node in nodes:
            ndist = shapely.distance(shapely.Point(nodes[node]['coords']), point)
            if distance is not None:
                if ndist < distance:
                    distance = ndist
                    closest = node
            else:
                distance = ndist
                closest = node
        return closest

    real_structural = []
    if structural is not None:
        for s in structural:
            if shapely.intersects(polygon, s):
                real_structural.append(s)
    
    # Remove holes from the polygon
    polygon = shapely.Polygon(polygon.exterior.coords)

    if densify is not None:
        polygon = shapely.Polygon(resample_line(polygon.boundary, densify, True))

    # Calculate Voronoi diagram
    voronoi = list(shapely.voronoi_polygons(polygon, only_edges=True, extend_to=polygon).geoms)

    lines = []
    for v in voronoi:
        l = shapely.intersection(polygon, v)
        if len(l.coords) > 0:
            lines.append({'geometry': l})

    graph = create_graph(gpd.GeoDataFrame(lines))
    nodes = graph.nodes

    spines = []

    if entries is None or len(entries) == 0:
        longest = None
        length = 0
        # Loop through each node (source)
        for source in nodes:
            # Loop again through each node (target)
            for target in nodes:
                if source != target:
                    # Check if the pair of nodes has a path
                    if nx.has_path(graph, source, target):
                        # Calculate the path
                        dijkstra = nx.shortest_path(graph, source=source, target=target, weight='weight')
                        # Calculate the length of the path
                        plength = nx.path_weight(graph, dijkstra, 'weight')
                        # Check if path length is longer than the longest path length
                        if plength > length:
                            longest = dijkstra
                            length = plength
        
        spines.append(shapely.LineString([ nodes[node]['coords'] for node in longest ]))
    else:
        if len(entries) == 1:
            longest = None
            length = 0
            source = __get_closest_node(nodes, entries[0])

            for target in nodes:
                if source != target:
                    # Check if the pair of nodes has a path
                    if nx.has_path(graph, source, target):
                        # Calculate the path
                        dijkstra = nx.shortest_path(graph, source=source, target=target, weight='weight')
                        # Calculate the length of the path
                        plength = nx.path_weight(graph, dijkstra, 'weight')
                        # Check if path length is longer than the longest path length
                        if plength > length:
                            longest = dijkstra
                            length = plength

            vertices = [entries[0]]
            for n in range(1, len(longest)):
                vertices.append(shapely.Point(nodes[longest[n]]['coords']))

            spines.append(shapely.LineString(vertices))

        # Here, there is more than one entry,
        # This is the quickest, shortest path
        # are only calculated between
        # provided entry points
        else:
            targets = []
            for e in entries:
                targets.append(__get_closest_node(nodes, e))
            source = targets.pop(0)

            paths = []
            # Loop again through each node (target)
            for target in targets:
                if source != target:
                    # Check if the pair of nodes has a path
                    if nx.has_path(graph, source, target):
                        # Calculate the path
                        paths.append(nx.shortest_path(graph, source=source, target=target, weight='weight'))
            
            # reconstruct the lines
            if len(paths) > 0:
                # Here there is more than one path
                if len(paths) > 1:
                    overlapped = []
                    for ip, path in enumerate(paths):
                        vertices = [ entries[0] ]
                        for n in range(1, len(path) - 1):
                            vertices.append(shapely.Point(nodes[path[n]]['coords']))
                        vertices.append(entries[ip + 1])
                        overlapped.append(shapely.LineString(vertices))

                    union = unary_union(overlapped)
                    merged = linemerge(union.geoms)

                    if merged.geom_type == 'MultiLineString':
                        merged = [ x for x in merged.geoms ]
                    
                    spines.extend(merged)
                # Here there is one path between two entries
                else:
                    vertices = [ entries[0] ]
                    for n in range(1, len(paths[0]) - 1):
                        vertices.append(shapely.Point(nodes[paths[0][n]]['coords']))
                    vertices.append(entries[1])

                    spines.append(shapely.LineString(vertices))

    result = []
    for s in spines:
        if sigma is not None:
            smoothed = None
            if len(real_structural) > 0:
                smoothed = __keep_structural(s, real_structural, sigma)
            else:
                smoothed = gaussian_smoothing(s, sigma)
            result.append(smoothed)
        else:
            result.append(s)

    return result

def spinalize_polygons(polygons, densify=None, sigma=None, entries=None, structural=None):
    """
    Collapse multiple polygons into one or multiple lines.

    This algorithm creates the spine of the polygons to generate a line.
    When polygons touches each other, the generated spine is the union
    of all the touching polygons individual spine. Entry points are
    calculated as the middle of the line where polygon touches to
    avoid high computation time.

    More information about the usage of this algorithm in
    Touya & Girres :footcite:p:`touya:2013` :footcite:p:`touya:2014`

    Parameters
    ----------
    polygons : list of Polygon
        The polygons to spinalize.
    densify : float, optional
        The densification step used to resample the polygon ring.
    sigma : float, optional
        Gaussian filter strength.
    entries : list of Point, optional
        A list of entry points to the polygon.
        When two polygons touches, an entry point is created, and it is
        the middle of the LineString representing the intersection
        between those two polygons.
    structural : list of Point, optional
        A list of point inside the polygon that will be kept by the algorithm.

    Returns
    -------
    list of LineString

    Warning
    -------
    Please read the documentation of the :func:`spinalize_polygon` function
    because some parameters can cause extreme computation time.

    See Also
    --------
    spinalize_polygon:
        Collapse multiple polgygons into one or multiple lines.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> polygons = [Polygon([(0, 0), (0, 10), (20, 10), (20, 0), (0, 0)]), Polygon([(20, 0), (20, 10), (40, 10), (40, 0), (20, 0)])]
    >>> spinalize_polygons(polygons)
    [<LINESTRING (0 5, 10 5, 20 5, 30 5, 40 5)>]
    """
    def __merging(entries, spines):
        lines = []

        def __recursive_merging(spine, spines):
            degree = 0
            connected = None
            for i, s in enumerate(spines):
                if shapely.intersects(spine, s):
                    add = True
                    for l in lines:
                        if shapely.intersects(spine, l):
                            add = False
                    if add:
                        degree += 1
                        connected = i

            if degree == 1:
                connected = spines.pop(connected)
                merged = merge_linestrings(spine, connected)
                __recursive_merging(merged, spines)
            else:
                lines.append(spine)

        for entry in entries:
            index = None
            for i, s in enumerate(spines):
                if shapely.intersects(entry, s):
                    index = i

            spine = spines.pop(index)
            __recursive_merging(spine, spines)

            if len(spines) == 0:
                break

        return lines

    groups = group_intersecting(polygons)
    result = []

    for group in groups:
        spines = []
        for i, ip1 in enumerate(group):
            p1 = polygons[ip1]
            polyentries = []

            if entries is not None:
                for e in entries:
                    if shapely.dwithin(p1, e, 0.00001):
                        polyentries.append(e)

            for j, ip2 in enumerate(group):
                p2 = polygons[ip2]
                if i != j:
                    if shapely.intersects(p1, p2):
                        intersection = shapely.intersection(p1, p2)
                        if intersection.geom_type == 'LineString':
                            polyentries.append(get_line_middle_point(intersection))

            spines += spinalize_polygon(p1, densify, sigma=None, entries=polyentries)

        if len(spines) > 1:
            spines = __merging(entries, spines)

        if sigma is not None:
            if structural is not None:
                intersecting = []
                for p1 in group:
                    for s in structural:
                        if shapely.intersects(polygons[p1], s):
                            intersecting.append(s)
                if len(intersecting) > 0:
                    spines = [ __keep_structural(line, intersecting, sigma) for line in spines ]
                else:
                    spines = [ gaussian_smoothing(line, sigma) for line in spines ]
            else:
                spines = [ gaussian_smoothing(line, sigma) for line in spines ]
        
        result += spines

    return result


def __keep_structural(line, structural, sigma):
    coords = list(line.coords)

    breaks = []
    for s in structural:
        distance = shapely.distance(s, shapely.Point(coords[0]))
        closest = 0
        for i in range(1, len(coords)):
            c = coords[i]
            dist = shapely.distance(s, shapely.Point(c))
            if dist < distance:
                distance = dist
                closest = i
        breaks.append(closest)

    broken = []
    current = []
    for i in range(0, len(coords)):
        c = coords[i]
        current.append(c)
        if i in breaks:
            broken.append(current)
            current = [c]
    broken.append(current)

    smoothed = gaussian_smoothing(shapely.LineString(broken[0]), sigma)
    
    for ib in range(1, len(broken)):
        gaussian = gaussian_smoothing(shapely.LineString(broken[ib]), sigma)
        smoothed = merge_linestrings(smoothed, gaussian)

    return smoothed