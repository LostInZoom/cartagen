from random import sample
import shapely
import geopandas as gpd
import networkx as nx

from cartagen.utils.network.graph import create_graph

def rural_betweeness(roads, sample_size=10, threshold=0, cost=None):
    """
    Detect central roads inside a network using betweeness centrality.

    This algorithm is used to detect central road inside the network.
    It is primarily used to eliminate less used roads in rural areas.
    The centrality of the road sections are determiined using the 
    betweeness centrality.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        The road network to analyze.
    sample_size : int, optional
        The percentage of nodes to keep for the betweenness calculation.
        Higher values means better results at the
        cost of exponential computation length.
    threshold : int, optional
        The minimum betweenness value for the road to be kept.
    cost : str, optional
        The name of the attribute giving the cost of the road section.
        Make sure the attribute is a number.
        Default to None, which means the length of the road is used as the cost.

    Returns
    -------
    GeoDataFrame of LineString

    See Also
    --------
    rural_traffic :
        Detect central roads inside a network using traffic simulation.

    Notes
    -----
    This algorithm is based on Touya :footcite:p:`touya:2010` but
    the traffic simulation is replaced by the calculation of
    the betweeness centrality, which is less resource-intensive.
    However, it is harder to find the right betweeness threshold to use,
    and some road sections might become unconnected to the
    rest of the network.

    References
    ----------
    .. footbibliography::
    """
    # Retrieve crs for output
    crs = roads.crs
    network = roads.to_dict('records')

    # If cost has been set
    if cost is not None:
        # Raise an exception if attribute does not exists
        if hasattr(roads, cost) == False:
            raise Exception('Selected cost attribute does not exists.')

    # Create the graph from the road network
    graph = create_graph(roads, cost)

    # Determinate the number of nodes to compute the centrality index
    k = round((len(list(graph.nodes)) * sample_size) / 100)

    # Calculate the edge betweenness centrality
    centrality = nx.edge_betweenness_centrality(graph, k=k, weight='weight')

    # Assign the value of betweenness for each road
    roads_centrality = [0] * len(network)
    for edge, value in centrality.items():
        roads_centrality[graph.edges[edge]['rid']] = value

    result = []
    # Loop through the roads
    for rid, road in enumerate(network):
        # Get the centrality of the road
        b = roads_centrality[rid]

        # Define the state
        state = False
        # If the associated traffic is above the min_traffic threshold
        if b > threshold:
            # Set the state to True
            state = True
        
        # Set attributes to the road
        road['rural'] = state
        road['betweeness'] = b

        result.append(road)

    return gpd.GeoDataFrame(result, crs=crs)
    

def rural_traffic(roads, min_traffic=1,
        attraction_points=None, max_distance=None,
        sample_size=80, export_samples=False, cost=None
    ):
    """
    Detect central roads inside a network using traffic simulation.

    This algorithm proposed by Touya :footcite:p:`touya:2010`
    detects the most used road inside a network by simulating
    traffic between attraction points. Thos points can be provided
    or they can randomly be extracted from the network.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        The road network to analyze.
    min_traffic : int, optional
        The minimum number of time a road must be used to be kept.
    attraction_points : GeoDataFrame of Point, optional
        The attraction points between which traffic will be calculated.
        If provided, the points doesn't need to be on the network, they
        will be snapped to the closest node (*i.e.* intersection) on the network.
        Default to None, meaning a random set a nodes will be taken.
    max_distance : int, optional
        If attraction points are provided, it is the maximum distance
        from which the provided attraction points are snapped to the road network nodes.
        If further from this distance, the point will not be used.
        Default to None, meaning the point is snapped no matter its distance from the network.
    sample_size : int, optional
        If no attraction points are provided, it is the number
        of nodes to take into account when calculating traffic.
        This means the calculation between all pairs of nodes, this
        can be really resource-intensive.
    export_samples : bool, optional
        If set to True, returns a tuple of two GeoDataFrame instead of a single GeoDataFrame. The second
        is a point GeoDataFrame of the samples used for the traffic calculation. If attraction points
        were provided, returns the snapped points on the nodes, else returns the random node sample.
    cost : str, optional
        The name of the attribute giving the cost of the road section. Make sure the attribute is a number.
        Default to None, which means the length of the road is used as the cost.

    Returns
    -------
    roads : GeoDataFrame of LineString
        The network with simulated traffic information.
    samples : GeoDataFrame of Point, optional
        The random sample if wanted, and if no attraction points were provided.

    Warning
    -------
    This algorithm can quickly become very intensive, especially
    when using a large dataset and a high ``sample_size`` parameter,
    but it preserves the connectivity of the network.

    See Also
    --------
    rural_betweeness :
        Detect central roads inside a network using traffic simulation.

    References
    ----------
    .. footbibliography::
    """

    # Retrieve crs for output
    crs = roads.crs

    # If cost has been set
    if cost is not None:
        # Raise an exception if attribute does not exists
        if hasattr(roads, cost) == False:
            raise Exception('Selected cost attribute does not exists.')            

    # Create the graph from the road network
    graph = create_graph(roads, cost)

    # Convert geodataframe to list of dicts
    roads = roads.to_dict('records')

    # Create a counter for the number of time a road is used
    count = [0] * len(roads)

    # If no attraction were provided
    if attraction_points is None:
        # Formula to convert the sample_size from percentage to size depending on nodes number
        # nb_sample = round((len(list(graph.nodes)) * sample_size) / 100)
        
        # Select random targets in the graph node
        targets = sample(list(graph.nodes), sample_size)
    else:
        # Convert attraction points to a list of records
        points = attraction_points.to_dict('records')

        # Create a list of all the nodes as points
        ncoordinates = []
        for n in graph.nodes:
            ncoordinates.append(shapely.Point(graph.nodes[n]['coords']))

        # Calculate the tree
        tree = shapely.STRtree(ncoordinates)

        # Loop through attraction points
        attraction = []
        for p in points:
            # Get its geometry
            attraction.append(p['geometry'])

        # Retrieve the closest node from the provided points
        if max_distance is None:
            nearest = tree.query_nearest(attraction, all_matches=False)
        else:
            nearest = tree.query_nearest(attraction, max_distance=max_distance, all_matches=False)

        # Set the list of node indexes as the targets
        targets = list(nearest[1])

    if len(targets) < 2:
        raise Exception("Too few or too far attraction points to calculate traffic.")

    # Loop through all targets to get the source
    for source in targets:
        # Loop through all targets to get the target
        for target in targets:
            # Make sure source and target is different
            if source != target:
                # Check if the pair of nodes has a path
                if nx.has_path(graph, source, target):
                    # Calculate the path
                    path = nx.shortest_path(graph, source=source, target=target, weight='weight')
                    # Loop through each node used
                    for inode, node in enumerate(path):
                        # Avoid last node
                        if inode < len(path) - 1:
                            # Get source and destination of the processed edge
                            src, dst = node, path[inode + 1]
                            # Check this edge exists
                            if graph.has_edge(src, dst):
                                # Increase the counter of the associated road
                                count[graph.edges[src, dst]['rid']] += 1

    result = []
    # Loop through the roads
    for rid, road in enumerate(roads):
        # Get the traffic of the road
        traffic = count[rid]

        # Define the state
        state = False
        # If the associated traffic is above the min_traffic threshold
        if traffic > min_traffic:
            # Set the state to True
            state = True
        
        # Set attributes to the road
        road['rural'] = state
        road['traffic'] = traffic

        # Add it to the results
        result.append(road)

    # If export samples is selected
    if export_samples:
        # Retrieve each source/target nodes
        samples = []
        for tid, target in enumerate(targets):
            samples.append({
                'tid': tid,
                'geometry': shapely.Point(graph.nodes[target]['coords'])
            })

        rsamples = None
        if len(samples) > 0:
            rsamples = gpd.GeoDataFrame(samples, crs=crs)
        if len(result) > 0:
            return gpd.GeoDataFrame(result, crs=crs), rsamples
        else:
            return None, rsamples
    else:
        if len(result) > 0:
            return gpd.GeoDataFrame(result, crs=crs)
        else:
            return gpd.GeoDataFrame()