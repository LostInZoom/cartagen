# This is an implementation of the random building displacement algorithm
import random, math
import geopandas as gpd
import shapely
from cartagen.utils.partitioning.network import partition_networks

def random_displacement(
        polygons, networks=None, polygon_distance=10.0, network_distance=10.0,
        max_trials=25, max_displacement=10.0, network_partitioning=None
    ):
    """
    Iteratively displace polygons overlapping each other and the provided network.

    Displace the provided buildings if they overlap each other
    or are closer than the width value to the provided networks.
    This method is not deterministic, at each iteration, the algorithm
    select a random angle and a random distance to displace the polygon
    and see if it improves the overlapping issues. Running twice the
    algorithm returns two different solutions.
    This algorithm is succinctly presented in Wabiński *et al.*
    :footcite:p:`wabinski:2022`

    Parameters
    ----------
    polygons : GeoDataFrame of Polygon
        The buildings to displace.
    networks : list of GeoDataFrame of LineString, optional
        A list of networks the polygons need to be moved away from.
        If left to None, polygons will only be moved away from each other.
    polygon_distance : float, optional
        The minimum acceptable distance between polygons.
    network_distance : float, optional
        The minimum acceptable distance between the polygons
        and the provided networks.
    max_trials : int, optional
        The maximum number of trials before stopping the iteration.
        A trial represent the movement of one polygon that did not
        lower the mean overlapping area between polygons and networks.
    max_displacement : float, optional
        The maximum allowed distance of displacement per iteration.
    network_partioning : GeoDataFrame of LineString, optional
        The network to partition the data with. If provided, each network
        face is treated individually, thus improving performance on larger dataset.

    Returns
    -------
    GeoDataFrame

    See Also
    --------
    network_faces:
        Calculates the faces of one or multiple networks.

    References
    ----------
    .. footbibliography::
    """

    POLYGONS = polygons.to_dict('records')

    def __loop(indexes, network):
        """
        Launch a loop to iteratively displace polygons randomly
        and checks the congestion before validating
        """
        # Calculating the mean rate of overlapping polygons with other polygons and the provided network
        # It represents the mean congestion of polygons within the polygon block
        rate_mean = __get_buildings_overlapping_rate_mean(indexes, network)

        # Starting trial count
        trial = 0

        # Launching the loop which will displace buildings randomly as
        # long as the rate mean is above 0 and the max trial count is not exceeded
        while rate_mean > 0 and trial <= max_trials:
            # Selecting a random building index
            random_index = random.randint(0, len(indexes) - 1)
            random_building = indexes[random_index]

            overlap = __get_building_overlap(random_building, indexes, network)

            # Checking if that building is overlapping
            if overlap > 0:
                # Selecting a random angle (0-360)
                random_angle = random.uniform(0, 360)
                # Selecting a random length (0-max displacement variable)
                random_length = random.uniform(0, max_displacement)
                # Calculate displacement for x and y
                dx = math.cos(random_angle) * random_length
                dy = math.sin(random_angle) * random_length

                # Translating the building with the random values
                untranslated = POLYGONS[random_building]['geometry']
                translated = shapely.affinity.translate(untranslated, dx, dy)
                POLYGONS[random_building]['geometry'] = translated

                # Calulcating the new rate mean
                new_rate_mean = __get_buildings_overlapping_rate_mean(indexes, network)

                # If the new rate mean is equal or higher, we cancel the translation
                if (new_rate_mean >= rate_mean):
                    POLYGONS[random_building]['geometry'] = untranslated
                    trial += 1
                # Else, resetting the trial count and updating the rate mean
                else:
                    rate_mean = new_rate_mean
                    trial = 0

    def __get_buildings_overlapping_rate_mean(indexes, network):
        """
        Calculate the overlapping mean rate.
        """
        # If no polygons are provided, returns 0
        if len(indexes) == 0:
            return 0

        total = 0
        nb = 0
        # For each polygon, calculating the area overlapping other polygons,
        # and the network depending on a distance value
        for i in indexes:
            total += __get_building_overlap(i, indexes, network)
            nb += 1
        # Calculating the mean area
        if nb > 0:
            mean = total/nb

        return mean

    def __get_building_overlap(index, indexes, network):
        """
        Calculate the area of overlapping between polygons and the network
        """
        overlap = None
        b = POLYGONS[index]['geometry']
        buffered = b.buffer(polygon_distance)

        # For each buildings...
        for i in indexes:
            # Checking if it's not the same building
            if i != index:
                overlap = __get_overlapping_geometries(buffered, POLYGONS[i]['geometry'], overlap)
            
        # For each network section
        for n in network:
            overlap = __get_overlapping_geometries(b, n, overlap)

        # Returning the area of the geometry if it exists
        if (overlap is None) or (overlap.is_empty == True) or (overlap.area == 0):
            return 0
        else:
            return overlap.area

    def __get_overlapping_geometries(geom1, geom2, overlap):
        """
        Calculate the geometry of the intersection between a geographic object and a building
        """
        # If the building intersects the object
        if shapely.intersects(geom1, geom2):
            # Creating the intersection between the building and the object
            intersection = shapely.intersection(geom1, geom2)
            # If the geometry is empty, return the intersection...
            if overlap is None:
                return intersection
            # Else, returning the union between the intersection and the existing geometry
            else:
                return overlap.union(intersection)
        else:
            return overlap

    # Retrieve crs and convert polygons to list of dict records
    crs = polygons.crs

    network = []
    # Unpack the provided networks
    if networks is not None:
        for line in networks:
            # Append the buffered lines to the list
            for l in line.geometry:
                network.append(l.buffer(network_distance))

    if network_partitioning is not None:
        # Create the partitions -> tuple ([polygon index], [partition polygon geometry])
        partitions = partition_networks(polygons, *network_partitioning)

        # Loop through each partition
        for i, partition in enumerate(partitions[0]):
            partition_network = []

            if networks is not None:
                # Create the spatial index for the buffered network
                btree = shapely.STRtree(network)
                # Retrieve network sections that intersects the considered network face
                intersects = list(btree.query(partitions[1][i], predicate='intersects'))                
                # Append the buffered networks only if it truly intersects
                for inter in intersects:
                    if shapely.intersects(network[inter], partitions[1][i]):
                        partition_network.append(network[inter])

            # Launch the random displacement for the current partition
            __loop(partition, partition_network)
    else:
        # Launch the random displacement with all polygons and all buffered networks
        __loop(list(range(0, len(POLYGONS))), network)

    return gpd.GeoDataFrame(POLYGONS, crs=crs)