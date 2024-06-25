# This is an implementation of the random building displacement algorithm
import random, math, numpy
import geopandas as gpd
import shapely
from cartagen4py.utils.partitioning.network import network_partition

def random_displacement(
        polygons, networks=None, polygon_distance=10, network_width=10,
        max_trials=25, max_displacement=10, network_partitioning=None
    ):
    """
    Iteratively displace polygons overlapping each other and the provided network.

    Displace the provided buildings if they overlap each other or are closer than the width value to the provided networks.

    Parameters
    ----------
    buildings : geopandas GeoDataFrame of Polygons.
        The buildings to displace.
    width : float.
        The width of the provided networks. A buffer is applied of the given width to the network
        and then is used to displace the buildings when overlapping.
    networks : geopandas GeoDataFrame of LineStrings, optional.
        The networks which will be used to displace the buildings.
        If no networks is provided, the building will only be moved if they intersects each other.

    Returns
    -------
    geopandas GeoDataFrame of Polygons
        The displaced buildings.
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

            # Checking if that building is overlapping
            if (__get_building_overlap(random_building, indexes, network) != 0):
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
        Calculate the overlapping mean rate
        """
        # If no polygons are provided, returns 0
        if len(indexes) == 0:
            return 0
        mean = 0
        nb = 0
        # For each polygon, calculating the area overlapping other polygons,
        # and the network depending on a distance value
        for i in indexes:
            mean += __get_building_overlap(i, indexes, network)
            nb += 1
        # Calculating the mean area
        if nb > 0:
            mean = mean/nb

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
            if i != indexes:
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

    buffer = []
    # Unpack the provided networks
    if networks is not None:
        for line in networks:
            # Append the buffered lines to the list
            for l in line.geometry:
                buffer.append(l.buffer(network_width))

    if network_partitioning is not None:
        # Create the partitions -> tuple ([polygon index], [partition polygon geometry])
        partitions = network_partition(polygons, *network_partitioning)

        # Loop through each partition
        for i, partition in enumerate(partitions[0]):
            partition_buffer = []

            if networks is not None:
                # Create the spatial index for the buffered network
                btree = shapely.STRtree(buffer)
                # Retrieve network sections that intersects the considered network face
                intersects = list(btree.query(partitions[1][i], predicate='intersects'))                
                # Append the buffered networks only if it truly intersects
                for inter in intersects:
                    if shapely.intersects(buffer[inter], partitions[1][i]):
                        partition_buffer.append(buffer[inter])

            # Launch the random displacement for the current partition
            __loop(partition, partition_buffer)
    else:
        # Launch the random displacement with all polygons and all buffered networks
        __loop(list(range(0, len(POLYGONS))), buffer)

    return gpd.GeoDataFrame(POLYGONS, crs=crs)


# class RandomDisplacement:
#     """
#     Iteratively displace buildings that overlap with each other and/or a provided network.

#     Parameters
#     ----------
#     max_trials : int optional
#         When a building intersects an other building or the network, a random displacement is apply. If the building still
#         intersects an other building or the network, this represents a trial.
#         Default to 25. This means a building will be displaced 25 times and if no solution has been found, the building is not moved.
#     max_displacement : float optional
#         The maximal displacement in meters allowed per iteration.
#         Default to 10.
#     network_partitioning : list of geopandas GeoDataFrame of LineStrings
#         A list of GeoDataFrame representing the networks which will be used for the partitioning.
#         Default value is set to None which doesn't apply any network partitioning.
    
#     """
#     def __init__(self, max_trials=25, max_displacement=10, network_partitioning=None):
#         self.MAX_TRIALS = max_trials
#         self.MAX_DISPLACEMENT = max_displacement
#         self.NETWORK_PARTITIONING = network_partitioning

#     def displace(self, buildings, width, *networks):
#         """
        
#         """  