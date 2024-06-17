# This is an implementation of the random building displacement algorithm
import random, math, numpy
import geopandas as gpd
import shapely
from cartagen4py.utils.partitioning.network import network_partition

class RandomDisplacement:
    """
    Iteratively displace buildings that overlap with each other and/or a provided network.

    Parameters
    ----------
    max_trials : int optional
        When a building intersects an other building or the network, a random displacement is apply. If the building still
        intersects an other building or the network, this represents a trial.
        Default to 25. This means a building will be displaced 25 times and if no solution has been found, the building is not moved.
    max_displacement : float optional
        The maximal displacement in meters allowed per iteration.
        Default to 10.
    network_partitioning : list of geopandas GeoDataFrame of LineStrings
        A list of GeoDataFrame representing the networks which will be used for the partitioning.
        Default value is set to None which doesn't apply any network partitioning.
    
    """
    def __init__(self, min_dist, max_trials=25, max_displacement=10, network_partitioning=None):
        self.MAX_TRIALS = max_trials
        self.MAX_DISPLACEMENT = max_displacement
        self.NETWORK_PARTITIONING = network_partitioning
        self.MIN_DIST = min_dist

    def displace(self, buildings, width, *networks):
        """
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

        # Retrieve crs and convert buildings to list of dict records
        crs = buildings.crs
        records = buildings.to_dict('records')
        self.__BUILDINGS = records

        buffer = []
        # Unpack the provided networks
        for line in networks:
            # Append the buffered lines to the list
            for l in line.geometry:
                buffer.append(l.buffer(width))

        if self.NETWORK_PARTITIONING is not None:
            # Create the partitions -> tuple ([buildings index], [partition polygon geometry])
            partitions = network_partition(buildings, *self.NETWORK_PARTITIONING)

            # Loop through each partition
            for i, partition in enumerate(partitions[0]):
                # Create the spatial index for the buffered network
                btree = shapely.STRtree(buffer)
                # Retrieve network sections that intersects the considered network face
                intersects = list(btree.query(partitions[1][i], predicate='intersects'))

                partition_buffer = []
                # Append the buffered networks only if it truly intersects
                for inter in intersects:
                    if shapely.intersects(buffer[inter], partitions[1][i]):
                        partition_buffer.append(buffer[inter])

                # Launch the random displacement for the current partition
                self.__random_displacement(partition, partition_buffer)
        else:
            # Launch the random displacement with all buildings and all buffered networks
            self.__random_displacement(list(range(0, len(records))), buffer)

        return gpd.GeoDataFrame(records, crs=crs)

    def __random_displacement(self, buildings, buffer):
        """
        Launch a loop to iteratively displace buildings randomly and checks the congestion before validating
        """
        # Calculating the mean rate of overlapping buildings with other buildings and the provided network
        # It represents the mean congestion of buildings within the building block
        rate_mean = self.__get_buildings_overlapping_rate_mean(buildings, buffer)

        # Starting trial count
        trial = 0

        # Launching the loop which will displace buildings randomly as long as the rate mean is above 0 and the max trial count is not exceeded
        while rate_mean > 0 and trial <= self.MAX_TRIALS:
            # Selecting a random building index
            random_building = buildings[random.randint(0, len(buildings) - 1)]

            # Checking if that building is overlapping
            if (self.__get_building_overlap(random_building, buildings, buffer) != 0):
                # Selecting a random angle (0-360)
                random_angle = random.uniform(0, 360)
                # Selecting a random length (0-max displacement variable)
                random_length = random.uniform(0, self.MAX_DISPLACEMENT)
                # Calculate displacement for x and y
                dx = math.cos(random_angle) * random_length
                dy = math.sin(random_angle) * random_length

                # Translating the building with the random values
                untranslated = self.__BUILDINGS[random_building]['geometry']
                translated = shapely.affinity.translate(untranslated, dx, dy)
                self.__BUILDINGS[random_building]['geometry'] = translated

                # Calulcating the new rate mean
                new_rate_mean = self.__get_buildings_overlapping_rate_mean(buildings, buffer)

                # If the new rate mean is equal or higher, we cancel the translation
                if (new_rate_mean >= rate_mean):
                    self.__BUILDINGS[random_building]['geometry'] = untranslated
                    trial += 1
                # Else, resetting the trial count and updating the rate mean
                else:
                    rate_mean = new_rate_mean
                    trial = 0


    def __get_buildings_overlapping_rate_mean(self, buildings, buffer):
        """
        Calculate the overlapping mean rate
        """
        # If no buildings are provided, returns 0
        if len(buildings) == 0:
            return 0
        mean = 0
        nb = 0

        # For each building, calculating the area overlapping other buildings, roads and rivers depending on a distance value
        for i in buildings:
            mean += self.__get_building_overlap(i, buildings, buffer)
            nb += 1

        # Calculating the mean area
        if nb > 0:
            mean = mean/nb

        return mean

    def __get_building_overlap(self, index, buildings, buffer):
        """
        Calculate the area of overlapping between buildings, roads and rivers
        """
        geometry = None
        b = self.__BUILDINGS.iloc[processed_building].geometry
        buffered = b.buffer(self.MIN_DIST)
        # For each buildings...
        for building in buildings:
            # Checking if it's not the same building
            if building != processed_building:
                geometry = self.__get_overlapping_geometries(buffered, self.__BUILDINGS.iloc[building].geometry, geometry, "building")

        # For each road section
        for road in roads:
            geometry = self.__get_overlapping_geometries(buffered, road, geometry, "road")

        # For each river section
        for river in rivers:
            geometry = self.__get_overlapping_geometries(buffered, river, geometry, "river")

        # Returning the area of the geometry if it exists
        if (geometry is None) or (geometry.is_empty == True) or (geometry.area == 0):
            return 0
        else:
            return geometry.area

    def __get_overlapping_geometries(self, geom1, geom2, overlap):
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

    def __gdf_to_geomlist(self, gdf):
        l = []
        for geom in gdf.geometry:
            l.append(geom)
        return l
