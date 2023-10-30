from cartagen4py.processes.AGENT.agents.abstract_agents import MesoAgent
from shapely.ops import unary_union
from cartagen4py.data_enrichment.building_measures import *

class BlockAgent(MesoAgent):

    triangulation = None

    def __init__(self, feature, components, sections, importance=1):
        super().__init__(feature,components)
        self.importance = importance
        self.initial_geom = feature['geometry']
        self.sections = sections

    def get_simulated_density(self, building_min_size, road_sizes):
        """
        Gets the simulated density of the block taking into account the road symbol and the future size 
        of the buildings once they are generalised.
        ----------
        building_min_size : float
            The minimum size in m² of the buildings at the target scale.
        road_sizes : list of floats
            a lists of symbol widths for all the sections around the block. 
            The list uses the same order as the sections around the block.
        """

        if len(self.components) == 0:
            return 0.0
        
        building_area = 0.0
        for building in self.components:
            if building.deleted:
                continue
            # if the building is no longer inside the block, do not count it
            block_geom = self.feature['geometry']
            building_geom = building.feature['geometry']
            if block_geom.intersects(building_geom) == False:
                continue

            building_area += max([building_min_size,building_geom.area])

        # now compute the area of road symbols inside the block
        road_id = 0
        polygons = []
        for road in self.sections:
            symbol_width = road_sizes[road_id]
            polygons.append(road['geometry'].buffer(symbol_width))
            road_id += 1
        merged_symbols = unary_union(polygons)
        inter = merged_symbols.intersection(block_geom)

        return (inter.area + building_area)/block_geom.area
    
    def compute_block_triangulation(self):
        buildings = []
        roads = []
        for building in self.components:
            buildings.append(building['geometry'])
        for road in self.sections:
            roads.append(road['geometry'])
        self.triangulation = block_triangulation(buildings, roads, 30)
        return self.triangulation
