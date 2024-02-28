from cartagen4py.processes.AGENT.actions.generalisation_action import GeneralisationAction
from cartagen4py.algorithms.buildings.random_displacement import BuildingDisplacementRandom
from cartagen4py.algorithms.blocks.building_elimination import *
import geopandas as gpd

class RandomBlockDisplacementAction(GeneralisationAction):

    def __init__(self, constraint, agent, weight, section_symbols):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.section_symbols = section_symbols
        self.name = "RandomDisplacement"

    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        buffered_sections = [self.agent.sections[i].buffer(self.section_symbols[i]) for i in range(len(self.section_symbols))]
        roads_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(buffered_sections))
        components = [component.feature['geometry'] for component in self.agent.components]
        buildings_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(components))
        displacement = BuildingDisplacementRandom(network_partitioning=False)
        displaced_gdf = displacement.displace(buildings_gdf, roads_gdf)
        for i in range(len(components)):
            component = self.agent.components[i]
            component.feature['geometry'] = displaced_gdf.iloc(i)['geometry']

class PromBlockEliminationAction(GeneralisationAction):
    goal_area = 0.0

    def __init__(self, constraint, agent, weight, goal_area):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.goal_area = goal_area
        self.name = "PrometheeElimination"

    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        # TODO
        roads = [section.feature['geometry'] for section in self.agent.sections]
        buildings = [component.feature['geometry'] for component in self.agent.components]
        corners, corner_areas = corner_buildings(buildings, roads)
        triangulation = block_triangulation(buildings,roads,30.0)
        congestion = []
        for building in buildings:
            congestion.append(building_congestion(building,triangulation,30.0))
        elimination = building_elimination_ranking_in_block(buildings, triangulation, congestion, 250.0, corners)
        for eliminated in elimination:
            # check if we have to eliminate even more
            building = buildings[eliminated[0][0]]
