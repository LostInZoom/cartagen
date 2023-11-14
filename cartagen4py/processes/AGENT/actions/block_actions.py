from cartagen4py.processes.AGENT.actions.generalisation_action import GeneralisationAction
from cartagen4py.algorithms.buildings.random_displacement import BuildingDisplacementRandom
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