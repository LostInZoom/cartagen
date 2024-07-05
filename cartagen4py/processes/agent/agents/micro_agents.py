from cartagen4py.processes.agent.agents.abstract_agents import Agent

class BuildingAgent(Agent):

    def __init__(self, feature, importance=1):
        super().__init__(feature)
        self.importance = importance
        self.initial_geom = feature['geometry']
