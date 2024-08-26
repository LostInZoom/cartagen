from cartagen.processes.agent.agents.abstract_agents import Agent

class BuildingAgent(Agent):
    """
    Agent for buildings generalisation.

    Create an object to generalise buildings within the AGENT process.

    Parameters
    ----------
    feature : GeoSeries of Polygon
        The building to create the agent from.
    importance : int, optional
        Importance of the agent within the process.

    See Also
    --------
    run_agents:
        Execute the AGENT process.
    """
    def __init__(self, feature, importance=1):
        super().__init__(feature)
        self.importance = importance
        self.initial_geom = feature['geometry']
