from cartagen4py.processes.AGENT.actions.generalisation_action import GeneralisationAction
from shapely import affinity

class EnlargementAction(GeneralisationAction):
    goal_area = 0.0

    def __init__(self, constraint, agent, weight, goal_area):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.goal_area = goal_area

    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        geom = self.agent.feature['geometry']
        factor = self.goal_area / geom.area
        self.agent.feature['geometry'] = affinity.scale(geom, xfact= factor, origin=geom.centroid)

class EnlargeToRectangleAction(GeneralisationAction):
    goal_area = 0.0

    def __init__(self, constraint, agent, weight, goal_area):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.goal_area = goal_area

    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        geom = self.agent.feature['geometry']
        ssr = geom.minimum_rotated_rectangle
        factor = self.goal_area / ssr.area
        self.agent.feature['geometry'] = affinity.scale(ssr, xfact= factor, origin=ssr.centroid)