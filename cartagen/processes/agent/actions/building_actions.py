from cartagen.processes.agent.actions.generalisation_action import GeneralisationAction
from cartagen.algorithms.buildings.simplification.ruas import simplify_building_ruas
from cartagen.algorithms.buildings.squaring.least_square import square_polygon_ls
from cartagen.utils.geometry.polygon import enclosing_rectangle
from shapely import affinity
from shapely.geometry import Polygon

class EnlargementAction(GeneralisationAction):
    """Building enlargement action."""
    def __init__(self, constraint, agent, weight, goal_area=0.0):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.goal_area = goal_area
        self.name = "Enlargement"

    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        geom = self.agent.feature['geometry']
        factor = self.goal_area / geom.area
        self.agent.feature['geometry'] = affinity.scale(geom, xfact= factor, origin=geom.centroid)

class EnlargeToRectangleAction(GeneralisationAction):
    """Building rectangle enlargement action."""
    def __init__(self, constraint, agent, weight, goal_area=0.0):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.goal_area = goal_area
        self.name = "Enlarge to rectangle"

    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        geom = self.agent.feature['geometry']
        ssr = enclosing_rectangle(geom, mode='input')
        factor = self.goal_area / ssr.area
        self.agent.feature['geometry'] = affinity.scale(ssr, xfact= factor, origin=ssr.centroid)

class SimplificationAction(GeneralisationAction):
    """Building simplification action."""
    def __init__(self, constraint, agent, weight, edge_threshold=0.0):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.edge_threshold = edge_threshold
        self.name = "Simplification"

    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        geom = self.agent.feature['geometry']
        new_geom = simplify_building_ruas(geom,self.edge_threshold)
        self.agent.feature['geometry'] = new_geom

class SquaringAction(GeneralisationAction):
    """Building squaring action."""
    def __init__(self, constraint, agent, weight):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.name = "Squaring"
    
    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        geom = self.agent.feature['geometry']
        new_geom = square_polygon_ls(geom)
        self.agent.feature['geometry'] = new_geom