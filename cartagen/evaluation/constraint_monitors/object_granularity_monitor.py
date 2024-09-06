from evaluation.constraint_monitors.constraint_monitor import ConstraintMonitor
from shapely.geometry import Point, Polygon
from evaluation.constraint_satisfaction import ConstraintSatisfaction
from cartagen.utils.geometry.line import get_shortest_edge_length

class ObjectGranularityMonitor(ConstraintMonitor):

    epsilon = 0.05

    # the goal value should be an area converted in meters, not in map mm, as it is usually defined
    # if the object is deleted by generalisation, the object variable used in this constructor should be None
    # this constraint is valid for line or polygon objects
    def __init__(self, goal_value, object: Polygon, object_ini: Polygon, importance: int):
        self.goal_value = goal_value
        self.object = object
        self.object_ini = object_ini
        self.importance = importance
        self.current_value = (0.0,0.0)

    def get_geometry(self) -> Point:
        """Get the point geometry locating the constraint monitor."""
        return self.object.centroid

    def get_importance(self) -> int:
        """Get the importance of the constraint monitor."""
        return self.importance

    def get_satisfaction(self) -> ConstraintSatisfaction:
        """Get the current satisfaction of the constraint."""
        satisfaction = ConstraintSatisfaction.PERFECT;
        if(self.building is None):
            return ConstraintSatisfaction.VERY_SATISFIED;

        current_value = self.compute_current_value()

        if(current_value[0] >= self.goal_value[0]):
            satisfaction = ConstraintSatisfaction.PERFECT;
        elif((current_value[0] > 3 * self.goal_value[0] / 4) and (current_value[1] < self.epsilon)):
            satisfaction = ConstraintSatisfaction.VERY_SATISFIED;
        elif(current_value[1] < self.epsilon):
            satisfaction = ConstraintSatisfaction.CORRECT;
        elif(current_value[0] > 3 * self.goal_value[0] / 4) and (current_value[1] < 2 * self.epsilon):
            satisfaction = ConstraintSatisfaction.ACCEPTABLE;
        elif(current_value[0] > self.goal_value[0] / 2) and (current_value[1] < 2 * self.epsilon):
            satisfaction = ConstraintSatisfaction.FAIR;
        elif(current_value[0] > self.goal_value[0] / 2) and (current_value[1] > 2 * self.epsilon):
            satisfaction = ConstraintSatisfaction.BARELY_SATISFIED;
        elif(current_value[0] > self.goal_value[0] / 4) and (current_value[1] > 2 * self.epsilon):
            satisfaction = ConstraintSatisfaction.NOT_SATISFIED;   
        else:
            satisfaction = ConstraintSatisfaction.UNACCEPTABLE;            

        return satisfaction

    def get_subject(self):
        """Get the subject of the constraint monitor: it can be a geometry, or a group of geometries."""
        return self.object

    def compute_current_value(self):
        vertex_density = 0.0
        min_length = 0.0
        if (self.object.geom_type == "Polygon"):
            minLength = get_shortest_edge_length(self.object.exterior)
        elif(self.object.geom_type == "LineString"):
            minLength = get_shortest_edge_length(self.object)
        
        nb_vertices = len(self.object.coords)
        length = self.object.length
        vertex_density = nb_vertices/length
        self.current_value = (min_length, vertex_density)

        
        