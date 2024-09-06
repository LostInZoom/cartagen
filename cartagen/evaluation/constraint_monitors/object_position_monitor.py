from evaluation.constraint_monitors.constraint_monitor import ConstraintMonitor
from shapely.geometry import Point, BaseGeometry
from evaluation.constraint_satisfaction import ConstraintSatisfaction

class BuildingPositionMonitor(ConstraintMonitor):

    epsilon = 1.0

    # the goal value should be an area converted in meters, not in map mm, as it is usually defined
    # if the building is deleted by generalisation, the building variable used in this constructor should be None
    def __init__(self, goal_value, object: BaseGeometry, object_ini: BaseGeometry, importance: int):
        self.goal_value = goal_value
        self.object = object
        self.object_ini = object_ini
        self.importance = importance

    def get_geometry(self) -> Point:
        """Get the point geometry locating the constraint monitor."""
        return self.object.centroid

    def get_importance(self) -> int:
        """Get the importance of the constraint monitor."""
        return self.importance

    def get_satisfaction(self) -> ConstraintSatisfaction:
        """Get the current satisfaction of the constraint."""
        satisfaction = ConstraintSatisfaction.PERFECT;
        if(self.object is None):
            return ConstraintSatisfaction.VERY_SATISFIED;

        current_value = self.compute_current_value()
        
        if(current_value < self.epsilon):
            satisfaction = ConstraintSatisfaction.PERFECT;
        elif(current_value< self.goal_value/4):
            satisfaction = ConstraintSatisfaction.VERY_SATISFIED;
        elif(current_value < self.goal_value/2):
            satisfaction = ConstraintSatisfaction.CORRECT;
        elif(current_value <= self.goal_value):
            satisfaction = ConstraintSatisfaction.ACCEPTABLE;
        elif(current_value < self.goal_value * 1.5):
            satisfaction = ConstraintSatisfaction.FAIR;
        elif(current_value < self.goal_value * 2):
            satisfaction = ConstraintSatisfaction.BARELY_SATISFIED;
        elif(current_value < self.goal_value * 2.5):
            satisfaction = ConstraintSatisfaction.NOT_SATISFIED;   
        else:
            satisfaction = ConstraintSatisfaction.UNACCEPTABLE;            

        return satisfaction

    def get_subject(self):
        """Get the subject of the constraint monitor: it can be a geometry, or a group of geometries."""
        return self.object
    
    def get_distance(self, geometry1: BaseGeometry, geometry2: BaseGeometry):
        if(geometry1.type() == "Point"):
            return geometry1.distance(geometry2)
        elif(geometry1.type() == "LineString"):
            return geometry1.hausdorff_distance(geometry2)
        else:
            return geometry1.centroid().distance(geometry2.centroid())
    
    def compute_current_value(self):
        dist = self.get_distance(self.object_ini,self.object)
        return dist