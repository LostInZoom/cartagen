from evaluation.constraint_monitors.constraint_monitor import ConstraintMonitor
from shapely.geometry import Point, Polygon
from evaluation.constraint_satisfaction import ConstraintSatisfaction

class BuildingSizeMonitor(ConstraintMonitor):

    epsilon = 5.0

    # the goal value should be an area converted in meters, not in map mm, as it is usually defined
    # if the building is deleted by generalisation, the building variable used in this constructor should be None
    def __init__(self, goal_value, building: Polygon, building_ini: Polygon, importance: int):
        self.target_area = goal_value
        self.building = building
        self.building_ini = building_ini
        self.goal_value = goal_value
        self.importance = importance

    def get_geometry(self) -> Point:
        """Get the point geometry locating the constraint monitor."""
        return self.building.centroid

    def get_importance(self) -> int:
        """Get the importance of the constraint monitor."""
        return self.importance

    def get_satisfaction(self) -> ConstraintSatisfaction:
        """Get the current satisfaction of the constraint."""
        satisfaction = ConstraintSatisfaction.PERFECT;
        if(self.building is None):
            return ConstraintSatisfaction.VERY_SATISFIED;

        current_value = self.building.area;
        self.compute_goal_value()

        if(abs(self.goal_value-current_value) < self.epsilon):
            satisfaction = ConstraintSatisfaction.PERFECT;
        elif(current_value>self.goal_value):
            satisfaction = ConstraintSatisfaction.VERY_SATISFIED;
        elif(current_value > 5 * self.goal_value / 6):
            satisfaction = ConstraintSatisfaction.CORRECT;
        elif(current_value > 3 * self.goal_value / 4):
            satisfaction = ConstraintSatisfaction.ACCEPTABLE;
        elif(current_value>self.goal_value / 2):
            satisfaction = ConstraintSatisfaction.FAIR;
        elif(current_value>self.goal_value / 3):
            satisfaction = ConstraintSatisfaction.BARELY_SATISFIED;
        elif(current_value>self.goal_value / 4):
            satisfaction = ConstraintSatisfaction.NOT_SATISFIED;   
        else:
            satisfaction = ConstraintSatisfaction.UNACCEPTABLE;            

        return satisfaction

    def get_subject(self):
        """Get the subject of the constraint monitor: it can be a geometry, or a group of geometries."""
        return self.building

    def compute_goal_value(self):
        initial_value = self.building_ini.area
        if(initial_value < self.target_area):
            self.goal_value=self.target_area
        elif (initial_value < 2 * self.target_area):
            a = (1.0 - 2.0 * self.epsilon) / self.target_area
            b = self.target_area - 1.0 + 3.0 * self.epsilon
            self.goal_value=a * initial_value + b
        else:
            self.goal_value=initial_value