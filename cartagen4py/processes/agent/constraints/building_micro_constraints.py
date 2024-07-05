from cartagen4py.processes.agent.constraints.generalisation_constraint import GeneralisationConstraint
from cartagen4py.processes.agent.actions.generalisation_action import DeletionAction
from cartagen4py.processes.agent.actions.building_actions import EnlargementAction, EnlargeToRectangleAction, SimplificationAction, SquaringAction
from cartagen4py.utils.geometry.segment import *
from cartagen4py.utils.geometry.line import geometry_flatten, geometry_len
import math
import numpy as np

class BuildingSizeConstraint(GeneralisationConstraint):
    elimination_area_threshold = 70.0
    building_min_area = 0.0

    def __init__(self, importance, agent, building_min_area, elimination_area_threshold=70.0):
        self.importance = importance
        self.agent = agent
        self.building_min_area = building_min_area
        self.elimination_area_threshold = elimination_area_threshold
        self.type = "Size"

    def compute_priority(self):
        """compute the priority of the constraint given its current state."""
        self.priority = 10
        
    def compute_current_value(self):
        """compute the current value of the constraint."""
        if self.agent.deleted:
            self.current_value = 0.0
        geom = self.agent.feature['geometry']
        self.current_value = geom.area

    def compute_goal_value(self):
        """compute the goal value of the constraint."""
        if self.agent.deleted:
            self.goal_value = 0.0
        
        area = self.agent.feature['geometry'].area
        if(area < self.elimination_area_threshold):
            self.goal_value = 0.0
        elif(area > self.building_min_area):
            self.goal_value = area
        else:
            self.goal_value = self.building_min_area

    def compute_satisfaction(self):
        """compute the satisfaction of the constraint according to the current and goal values."""
        if(self.agent.deleted):
            self.satisfaction = 100.0
            return
        
        self.compute_current_value()
        self.compute_goal_value()
        
        if(self.goal_value == 0.0):
            self.satisfaction = 0.0
        elif(self.current_value > self.goal_value):
            self.satisfaction = 100.0
        else:
            self.satisfaction = 100 - 100 * abs((self.current_value-self.goal_value)/self.goal_value)


    def compute_actions(self):
        self.actions.clear()
        if(self.goal_value == 0.0):
            # propose to delete the agent
            action = DeletionAction(self,self.agent,1)
            self.actions.append([action, self, 1])
        else:
            # propose the enlargement action
            action = EnlargementAction(self,self.agent,2,self.goal_value)
            self.actions.append([action, self, 2])
            # then propose the enlarge-to-rectangle action
            action = EnlargeToRectangleAction(self,self.agent,1,self.goal_value)
            self.actions.append([action, self, 1])

# A constraint on building granularity that asks for simplification if the building is too detailed for the scale
class BuildingGranularityConstraint(GeneralisationConstraint):
    min_length_granularity = 0 # the minimum length of edges in the building, expressed in meters

    def __init__(self, importance, agent, min_length_granularity):
        self.importance = importance
        self.agent = agent
        self.min_length_granularity = min_length_granularity
        self.type = "Granularity"

    def compute_priority(self):
        """compute the priority of the constraint given its current state."""
        self.priority = 8
        
    def compute_current_value(self):
        """compute the current value of the constraint."""
        if self.agent.deleted:
            self.current_value = 0.0
        geom = self.agent.feature['geometry']
        small_segments = []
        segments = get_segment_list(geom)
        for segment in segments:
            if(segment.length() < self.min_length_granularity):
                small_segments.append(segment)
        
        small_segment_deficit = 0
        if(len(small_segments)!=0):
            for segment in small_segments:
                small_segment_deficit += self.min_length_granularity - segment.length()
            small_segment_deficit = small_segment_deficit / len(small_segments)
        
        nb = geometry_len(geom) - 1
        tooShortEdgesRatio = len(small_segments)/nb
        self.current_value = math.sqrt(tooShortEdgesRatio*small_segment_deficit/self.min_length_granularity)

    def compute_goal_value(self):
        """compute the goal value of the constraint."""
        # in this case, do nothing
        return

    def compute_satisfaction(self):
        """compute the satisfaction of the constraint according to the current and goal values."""
        if(self.agent.deleted):
            self.satisfaction = 100.0
            return
        
        self.compute_current_value()
        
        if(self.current_value==0.0):
            self.satisfaction = 100.0
        else:
            self.satisfaction = round(100.0 - 100.0 * self.current_value,2)
            if(self.satisfaction < 0):
                self.satisfaction = 0.0


    def compute_actions(self):
        self.actions.clear()
        # propose the simplification action
        action = SimplificationAction(self,self.agent,5,self.min_length_granularity)
        self.actions.append([action, self, 5])
        # then propose the enlarge-to-rectangle action
        area = self.agent.feature['geometry'].area
        action = EnlargeToRectangleAction(self,self.agent,1,area)
        self.actions.append([action, self, 1])


# A constraint on building squareness that asks for squaring if the building angles are not exactly square
class BuildingSquarenessConstraint(GeneralisationConstraint):
    almostRightAnglesNumber = 0 # the number of angles that are almost right
    almostFlatAnglesNumber = 0 # the number of angles that are almost flat
    nbVertices = 0
    angle_tolerance = 15.0
    # threshold values used to compute satisfaction
    ALMOST_RIGHT_ANGLE_NUMBER_THRESHOLD = 6
    MINIMUM_THRESHOLD_TOLERANCE_ANGLE = 5
    NB_VERTICES_SIMPLE_POLYGONS = 6


    def __init__(self, importance, agent, angle_tolerance = 15.0):
        self.importance = importance
        self.agent = agent
        self.angle_tolerance = angle_tolerance
        self.type = "Squareness"

    def compute_priority(self):
        """compute the priority of the constraint given its current state."""
        self.priority = 7
        
    def compute_current_value(self):
        """compute the current value of the constraint."""
        if self.agent.deleted:
            self.current_value = 0.0
        geom = self.agent.feature['geometry']
        self.nbVertices = sum(len(g.coords) for g in geometry_flatten(geom))
        self.almostRightAnglesNumber = 0
        self.almostFlatAnglesNumber = 0
        self.count_angles(geom.exterior)

        for interior in geom.interiors:
            self.count_angles(interior)

    def compute_goal_value(self):
        """compute the goal value of the constraint."""
        # in this case, do nothing
        return

    def compute_satisfaction(self):
        """compute the satisfaction of the constraint according to the current and goal values."""
        if(self.agent.deleted):
            self.satisfaction = 100.0
            return
        
        self.compute_current_value()
        
        if(self.almostRightAnglesNumber==0.0):
            self.satisfaction = 100.0
        elif(self.almostRightAnglesNumber >= self.ALMOST_RIGHT_ANGLE_NUMBER_THRESHOLD):
            self.satisfaction = 0.0
        else:
            self.satisfaction = 100 - (100 / self.ALMOST_RIGHT_ANGLE_NUMBER_THRESHOLD)* self.almostRightAnglesNumber


    def compute_actions(self):
        self.actions.clear()
        # propose the squaring action
        action = SquaringAction(self,self.agent,5)
        self.actions.append([action, self, 5])

    def measure_angle(self,point1, point2, point3):
        # Calculation of the coordinate of point2-point1 vector
        x21 = point1.x-point2.x
        y21 = point1.y-point2.y
        # Calculation of the coordinate of point2-point3 vector
        x23 = point3.x-point2.x
        y23 = point3.y-point2.y
        # computation of the angle
        angle = abs(np.arctan2(x21 * y23 - y21 * x23, x21 * x23 + y21 * y23))

        # test if angle is almost flat
        diffFlat = abs(angle - math.pi) * 180 / math.pi
        if(diffFlat < self.angle_tolerance and diffFlat >= self.MINIMUM_THRESHOLD_TOLERANCE_ANGLE):
            self.almostFlatAnglesNumber+= 1
            return
        
        # test if angle is almost right
        difference = abs(angle - math.pi / 2) * 180 / math.pi
        if(difference <= self.angle_tolerance and difference >= self.MINIMUM_THRESHOLD_TOLERANCE_ANGLE):
            self.almostRightAnglesNumber+= 1
            return
        return
    
    def count_angles(self, ring):
        coords = ring.coords
        c0 = coords[0]
        c1 = coords[1]
        for i in range(2, len(coords)):
            c2 = coords[i]
            self.measure_angle(Point(c0), Point(c1), Point(c2))
            c0 = c1
            c1 = c2
        
        return