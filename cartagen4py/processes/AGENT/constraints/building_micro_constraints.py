from cartagen4py.processes.AGENT.constraints.generalisation_constraint import GeneralisationConstraint
from cartagen4py.processes.AGENT.actions.generalisation_action import DeletionAction
from cartagen4py.processes.AGENT.actions.enlargement_actions import EnlargementAction, EnlargeToRectangleAction
from cartagen4py.util.segment import *
import math

class BuildingSizeConstraint(GeneralisationConstraint):
    elimination_area_threshold = 70.0
    building_min_area = 0.0

    def __init__(self, importance, agent, building_min_area, elimination_area_threshold=70.0):
        self.importance = importance
        self.agent = agent
        self.building_min_area = building_min_area
        self.elimination_area_threshold = elimination_area_threshold

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

class BuildingGranularityConstraint(GeneralisationConstraint):
    min_length_granularity = 0 # the minimum length of edges in the building, expressed in meters

    def __init__(self, importance, agent, min_length_granularity):
        self.importance = importance
        self.agent = agent
        self.min_length_granularity = min_length_granularity

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
        
        nb = len(geom.coords) - 1
        tooShortEdgesRatio = len(small_segments)/nb
        self.granularity = math.sqrt(tooShortEdgesRatio*small_segment_deficit/self.min_length_granularity)

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
        # propose the simplification action
        # TODO
        action = EnlargementAction(self,self.agent,5,self.goal_value)
        self.actions.append([action, self, 5])
        # then propose the enlarge-to-rectangle action
        area = self.agent.feature['geometry'].area
        action = EnlargeToRectangleAction(self,self.agent,1,area)
        self.actions.append([action, self, 1])