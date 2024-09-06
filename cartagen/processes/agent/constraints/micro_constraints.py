import math
from shapely import Point
import numpy as np

from cartagen.processes.agent.constraints.abstract import GeneralisationConstraint
from cartagen.processes.agent.actions.generalisation_action import DeletionAction
from cartagen.processes.agent.actions.building_actions import EnlargementAction, EnlargeToRectangleAction, SimplificationAction, SquaringAction
from cartagen.utils.geometry.segment import get_segment_list
from cartagen.utils.geometry.line import geometry_flatten, geometry_len

class BuildingSizeConstraint(GeneralisationConstraint):
    """
    Agent onstraint to handle building size.

    This constraint asks for the building elimination if its area is
    too small for the wanted scale.

    Parameters
    ----------
    agent : BuildingAgent
        The building agent object to apply the constraint to.
    importance : int
        The Importance of the constraint.
    min_area : float, optional
        Area in square meters below which the building will be enlarge
        to this same value.
    elimination_area_threshold : float, optional
        Area in square meters below which the building will be eliminated.

    See Also
    --------
    run_agents:
        Execute the AGENT process.
    BuildingGranularityConstraint :
        Agent constraint to handle building granularity.
    BuildingSquarenessConstraint :
        Agent constraint to handle building squareness.
    """

    def __init__(self, agent, importance, min_area=0.0, area_threshold=70.0):
        self.agent = agent
        self.importance = importance
        self.building_min_area = min_area
        self.elimination_area_threshold = area_threshold
        self.type = "Size"

    # Compute the priority of the constraint given its current state.
    def compute_priority(self):
        self.priority = 10
        
    # Compute the current value of the constraint.
    def compute_current_value(self):
        if self.agent.deleted:
            self.current_value = 0.0
        geom = self.agent.feature['geometry']
        self.current_value = geom.area

    # Compute the goal value of the constraint.
    def compute_goal_value(self):
        if self.agent.deleted:
            self.goal_value = 0.0
        
        area = self.agent.feature['geometry'].area
        if(area < self.elimination_area_threshold):
            self.goal_value = 0.0
        elif(area > self.building_min_area):
            self.goal_value = area
        else:
            self.goal_value = self.building_min_area

    # Compute the satisfaction of the constraint according to the current and goal values.
    def compute_satisfaction(self):
        if(self.agent.deleted):
            self.satisfaction = 100.0
            return self.satisfaction
        
        self.compute_current_value()
        self.compute_goal_value()
        
        if(self.goal_value == 0.0):
            self.satisfaction = 0.0
        elif(self.current_value > self.goal_value):
            self.satisfaction = 100.0
        else:
            self.satisfaction = 100 - 100 * abs((self.current_value - self.goal_value) / self.goal_value)
        return self.satisfaction

    # Compute the actions of the constraint according to the current and goal values.
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
    """
    Agent constraint to handle building granularity.
    
    This constraint asks for the building simplification if its geometry
    is too detailed for the wanted scale.

    Parameters
    ----------
    agent : BuildingAgent
        The building agent object to apply the constraint to.
    importance : int
        The Importance of the constraint.
    min_length : float, optional
        The minimum segment lengths to consider the
        building too detailed.

    See Also
    --------
    run_agents:
        Execute the AGENT process.
    BuildingSizeConstraint :
        Agent constraint to handle building size.
    BuildingSquarenessConstraint :
        Agent constraint to handle building squareness.
    """

    def __init__(self, agent, importance, min_length=0.0):
        self.importance = importance
        self.agent = agent
        self.min_length_granularity = min_length
        self.type = "Granularity"

    # Compute the priority of the constraint given its current state.
    def compute_priority(self):
        self.priority = 8
        
    # Compute the current value of the constraint.
    def compute_current_value(self):
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
    
    # Compute the goal value of the constraint.
    def compute_goal_value(self):
        # in this case, do nothing
        return

    # Compute the satisfaction of the constraint according to the current and goal values.
    def compute_satisfaction(self):
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
    """
    Agent constraint to handle building squareness.
    
    This constraint asks for the building squaring if its angles
    are not exactly square.

    Parameters
    ----------
    agent : BuildingAgent
        The building agent object to apply the constraint to.
    importance : int
        The Importance of the constraint.
    angle_tolerance : float, optional
        The angle tolerance in degree to consider
        an angle right or flat.

    See Also
    --------
    run_agents:
        Execute the AGENT process.
    BuildingGranularityConstraint :
        Agent constraint to handle building granularity.
    BuildingSizeConstraint :
        Agent constraint to handle building size.
    """
    almostRightAnglesNumber = 0 # the number of angles that are almost right
    almostFlatAnglesNumber = 0 # the number of angles that are almost flat
    nbVertices = 0
    angle_tolerance = 15.0
    # threshold values used to compute satisfaction
    ALMOST_RIGHT_ANGLE_NUMBER_THRESHOLD = 6
    MINIMUM_THRESHOLD_TOLERANCE_ANGLE = 5
    NB_VERTICES_SIMPLE_POLYGONS = 6


    def __init__(self, agent, importance, angle_tolerance=15.0):
        self.importance = importance
        self.agent = agent
        self.angle_tolerance = angle_tolerance
        self.type = "Squareness"

    # Compute the priority of the constraint given its current state.
    def compute_priority(self):
        self.priority = 7
        
    # Compute the current value of the constraint.
    def compute_current_value(self):
        if self.agent.deleted:
            self.current_value = 0.0
        geom = self.agent.feature['geometry']
        self.nbVertices = sum(len(g.coords) for g in geometry_flatten(geom))
        self.almostRightAnglesNumber = 0
        self.almostFlatAnglesNumber = 0
        self.count_angles(geom.exterior)

        for interior in geom.interiors:
            self.count_angles(interior)

    # Compute the goal value of the constraint.
    def compute_goal_value(self):
        # in this case, do nothing
        return

    # Compute the satisfaction of the constraint according to the current and goal values.
    def compute_satisfaction(self):
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

    # Measure the given angle.
        
    # This method increment a flat or a right angle
    # if the provided angle is flat or right.

    # Parameters
    # ----------
    # point1 : Point
    #     First point to calculate the angle.
    # point2 : Point
    #     Second point to calculate the angle.
    # point3 : Point
    #     Third point to calculate the angle.
    def measure_angle(self, point1, point2, point3):
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