from cartagen.processes.agent.constraints.abstract import GeneralisationConstraint
from cartagen.processes.agent.actions.generalisation_action import MesoMicroActivation
from cartagen.processes.agent.actions.block_actions import PromBlockEliminationAction, RandomBlockDisplacementAction

class ComponentsSatisfactionConstraint(GeneralisationConstraint):

    def __init__(self, importance, agent):
        self.importance = importance
        self.agent = agent
        self.type = "ComponentsSatisfaction"

    def compute_priority(self):
        """compute the priority of the constraint given its current state."""
        self.priority = 10
        
    def compute_current_value(self):
        """compute the current value of the constraint."""
        self.current_value = self.agent.get_components_satisfaction()

    def compute_goal_value(self):
        """compute the goal value of the constraint."""
        # do nothing
        return

    def compute_satisfaction(self):
        """compute the satisfaction of the constraint according to the current and goal values."""
        if(self.agent.deleted):
            self.satisfaction = 100.0
            return
        
        self.compute_current_value()
        # the satisfaction is directly the components statisfaction value
        self.satisfaction = self.current_value

    def compute_actions(self):
        self.actions.clear()
        # propose the action that activates the micro agents
        action = MesoMicroActivation(self,self.agent,1)
        self.actions.append([action, self, 1])

class BlockDensityConstraint(GeneralisationConstraint):
    build_min_size = 200
    density_ratio = 1.0
    road_sizes = []
    goal_density = -1.0
    initial_density = -1.0
    
    def __init__(self, importance, agent, build_min_size, density_ratio, road_sizes):
        self.importance = importance
        self.agent = agent
        self.density_ratio = density_ratio
        self.build_min_size = build_min_size
        self.road_sizes = road_sizes
        self.type = "BlockDensity"

    def compute_priority(self):
        """compute the priority of the constraint given its current state."""
        self.priority = 9
        
    def compute_current_value(self):
        """compute the current value of the constraint."""
        self.current_value = self.agent.get_simulated_density(self.build_min_size, self.road_sizes)

    def compute_goal_value(self):
        """compute the goal value of the constraint."""
        self.initial_density = self.agent.get_initial_density(self.road_sizes)
        self.goal_density = self.density_ratio * self.initial_density
        if(self.goal_density > 1):
            self.goal_density = 1.0
        return

    def compute_satisfaction(self):
        """compute the satisfaction of the constraint according to the current and goal values."""
        if(self.agent.deleted):
            self.satisfaction = 100.0
            return

        self.compute_current_value()
        if self.goal_density == -1.0:
            self.compute_goal_value()

        # the satisfaction is directly the components statisfaction value
        density = abs(self.current_value - self.goal_density)
        self.satisfaction = 100 - round(density/0.003)
        if self.satisfaction < 0:
            self.satisfaction = 0.0

    def compute_actions(self):
        self.actions.clear()
        if self.current_value > self.initial_density:
            nb_elim = 1
            if self.satisfaction < 75 and len(self.agent.components) <= 10:
                nb_elim = 3
            elif self.satisfaction < 75 and len(self.agent.components) < 15:
                nb_elim = 4
            elif self.satisfaction < 75:
                nb_elim = 6
            elif self.satisfaction < 85:
                nb_elim = 2
            action = PromBlockEliminationAction(self,self.agent,2, nb_elim)
            self.actions.append([action, self, 2])

class BlockProximityConstraint(GeneralisationConstraint):
    min_sep = 0.0
    road_sizes = []

    def __init__(self, importance, agent, min_sep, road_sizes):
        self.importance = importance
        self.agent = agent
        self.type = "BlockProximity"
        self.min_sep = min_sep
        self.road_sizes = road_sizes

    def compute_priority(self):
        """compute the priority of the constraint given its current state."""
        self.priority = 2
        
    def compute_current_value(self):
        """compute the current value of the constraint."""
        self.current_value = self.agent.get_mean_overlapping_rate(self.min_sep, self.road_sizes)

    def compute_goal_value(self):
        """compute the goal value of the constraint."""
        # do nothing
        return

    def compute_satisfaction(self):
        """compute the satisfaction of the constraint according to the current and goal values."""
        if(self.agent.deleted):
            self.satisfaction = 100.0
            return
        
        self.compute_current_value()
        # the satisfaction is directly the components statisfaction value
        self.satisfaction = 100 - round(self.current_value / 0.005)
        if self.satisfaction < 0:
            self.satisfaction = 0.0

    def compute_actions(self):
        self.actions.clear()
        # proposes displacement and then elimination
        action = RandomBlockDisplacementAction(self,self.agent,3,self.road_sizes, self.min_sep)
        self.actions.append([action, self, 3])
        action2 = PromBlockEliminationAction(self,self.agent,1, 1)
        self.actions.append([action2, self, 1])