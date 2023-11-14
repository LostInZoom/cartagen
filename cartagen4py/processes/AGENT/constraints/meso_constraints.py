from cartagen4py.processes.AGENT.constraints.generalisation_constraint import GeneralisationConstraint
from cartagen4py.processes.AGENT.actions.generalisation_action import MesoMicroActivation

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
        # propose the action that activates the micro agents
        action = MesoMicroActivation(self,self.agent,1)
        self.actions.append([action, self, 1])