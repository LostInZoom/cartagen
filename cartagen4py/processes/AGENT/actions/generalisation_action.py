from cartagen4py.processes.AGENT.core.agent_scheduler import *

class GeneralisationAction():
    weight = 0
    constraint = None
    agent = None
    restriction = -1
    name = ""

    def __init__(self, constraint, agent, weight):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint

    def clean(self):
        self.agent = None
        self.constraint = None
    
    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        pass

class DeletionAction(GeneralisationAction):
    
    def __init__(self, constraint, agent, weight):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.name="Deletion"
    
    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        self.agent.deleted = True

class MesoMicroActivation(GeneralisationAction):

    def __init__(self, constraint, agent, weight):
        self.weight = weight
        self.agent = agent
        self.constraint = constraint
        self.name="MicroActivation"
    
    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        components_to_trigger = []
        for component in self.components:
            if component.deleted == False:
                components_to_trigger.append(component)
        
        nb = len(components_to_trigger)
        i = 1
        while(len(components_to_trigger)>0):
            best = self.agent.get_best_component_to_activate(components_to_trigger)
            run_agents([best])
            self.agent.manage_internal_side_effects(best)
            components_to_trigger.remove(best)
    