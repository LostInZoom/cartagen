
class GeneralisationAction():
    weight = 0
    constraint = None
    agent = None
    restriction = -1

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
    
    def compute(self):
        """Compute the action, i.e. triggers the algorithm."""
        self.agent.deleted = True
    