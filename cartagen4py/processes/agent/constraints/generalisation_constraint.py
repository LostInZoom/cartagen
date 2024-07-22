
class GeneralisationConstraint:
    """
    Abstract generalisation constraint object for the AGENT process.
    """

    actions = []
    current_value = None
    goal_value = None
    priority = 0
    satisfaction = 100.0

    def __init__(self, importance, agent):
        self.__importance = importance
        self.__agent = agent

    # def compute_priority(self):
    #     """compute the priority of the constraint given its current state."""
    #     pass
        
    # def compute_current_value(self):
    #     """compute the current value of the constraint."""
    #     pass

    # def compute_goal_value(self):
    #     """compute the goal value of the constraint."""
    #     pass

    # def compute_satisfaction(self):
    #     """compute the satisfaction of the constraint according to the current and goal values."""
    #     pass

    # def compute_actions(self):
    #     """compute the actions of the constraint according to the current and goal values."""
    #     pass