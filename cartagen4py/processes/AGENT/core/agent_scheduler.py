from cartagen4py.processes.AGENT import *

def run_agents(agents,lifecycle='basic',store_states=False):
    
    while (len(agents) != 0):
        agent = agents.pop()

        if(lifecycle == 'basic'):
            # run the basic lifecycle on the current agent
            __activate_agent_basic(agent)

def __activate_agent_basic(agent, store_states=False, validity_satisfaction = 0.5):
    # clean the possibly previously created and stored states
    if(store_states):
        agent.clean_states()

    # compute the satisfaction of the agent
    agent.compute_satisfaction()

    # get actions from constraints
    agent.update_action_proposals()

    # if we store states, we create the current state as the root state
    if(store_states):
        # TODO
        current_state = None
    
    # test if the agent is satisfied
    if(agent.satisfaction >= 100.0 - validity_satisfaction):
        return 
    
    # the agent is not satisfied: try to improve its satisfaction by triggering some actions
    while(True):
        # test if there are actions to try
        if(len(agent.actions_to_try)==0):
            # all possible action have been tried: the best possible state as been reached
            break

        # take the next action to try from the list
        action_proposal = agent.get_best_action_proposal()
        action = action_proposal[0]

        # store the geometry of the agent before the action
        previous_geom = agent.feature['geometry']
        # trigger the action
        action.compute()

        # compute the new satisfaction
        previous_satisfaction = agent.satisfaction
        agent.compute_satisfaction()
        # get actions from constraints
        agent.update_action_proposals()
        # remove the action that we just tested to avoid trying the same action if there is a backtrack
        agent.remove_action(action_proposal)

        # test if the agent is satisfied
        if(agent.satisfaction >= 100.0 - validity_satisfaction):
            return 
        
        # check if the new state is valid
        if (agent.satisfaction - previous_satisfaction > validity_satisfaction):
            # the agent has improved his state, but other actions are necessary
            continue
        else:
            # the state is not valid and the agent is backtracked to its previous state
            agent.feature['geometry'] = previous_geom



