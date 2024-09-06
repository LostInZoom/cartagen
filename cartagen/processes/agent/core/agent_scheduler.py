def run_agents(agents, lifecycle='basic', store_states=False, verbose=0):
    """
    Execute the generalisation process on the given agents.

    This function executes the AGENT process on the given list of agents.

    Parameters
    ----------
    agents : list of Agent
        The agents to use for generalisation.
    lifecycle : str, optional
        Type of life cycle to apply on the agents.
    store_states : bool, optional
        If True, the function will output all intermediate states
        of the agents.
    verbose : int, optional
        Verbose level.
    """
    accepted_lifecycles = [
        'basic'
    ]

    if lifecycle not in accepted_lifecycles:
        raise Exception('Lifecycle type not handled: {0}'.format(lifecycle))

    while (len(agents) != 0):
        agent = agents.pop()
        
        if(verbose > 0):
            print("agent {} is processed by the scheduler.".format(agent.id))

        if(lifecycle == 'basic'):
            # run the basic lifecycle on the current agent
            __activate_agent_basic(agent, verbose=verbose)

def __activate_agent_basic(agent, store_states=False, validity_satisfaction = 0.5, verbose=0, iter_max = 20):
    # clean the possibly previously created and stored states
    if(store_states):
        agent.clean_states()

    # compute the satisfaction of the agent
    agent.compute_satisfaction()

    if verbose > 0:
        print("satisfaction value = {}".format(agent.satisfaction))

    # get actions from constraints
    agent.update_action_proposals()

    if verbose > 2:
        print("proposed actions:")
        for action in agent.actions_to_try:
            print(action)

    # if we store states, we create the current state as the root state
    if(store_states):
        # TODO
        current_state = None
    
    # test if the agent is satisfied
    if(agent.satisfaction >= 100.0 - validity_satisfaction):
        return 
    
    # the agent is not satisfied: try to improve its satisfaction by triggering some actions
    i = 0
    while(i < iter_max):
        # test if there are actions to try
        if(len(agent.actions_to_try)==0):
            # all possible action have been tried: the best possible state as been reached
            break
        
        # take the next action to try from the list
        action_proposal = agent.get_best_action_proposal()
        
        if(action_proposal == None):
            if verbose > 0:
                print("agent {} has no more action to try".format(agent.id)) 
            return
        action = action_proposal[0]
        #agent.remove_action_to_try(action_proposal)

        if verbose > 0:
            print("selected action: {}".format(action_proposal))

        # store the geometry of the agent before the action
        previous_geom = agent.feature['geometry']
        # trigger the action
        action.compute()

        # compute the new satisfaction
        previous_satisfaction = agent.satisfaction
        agent.compute_satisfaction()

        if verbose > 0:
            print("new satisfaction after the action: {}".format(agent.satisfaction))

        # get actions from constraints
        # agent.update_action_proposals()
        # remove the action that we just tested to avoid trying the same action if there is a backtrack
        # agent.remove_action(action_proposal)

        if verbose > 2:
            print("new proposed actions:")
            for action in agent.actions_to_try:
                print(action)

        # test if the agent is satisfied
        if(agent.satisfaction >= 100.0 - validity_satisfaction):
            if verbose > 0:
                print("agent {} is satisfied".format(agent.id))
            return 
        
        # check if the new state is valid
        if (agent.satisfaction - previous_satisfaction > validity_satisfaction):
            # the agent has improved his state, but other actions are necessary
            if verbose > 0:
                print("agent {} improved but is still unsatisfied".format(agent.id)) 
            i += 1
            continue
        else:
            # the state is not valid and the agent is backtracked to its previous state
            if verbose > 0:
                print("agent {} did not improve so it is backtracked to previous state".format(agent.id)) 
            agent.feature['geometry'] = previous_geom
            i += 1



