.. _processes:

#########
Processes
#########

*****
AGENT
*****

Agent-based generalization follows a step by step, knowledge based approach,
which arises from the seminal model from Ruas & Plazanet :footcite:p:`ruas:1996`.
The aim is not to find a global solution to generalize the whole map at once,
but to identify local situations that require generalization,
and try generalization algorithms to improve the legibility of this situation.

This process is a port of the version present in the
`CartAGen Java application <https://ignf.github.io/CartAGen/>`_.
For a better understanding of the AGENT process, you can refer to
Barrault *et al.* :footcite:p:`barrault:2001`, Ruas & Duchêne :footcite:p:`ruas:2007`
and Duchêne *et al.* :footcite:p:`duchene:2018`


.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    run_agents

Micro-agents
------------

Micro-agents are the simpler kind of agents as they aim at generalising single
cartographic features without consideration for their surroundings.

Buildings
^^^^^^^^^

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    BuildingAgent
    BuildingSizeConstraint
    BuildingGranularityConstraint
    BuildingSquarenessConstraint

Meso-agents
-----------

Building blocks
^^^^^^^^^^^^^^^

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    BlockAgent

.. footbibliography::