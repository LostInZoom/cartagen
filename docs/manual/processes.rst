Processes
=========

AGENT
~~~~~

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

Micro-agents
------------

Micro-agents are the simpler kind of agents as they aim at generalising single
cartographic features without consideration for their surroundings.

Assuming the variable :code:`buildings` is a GeoDataFrame of Polygon, we
need to loop through buildings and add them as individual agent before
generalisation:

.. code-block:: Python

    import cartagen as c4
	
    agents = []
    for i, building in buildings.iterrows():
        # Create the agent
        agent = c4.BuildingAgent(building)

        # Adding a size constraint of 250 square meters to enlarge the building if needed
        size = c4.BuildingSizeConstraint(agent, importance=1, min_area=250)

        # Adding a granularity constraint that will simplify the building if
        # an edge is above the minimum allowed length
        granularity = c4.BuildingGranularityConstraint(agent, importance=1, min_length=6)

        # Adding a squareness constraint to square the building angles if needed
        squareness = c4.BuildingSquarenessConstraint(agent, importance=1)
        
        agent.constraints.append(size)
        agent.constraints.append(squareness)
        agent.constraints.append(granularity)
        agents.append(agent)
    
    c4.run_agents(agents)

Below is a table representing the different constraints available for the building generalisation
using AGENT along with their associated actions.

.. list-table:: Building micro constraints
    :widths: 50 20 50
    :header-rows: 1

    * - Name
      - Property
      - Actions
    * - BuildingSizeConstraint
      - area
      - enlarge, delete
    * - BuildingGranularityConstraint
      - granularity
      - simplify, simplify to rectangle
    * - BuildingSquarenessConstraint
      - squareness
      - squaring

.. plot:: code/manual/process_agent_micro.py

	Resulting generalisation using the micro-agents

Least squares-based method
~~~~~~~~~~~~~~~~~~~~~~~~~~

CartAGen also propose a generalisation method based on the method of least squares. This
method was proposed by Harrie :footcite:p:`harrie:1999` and relies on the least squares regression
to displace every vertexes of every provided objects (points, lines and polygons). Those vertexes
can move but they need to respect user-provided constraints.
In practice, it iteratively tries to resolve matrice equations until a threshold norm
is reached keeping every constraints as satisfied as possible.

Setup
-----

As en example, let's start with a set of three geographic objects (example in the image below), namely
some buildings (in gray), a simple road network (in red) and a river (in blue).

Let's consider those objects are too close to each other at the scale we want
to represent them. We can start by instanciating the least squares method class with default parameters:

.. code-block:: Python

    ls = cartagen.LeastSquaresMethod()

Constraints
-----------

But we also want to specify constraints to maintain some control over the movement of vertexes.
Using the following table, we can specify different constraints for the different objects:

.. list-table:: Constraints description used by the method of least squares
    :widths: 20 30 50
    :header-rows: 1

    * - Constraint
      - Geometry type
      - Impact
    * - movement
      - Point, LineString, Polygon
      - The object should move as little as possible
    * - stiffness
      - LineString, Polygon
      - The internal geometry should be invariant, `i.e.` the vertexes movement within the same object
        will try not to move closer or away from each other.
    * - curvature
      - LineString, Polygon
      - The curvature of a line or a polygon border should not change, `i.e.` the angle formed by
        two connected segments will try not to change.

Having that in mind, we can add the different objects with the constraints properly set. The numbers represent the weight of the
constraint. A higher value means the constraint will be more respected during generalisation. So, in the following
example, the buildings have a high stiffness which means there shape will tend to stay the same:

.. code-block:: Python

    ls.add(buildings, movement=2, stiffness=10)
    ls.add(roads, movement=2, curvature=5)
    ls.add(rivers, movement=2, curvature=5)

Spatial relationships
---------------------

We decided to let users more freedom when setting the spatial conflict constraints. Thus, we can decide
the distance and weight of spatial conflicts between each pair of objects we already added.
So, we need to provide two matrices: one for the distances between pairs of objects, one for the
weight of those spatial constraints.

First, we retrieve the number of objects added for generalisation. This should be 3: buildings, roads, rivers:

.. code-block:: Python

    d = ls.get_objects_number()

We can then create two numpy arrays of the right shape:

.. code-block:: Python

    distances = np.zeros((d, d))
    spatial_weights = np.zeros((d, d))

We then can populate the new arrays with default distance and weight values:

.. code-block:: Python

    for i in range(d):
        for j in range(d):
            distances[i][j] = 20
            spatial_weights[i][j] = 5

Those two arrays represent the spatial conflicts and the weight of those conflicts between each pairs of objects. So, by setting a default distance,
we will force objects to have at least 20 meters between each other.

But maybe we want to allow buildings to be closer to each other and further away from the rivers.
To change those distances, we can use :code:`ls.distances[n][m]` where :code:`n` and :code:`m` are
respectively the index of the pair of object we want to modify in the order we added them.
Let's say we want buildings to be at least at 25 meters from each other, 28 meters from the roads and 30 meters
from the rivers, we can do:

.. code-block:: Python

    distances[0][0] = 25
    distances[0][1] = 28
    distances[0][2] = 30

The :code:`distances` variable now return:

.. code-block:: Python

    [[25. 28. 30.]
     [20. 20. 20.]
     [20. 20. 20.]]

And if we want to increase the weight of the spatial conflict between the roads and the rivers, we can change it
this way:

.. code-block:: Python

    spatial_weights[1][2] = 8

The :code:`spatial_weights` variable now return:

.. code-block:: Python

    [[5. 5. 5.]
     [5. 5. 8.]
     [5. 5. 5.]]

Finally, we can add the spatial relationships to the least squares method:

.. code-block:: Python

    ls.add_spatial_conflicts(distances, spatial_weights)

Results
-------

Finally, we can launch the generalisation that will return a tuple with as many as objects we added:

.. code-block:: Python

    buildings, roads, rivers = ls.generalise()

.. plot:: code/manual/process_ls_original.py
    
    Geographic objects before generalisation

.. plot:: code/manual/process_ls_conflicts.py
    
    Spatial conflicts detected

.. plot:: code/manual/process_ls_results.py
    
    Results of the generalisation