Networks
========

CartAGen also proposes network-related algorithms intended for the generalisation
of different type of artificial or natural linear networks. 

Roads
~~~~~

Road networks algorithms are designed to detect or collapse specific road network features
such as roundabouts, branching crossroads or dual carriageways. They can be dependent on each
other, for example, collapsing roundabouts depends on both roundabout detection and branching
crossroads detection to yield better results.

Having that in mind, the following section illustrates an example of a network generalisation
workflow using network enrichment functions and collapsing algorithms. Let's consider the following
road network as the base data, which has too much detail for the scale we want.

.. plot:: code/manual/network_original.py

    The original road network.

Roundabouts and branching crossroads
------------------------------------

We can start by using :func:`detect_roundabouts <cartagen.detect_roundabouts>` followed by
:func:`detect_branching_crossroads <cartagen.detect_branching_crossroads>` to respectively
detect roundabouts and then branching crossroads. We can specify detected roundabouts as
an input when detecting branching crossroads, that way, when those two objects are connected,
they can be treated as one object when collapsing.

.. code-block:: Python

    # The variable network is a GeoDataFrame of LineString representing
    # the network to generalise
    roundabouts = detect_roundabouts(network)
    crossroads = detect_branching_crossroads(network, roundabouts)

Then, we can apply :func:`collapse_branching_crossroads <cartagen.collapse_branching_crossroads>`
followed by :func:`collapse_roundabouts <cartagen.collapse_roundabouts>` to collapse those
shapes and simplify the network. The order is important as the first will collapse branching
crossroads that are not connected to roundabouts, and then collapse roundabouts and mixed
objects (roundabouts with incoming branching crossroads).

.. code-block:: Python

    network = collapse_branching_crossroads(network, crossroads)
    network = collapse_roundabouts(network, roundabouts, crossroads)

.. plot:: code/manual/network_phase1.py

    Collapsing roundabouts (in red) and branching crossroads (in blue) inside the road network.

Dual carriageways
-----------------

We can now detect dual carriageways using
:func:`detect_dual_carriageways <cartagen.detect_dual_carriageways>`...

.. code-block:: Python

    carriageways = detect_dual_carriageways(network)

...and collapse them using
:func:`collapse_dual_carriageways <cartagen.collapse_dual_carriageways>`.

.. code-block:: Python

    network = collapse_dual_carriageways(network, carriageways)

.. plot:: code/manual/network_phase2.py

    Collapsing dual carriageways inside the road network.

Dead ends
---------

Now, we can detect dead-ends by using :func:`detect_dead_ends <cartagen.detect_dead_ends>`.
This stage can be deployed at the discretion of the cartographer as it can be used to
clean up the network. Here, we used this after collapsing branching crossroads as the 
dead-end located at the center of the current network would not have been detected.
You can read more about how dead-ends are detected in their API reference page.

.. code-block:: Python

    network = detect_dead_ends(network, True)

Finally, we can eliminate dead-end groups using :func:`eliminate_dead_ends <cartagen.eliminate_dead_ends>`
that are shorter than 200 meters and simplify long dead-end groups by keeping only the longest path.

.. code-block:: Python

    network = eliminate_dead_ends(network, 200, keep_longest=True)

.. plot:: code/manual/network_phase3.py

    Eliminate or simplify dead-end groups (in red) inside the network.

We can then compare our input network and our result to see how the network generalisation
simplified the roads layout details in order to represent them at a larger scale.
Keep in mind that generalisation of a road network is a complex process and sometimes
algorithms can have a hard time finding specific features. For example, you could have
a branching crossroad that is not triangular, and thus, it won't be detected as one.
This is a reminder that good generalisation relies on the cartographer's validation
at the end of the process.

.. plot:: code/manual/network_comparison.py

    The original network (left) and the network after generalisation (right).