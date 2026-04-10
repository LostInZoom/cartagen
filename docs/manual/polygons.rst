Polygons
========

Buildings
~~~~~~~~~

CartAGen also contains algorithms that process any type of polygons,
and others specific to some types of map polygons, such as buildings.
This section contains only informations related to algorithms that process
one polygon at a time. If your buildings seems to be adapted to your current
scale, you can still run a squaring algorithm to correct angles that are almost
flat, right, or at 45 degrees, using:

- :func:`Least square squaring  <cartagen.square_polygon_ls>` :footcite:p:`touya:2016`
- :func:`Naive squaring  <cartagen.square_polygon_naive>` :footcite:p:`touya:2016`

.. plot:: code/manual/buildings_squaring.py

    Building squaring algorithms

If you want to modify further the geometry of the building, you can use simplification algorithms,
such as:

- :func:`Building simplification <cartagen.simplify_building>` :footcite:p:`ruas:1999`

.. plot:: code/manual/buildings_simplification.py

    Building simplification

The building simplification algorithm is a simple algorithm design to reduce the complexity
of buildings by removing edges.

CartAGen contains algorithms designed to process buildings computed from segmentation of
aerial or satellite images. Those algorithms can drastically 

- :func:`Rectangle transformation <cartagen.regularize_building_rectangle>`
- :func:`Recursive regression <cartagen.regularize_building_regression>` :footcite:p:`bayer:2009` :footcite:p:`yang:2024`
- :func:`Feature Edge Reconstruction (FER) <cartagen.regularize_building_fer>` :footcite:p:`yang:2024`

.. plot:: code/manual/buildings_regularization.py

    Building regularization algorithms

Recursive regression can be quite destructive and is will always return 45 degrees angles.
On the other hand, the rectangle transformation is mostly used inside more complex algorithms, such as AGENT
as a possible branch inside a decision tree. Feature edge reconstruction is a complex algorithm
with a decision tree that is design to keep the overall shape of the building.

CartAGen also contains algorithms that process several buildings at a time,
usually to create amalgamated representation of buildings.

- :func:`Building amalgamation <cartagen.amalgamate_buildings>` :footcite:p:`regnauld_generalisation_1998`

.. plot:: code/reference/amalgamate_buildings.py

    Building amalgamation algorithm

Building blocks
~~~~~~~~~~~~~~~

A building block can be defined as several buildings inside the face of a road
network, and usually in an urban area to have a high enough density to require generalisation.
Different operations can then be applied depending on factors such as the density of buildings and the
final scale of the map.

At large scale, maybe we need to keep the individual buildings without modifying their geometry and the
space inside the building block is enough for the buildings to be moved. Thus,
we can use the :func:`random displacement <cartagen.random_displacement>` algorithm.

.. plot:: code/manual/buildings_displacement.py

    Random displacement of buildings inside the building block

Still at large scale, sometimes there is not enough space inside the building block to move
the buildings around. It often is the case in urban areas, where buildings are touching each other.
So, you may want to aggregate those buildings into a representation of the building blocks that takes
into consideration the shapes of the individual buildings. For that purpose, you can use the
:func:`morphological amalgamation <cartagen.amalgamate_buildings_morphological>` :footcite:p:`damen:2008`
and recreate the building blocks.

.. plot:: code/manual/buildings_amalgamation.py

    Morphological amalgamation of buildings into building blocks

In dense suburban areas, you might want to use :func:`building typification <cartagen.typify_buildings_matching>` :footcite:p:`li:2005`
to enhance the legibility of the map.

.. plot:: code/manual/buildings_typification.py

    Building typification in dense suburban areas

Urban areas
~~~~~~~~~~~

Sometimes, you will need to generate a representation of an urban area.
This can be achieved by representing the extent of the buildings if you
consider your buildings as representative of the urban area.
This can be done using the algorithm to calculate
:func:`Boffet area. <cartagen.boffet_areas>` :footcite:p:`boffet:2003`

.. plot:: code/manual/buildings_boffet.py

    Boffet areas with a small gaussian smoothing applied to the result