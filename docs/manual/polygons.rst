Polygons
========

Buildings
~~~~~~~~~

CartAGen also contains algorithms that process any type of polygons,
and others specific to some types of map polygons, such as buildings.
This section contains only informations related to algorithms that process
one polygon at a time, including:

- :func:`Building simplification <cartagen.simplify_building>` :footcite:p:`ruas:1999`
- :func:`Least square squaring  <cartagen.square_polygon_ls>` :footcite:p:`touya:2016`
- :func:`Rectangle transformation <cartagen.rectangle_transformation>`
- :func:`Recursive regression <cartagen.recursive_regression>` :footcite:p:`bayer:2009` :footcite:p:`yang:2024`

.. plot:: code/manual/buildings_simplification.py

    Building generalisation algorithms

As you can see on the image above, those different algorithms can be used to
address different issues. For example, the building simplification is used to
reduce the complexity of the buildings, the least square squaring algorithm
is used to correct angles by forcing them to go flat or right. The rectangle
transformation is mostly used inside more complex algorithms, such as AGENT
as a possible branch inside a decision tree. The recursive regression can be quite destructive
and is mainly used to correct buildings automatically created from aerial or
satellite imagery.

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
:func:`morphological amalgamation <cartagen.morphological_amalgamation>` :footcite:p:`damen:2008`
and recreate the building blocks.

.. plot:: code/manual/buildings_amalgamation.py

    Morphological amalgamation of buildings into building blocks

Urban areas
~~~~~~~~~~~

Sometimes, you will need to generate a representation of an urban area.
This can be achieved by representing the extent of the buildings if you
consider your buildings as representative of the urban area.
This can be done using the algorithm to calculate
:func:`Boffet area. <cartagen.boffet_areas>` :footcite:p:`boffet:2003`

.. plot:: code/manual/buildings_boffet.py

    Morphological amalgamation of buildings into building
    blocks with a small gaussian smoothing applied to the result