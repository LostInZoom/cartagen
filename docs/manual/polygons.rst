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
transformation is mostly used inside more complex algorithms, such as AGENT :ref:`AGENT`
as a possible branch inside a decision tree. The recursive regression can be quite destructive
and is mainly used to correct buildings automatically created from aerial or
satellite imagery.

Building blocks
~~~~~~~~~~~~~~~

Building blocks are generally defined as multiple buildings inside the face of a road
network, and usually in an urban area to have a high enough density to require generalisation.
