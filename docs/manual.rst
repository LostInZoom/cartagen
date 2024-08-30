.. _manual:

===========
User Manual
===========

:Authors: Guillaume Touya, Justin Berli, Azelle Courtial
:Copyright:
  This work is licensed under a `European Union Public Licence v1.2`__.

.. __: https://eupl.eu/

.. _intro:

Introduction
============

The user manual aims at providing an **overview** of the cartographic generalisation algorithms proposed by Cartagen.
Some complex algorithms and processes are more deeply explained here than in the API Reference, in particular for
those having several steps.

This section is also designed to help you decide **which algorithm best fit your needs**. Indeed, several generalisation
methods are developped to address the same issue, for example to avoid polygon cluttering or to simplify lines and polygons.
Those issues can also depend on the type of geographic objects you are trying to generalise. For example, maybe you
don't want to simplify rivers and roads the same way, as their representation on the map doesn't have the same constraints.

With that in mind, we tried here to offer here an overview of the algorithms offered by Cartagen in order to help you make
a choice in your cartographic endeavors. Please keep in mind that **this is not a lecture on cartographic generalisation**
nor a full explanation of the generalisation process. We have designed several **Jupyter notebooks** to serve as learning
material, they are accessible `here. <https://github.com/LostInZoom/cartagen-notebooks>`_

.. include:: manual/points.rst
  
.. include:: manual/lines.rst

Polygons
========

This library also contains algorithms that process any type of polygons,
and others specific to some types of map polygons, such as buildings.
This section contains only informations related to algorithms that process
one polygon at a time, including:

- :func:`Simplify building <cartagen.simplify_building>`
- :func:`Square polygon  <cartagen.square_polygon_ls>`

.. plot:: code/manual/simplification_buildings.py

  Building simplification algorithms

Groups of objects
=================

.. method:: morphological_amalgamation(buildings, buffer_size, edge_length)

    Amalgamates a group of building polygons using morphological operators, with the algorithm presented in `(Damen et al., 2008) <https://www.semanticscholar.org/paper/High-Quality-Building-Generalisation-by-Extending-Damen-Kreveld/b64618584b3ae3725da7eeb5a545d1580e5f2113>`_. 
    The algorithm chains morphological opening and closing to amalgamate close buildings into a larger building area.
    The 'buffer_size' is parameter used for the opening and closing operations. The 'edge_length' gives the length of edges that are later simplified in the amalgamated polygons.

.. code-block:: pycon

  >>> buildings = [Polygon([(1, 0), (9, 0), (9, 6), (1, 6), (1, 0)]),Polygon([(10, 0), (17, 0), (17, 6), (10, 6), (10, 0)])]
  >>> morphological_amalgamation(buildings, 1.0, 1.0)
  <POLYGON ((1.207 1.983, 2.547 5.885, 16.768 4.282, 15.42 0.148, 1.207 1.983))>

.. plot:: code/building_amalgamation.py

Figure 6. Buildings amalgamated using the algorithm from Damen et al. (2008).

.. method:: reduce_points_kmeans(points, shrink_ratio, centroid_option = False)

    This algorithm reduces a set of points to a smaller set of points that is representative of the initial set. The algorithm uses a K-Means clustering to reduce the set to a number of clusters that corresponds to the shrinking ratio parameter.
    The 'shrink_ratio' parameter can vary between 0 (all points are removed) and 1 (all points are kept).
    Two options are possible: either keeping one of the initial points to replace a cluster (default option) or replace the cluster by its centroid.

.. code-block:: pycon

  >>> points = [Point(1,1), Point(1,2), Point(0,1), Point(2,1), Point(2,2), Point(5,5), Point(8,10), Point(10,10), Point(10,8), 
              Point(16,10), Point(16,9), Point(14,11)]
  >>> reduce_points_kmeans(points, 0.25)
  [<POINT (2.0 2.0)>, <POINT (10.0 10.0)>, <POINT (16.0 10.0)>]

.. plot:: code/kmeans_reduction.py

Figure 8. A set of points reduced to 25% of its initial amount, with the K-Means reduction algorithm.

.. method:: reduce_points_quadtree(points, depth, mode='simplification', attribute = "")

    Algorithm to reduce a point set based on a quadtree. The algorithm was proposed by Bereuter & Weibel (2012). The algorithm uses a quadtree that divdes itself until there is only one point feature inside the cell.
    The 'depth' parameter can vary between 0 (all points are removed) and the maximum depth of the quadtree (all points are kept). If depth is 2, the algorithm only retains 1 point for each cell with depth <= 2.
    Three options are possible to choose how the point is retained in the cell: 
    - mode = 'selection' means that for one cell, the algorithm retains the point with the largest value in the chosen attribute, weighted by the depth of the point. 
    - mode = 'simplification' means that the point retained in the cell is the closest to the center of the cell
    - mode = 'aggregation' means that the points are all aggregated to the centroid of the points.
    The algorithm returns the list of tuples (geometry, index, nb_of_points) where index is the index of the point in the initial Geodataframe (-1 if the point was created), and nb_of_points gives the amount of initial points replaced (which can be used to weight the size of the symbol of this point). 

.. code-block:: pycon

  >>> points = [Point(1,1), Point(1,2), Point(0,1), Point(2,1), Point(2,2), Point(5,5), Point(8,10), Point(10,10), Point(10,8), 
              Point(16,10), Point(16,9), Point(14,11)]
  >>> p1 = gpd.GeoSeries(points)
  >>> gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(p1))
  >>> reduce_points_quadtree(points, 0.25)
  [<POINT (1 2)>, <POINT (10 8)>, <POINT (10 10)>, <POINT (14 11)>]

.. plot:: code/quadtree_reduction.py

Figure 9. A set of points reduced to depth 2 of the quadtree, with the selection mode. The selected points are displayed in red.

Enrich your data prior to map generalisation
--------------------------------------------

Since the beginning of research on the automation of map generalisation, the necessity for enrichment has been clear. There are properties, structures, which are implicit in the spatial arrangement of geometries in the map. These properties, structures are necessary to make the best decision when generalising the map, and this data enrichment step helps by making these properties, these structures explicit cartographic data.

Extracting implicit geographic structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: boffet_areas(buildings, dilation_size, erosion_size, simplification_distance = 2)

    Computes the urban/built-up areas from a set of buildings, using a method from Boffet (2000). The algorithm computes buffers around each building ('dilation_size') and then merges all buffers.
    The merged areas are then further refined with an erosion ('erosion_size') and a Douglas & Peucker simplification ('simplification_distance').

.. code-block:: pycon

  # compute the built-up areas with a 25 m buffer and a 10 m erosion
  urbanareas = boffet_areas(polygons, 25.0, 10.0)

.. plot:: code/boffet_areas.py

Figure 10. Building polygons converted into built-up areas using the Boffet algorithm.


Measures on map features
^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: polygon_compactness(polygon)

    Returns the compactness of a ''Polygon'' using the Miller index, i.e. 4.Pi.area / perimeter². This index gives a maximum value of 1 for circles.

.. code-block:: pycon

  >>> polygon = Polygon([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
  >>> polygon_compactness(polygon)
  0.7853981633974483

.. method:: polygon_elongation(polygon)

    Returns the elongation of a ''Polygon'' using the measure from the AGENT project, i.e. the longest edge of the minimum bounding rectangle (MBR) divided by the shortest edge of the MBR.

.. code-block:: pycon

  >>> polygon = Polygon([(0, 0), (0, 10), (5, 10), (5, 0), (0, 0)])
  >>> polygon_elongation(polygon)
  2.0

.. method:: building_min_width(building)

    Returns the minimum width inside a building. The minimum width is the minimum distance between two edges of the buildings that are not adjacent. 
    The measure was proposed during the AGENT project. 'building' should be a shapely ''Polygon''.

.. code-block:: pycon

  >>> polygon = Polygon([(0, 0), (0, 10), (2, 10), (2, 6), (5, 6), (5, 0), (0, 0)])
  >>> building_min_width(polygon)
  2.0

Apply map generalisation complex processes
------------------------------------------

AGENT model
^^^^^^^^^^^
This user guide is not meant to fully explain the principles of the AGENT model, and how it works. If you are not familiar with the AGENT model, please read the scientific papers describing this model:
- `<http://icaci.org/files/documents/ICC_proceedings/ICC2001/icc2001/file/f13041.pdf>`_
- `<http://dx.doi.org/10.1016/b978-008045374-3/50016-8>`_
- `<https://hal.inria.fr/IFSTTAR/hal-01682131v1>`_
Though it is a tutorial for the JAVA CartAGen platform, this `webpage <https://ignf.github.io/CartAGen/docs/tuto_agents.html>`_ also contains complementary information on how to trigger agent-based map generalisation.

Micro agents
=============
You may want to use micro agents only, i.e. one cartographic feature such as a
building that generalises itself without consideration for the other cartographic features.
To generalise micro agents, you have to follow these steps:

  1. Create agent objects from your Geopandas features

  .. code-block:: pycon

    geoms = [loads('Polygon Z ((395038.7 6272970.9, 395030.4 6272984, 395025.3 6272982, 395023.2 6272983.7, 395020 6272981.3, 395016.9 6272985.9, 395021.8 6272990.7, 395020.6 6272993.7, 395024.7 6272997.2, 395028.5 6272994.5, 395032.8 6272988.2, 395038.1 6272991.6, 395044.9 6272979.1, 395047.1 6272980.4, 395049.5 6272976.8, 395038.7 6272970.9))'),
      loads('Polygon Z ((394999.5 6272975, 395006.7 6272962.4, 395010.6 6272957.5, 394996.6 6272944.4, 394991 6272949, 394999.2 6272956.3, 394996.1 6272959.7, 394998.3 6272961.3, 394992 6272969.4, 394999.5 6272975))'),
      loads('Polygon Z ((395007.3 6272975.8, 395013.2 6272981, 395021.2 6272969.6, 395024.2 6272971.9, 395031 6272963.8, 395020.8 6272957.4, 395007.3 6272975.8))'),
      loads('Polygon Z ((395082.3 6272967.4, 395089.9 6272958, 395071.9 6272945.9, 395068.4 6272950.6, 395066 6272949, 395056.3 6272962, 395058.5 6272963.5, 395056.40000000002328306 6272966.8, 395059.4 6272969.9, 395056.9 6272972.6, 395054.5 6272968.3, 395049.6 6272973.4, 395058.4 6272981.6, 395073.6 6272962.5, 395082.3 6272967.4))')
      ]
    envgdf = geopandas.GeoDataFrame(geometry=geopandas.GeoSeries(geoms))
    for index, feature in envgdf.iterrows():
      agent = BuildingAgent(feature)
  
  2. Add constraints to your agents. You can pick among the constraints provided by the library, but you can also code new constraints and add them to your agents. The list of default constraints is provided in the table below.

  .. list-table:: Building micro constraints
   :widths: 50 20 50
   :header-rows: 1

   * - name
     - property
     - actions
   * - BuildingSizeConstraint
     - area
     - enlarge, delete
   * - BuildingGranularityConstraint
     - granularity
     - simplify, simplify to rectangle
   * - BuildingSquarenessConstraint
     - squareness
     - squaring

  .. code-block:: pycon

    squareness = BuildingSquarenessConstraint(1,agent)
    size = BuildingSizeConstraint(1, agent, 250)
    granularity = BuildingGranularityConstraint(1, agent, 6)
    agent.constraints.append(size)
    agent.constraints.append(squareness)
    agent.constraints.append(granularity)

  3. Run the agents, i.e. start their life cycle iteratively. To run the agents, you have to use: 
  .. method:: run_agents(agents, lifecycle='basic', store_states=False, verbose=0)

    - ''agents'' is a list of agents to run.
    - ''lifecycle'' chooses the type of life cycle to apply on the agents (only the basic is implemented for now)
    - ''store_states'' is True if you want to get all the intermediate states of the agents as output of the function.
    - ''verbose'' defines how much detail is logged during the agents life cycle.
  
  .. code-block:: pycon

    run_agents(agents_to_run)

Meso agents
=============
The implementation of the meso agents is not yet completed. 

Bibliography
============

.. footbibliography::