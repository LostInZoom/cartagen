.. _user-guide:

==========
User Guide
==========

Please note that this user guide is still under construction and some functions of the library are still not enough documented or not documented at all. 
Feel free to contact us if you need documentation for a specific function.

Apply map generalisation operations
-----------------------------------

Operations for lines
^^^^^^^^^^^^^^^^^^^^

.. method:: visvalingam_whyatt(line, area_tolerance)

    Returns a simplified version of the line using the Visvalingam-Whyatt algorithm `(Visvalingam & Whyatt, 1993) <https://www.tandfonline.com/doi/abs/10.1179/000870493786962263?journalCode=ycaj20>`_.
    The 'area_tolerance' is the minimum area of the triangle formed by three consecutive vertices, to keep the middle vertex in the simplified line.

.. code-block:: pycon

  >>> line = LineString([(2, 0), (2, 4), (3, 4), (3, 5), (5, 7)])
  >>> visvalingam_whyatt(line, 1)
  <LINESTRING (2 0, 2 4, 3 5, 5 7)>

.. plot:: code/visvalingam.py

Figure 1. Two polylines simplified with the Visvalingam-Whyatt algorithm.


.. method:: raposo_simplification(line, initial_scale, final_scale, centroid=True, tobler=False)

    Returns a simplified version of the line using the Raposo algorithm `(Raposo, 2013) <http://dx.doi.org/10.1080/15230406.2013.803707>`_.
    The algorithm uses an hexagonal tessallation, with a size related to the final scale, and it only retains one vertex per hexagonal cell.
    Be careful, it uses the scale as parameter. If the 'centroid' parameter is ''True'', the vertices inside an hexagon cell are replaced by the centroid of the cell; if it is ''False'', they are replaced by the nearest vertex to the centroid of the cell.
    The Raposo algorithm is dedicated to the simplification of natural lines such as rivers, lakes or forests.

.. code-block:: pycon

  line = LineString([(2, 0), (2, 4), (3, 4), (3, 5), (5, 7)])
  simplified_line = raposo_simplification(line, 10000.0, 50000.0)

.. plot:: code/raposo.py

Figure 2. Two polylines simplified with the Raposo algorithm.

Operations for polygons
^^^^^^^^^^^^^^^^^^^^^^^

cartagen4py contains algorithms that process any type of polygons, and others specific to some types of map polygons, such as buildings. Only the algorithms that process one polygon at a time are documented in this section.

.. method:: building_simplification_ruas(building, edge_threshold, parallel_limit = 20 * pi / 180, orthogonal_limit = 20 * pi / 180)

    Returns a simplified version of the building polygon using the algorithm from Anne Ruas (1988). The algorithm was later used to simplify buildings in the AGENT project.
    The 'edge_threshold' is the minimum length of an edge between two vertices of the building to be removed. The 'parallel_limit' and 'orthogonal_limit' parameters define to what extent two edges are considered parallel or orthogonal.

.. code-block:: pycon

  >>> building = Polygon([(0, 0), (0, 10), (2, 10), (2, 9), (10, 9), (10, 0), (0, 0)])
  >>> building_simplification_ruas(building, 2.5)
  <POLYGON ((0 0, 0 9.5, 10 9.5, 10 0, 0 0))>

.. plot:: code/building_simplification.py

Figure 3. Four buildings simplified with the Ruas algorithm.

.. method:: square_polygon(polygons, max_iteration=1000, norm_tolerance=0.05, right_tolerance=10, flat_tolerance=10, fixed_weight=5, right_weight=100, flat_weight=50)

    Squares the angles of a polygon using the algorithm from `(Lokhat & Touya, 2016) <https://josis.org/index.php/josis/article/view/72>`_. The algorithm is based on a least squares adjustment, where angles that are almost right or almost flat are adjusted to be exactly right, or exactly flat.

.. code-block:: pycon

  >>> building = Polygon([(0, 0), (0, 10), (9.8, 9.8), (10, 0), (0, 0)])
  >>> square_polygon([building])
  [POLYGON((-0.00002213 -0.00002213, 0.00000159  9.90002291, 9.89999763  9.89999763, 9.90002291  0.00000159, -0.00002213 -0.00002213))]

.. plot:: code/building_squaring.py

Figure 4. Two buildings squared with the algorithm from (Lokhat & Touya, 2016).

Operations for groups of objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: morphological_amalgamation(buildings, buffer_size, edge_length)

    Amalgamates a group of building polygons using morphological operators, with the algorithm presented in `(Damen et al., 2008) <https://www.semanticscholar.org/paper/High-Quality-Building-Generalization-by-Extending-Damen-Kreveld/b64618584b3ae3725da7eeb5a545d1580e5f2113>`_. 
    The algorithm chains morphological opening and closing to amalgamate close buildings into a larger building area.
    The 'buffer_size' is parameter used for the opening and closing operations. The 'edge_length' gives the length of edges that are later simplified in the amalgamated polygons.

.. code-block:: pycon

  >>> buildings = [Polygon([(1, 0), (9, 0), (9, 6), (1, 6), (1, 0)]),Polygon([(10, 0), (17, 0), (17, 6), (10, 6), (10, 0)])]
  >>> morphological_amalgamation(buildings, 1.0, 1.0)
  <POLYGON ((1.207 1.983, 2.547 5.885, 16.768 4.282, 15.42 0.148, 1.207 1.983))>

.. plot:: code/building_amalgamation.py

Figure 5. Buildings amalgamated using the algorithm from Damen et al. (2008).


.. class:: BuildingDisplacementRandom(max_trials=25, max_displacement=10, network_partitioning=False, verbose=False)

    This algorithm displaces buildings that overlap with each other and/or other features. The algorithm was never published but was available in CartAGen. It is an iterative process that selects the building with the most overlaps, and then pushes slightly the building in a random direction. If the overlaps are reduced, the displacement is commited and a new iteration starts. If the overlaps are worse, the displacement is backtracked, and another one is tried.
    The 'max_trials' parameter gives the maximum number of random displacements tried on one building. The 'max_displacement' parameter is the maximum distance a building is allowed to move. For large datasets, the algorithm can work on smaller partitions, using the 'network_partitioning' parameter.
    The name of the class mentions buildings but other objects can be similarly displaced, as long as GeoDataframe with polygons is provided.

.. method:: displace(self, buildings, roads, rivers, *networks)

    This method displaces the buildings with roads and rivers acting as obstacles for the buildings, i.e. the algorithm minimises the overlaps between buildings and with the geometries contained in those two collections.
    'buildings', 'roads', and 'rivers' are geopandas GeoDataframes, not arrays of geometries. If you want to avoid overlaps with road and river symbols, you need to provide polygons as the main geometry of these GeoDataframes, i.e. buffer the road and river lines.
    'networks' contained the lines that are used to partition space in case of a large dataset. The lines may be the same as the ones used as obstacles, or not.;
    The algorithm returns a geopandas GeoDataframe.

.. code-block:: pycon

  displacement = BuildingDisplacementRandom(network_partitioning=False)
  displaced_gdf = displacement.displace(building_gdf, road_gdf, rivers_gdf)

.. plot:: code/random_displacement.py

Figure 6. A block with buildings displaced because of the width of the road symbol, using the Random Displacement algorithm.

Enrich your data prior to map generalisation
--------------------------------------------

Since the beginning of research on the automation of map generalisation, the necessity for enrichment has been clear. There are properties, structures, which are implicit in the spatial arrangement of geometries in the map. These properties, structures are necessary to make the best decision when generalising the map, and this data enrichment step helps by making these properties, these structures explicit cartographic data.

Extracting implicit geographic structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: compute_boffet_urban_areas(buildings, dilation_size, erosion_size, simplification_distance = 2)

    Computes the urban/built-up areas from a set of buildings, using a method from Boffet (2000). The algorithm computes buffers around each building ('dilation_size') and then merges all buffers.
    The merged areas are then further refined with an erosion ('erosion_size') and a Douglas & Peucker simplification ('simplification_distance').

.. code-block:: pycon

  # compute the built-up areas with a 25 m buffer and a 10 m erosion
  urbanareas = compute_boffet_urban_areas(polygons, 25.0, 10.0)

.. plot:: code/compute_boffet_urban_areas.py

Figure 7. Building polygons converted into built-up areas using the Boffet algorithm.


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
