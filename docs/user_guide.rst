.. _user-guide:

==========================
User Guide for cartagen4py
==========================

Please note that this user guide is still under construction and some functions of the library are still not enough documented or not documented at all. 
Feel free to contact us if you need documentation for a specific function.

Apply map generalisation operations
-----------------------------------

Operations for lines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: visvalingam_whyatt(line, area_tolerance)

    Returns a simplified version of the line using the Visvalingam-Whyatt algorithm '(Visvalingam & Whyatt, 1993) <https://www.tandfonline.com/doi/abs/10.1179/000870493786962263?journalCode=ycaj20>'_.
    The 'area_tolerance' is the minimum area of the triangle formed by three consecutive vertices, to keep the middle vertex in the simplified line.

.. code-block:: pycon

  >>> line = LineString([(2, 0), (2, 4), (3, 4), (3, 5), (5, 7)])
  >>> visvalingam_whyatt(line, 1)
  <LINESTRING (2 0, 2 4, 3 5, 5 7)>

.. plot:: code/visvalingam.py

Figure 1. Two polylines simplified with the Visvalingam-Whyatt algorithm.


.. method:: raposo_simplification(line, initial_scale, final_scale, centroid=True, tobler=False)

    Returns a simplified version of the line using the Raposo algorithm '(Raposo, 2013) <http://dx.doi.org/10.1080/15230406.2013.803707>'_.
    The algorithm uses an hexagonal tessallation, with a size related to the final scale, and it only retains one vertex per hexagonal cell.
    Be careful, it uses the scale as parameter. If the 'centroid' parameter is ''True'', the vertices inside an hexagon cell are replaced by the centroid of the cell; if it is ''False'', they are replaced by the nearest vertex to the centroid of the cell.
    The Raposo algorithm is dedicated to the simplification of natural lines such as rivers, lakes or forests.

.. code-block:: pycon

  line = LineString([(2, 0), (2, 4), (3, 4), (3, 5), (5, 7)])
  simplified_line = raposo_simplification(line, 10000.0, 50000.0)

.. plot:: code/raposo.py

Figure 2. Two polylines simplified with the Raposo algorithm.

Operations for polygons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^:

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

Operations for groups of objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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

Figure 4. Building polygons converted into built-up areas using the Boffet algorithm.


Measures on map features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
