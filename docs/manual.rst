.. _manual:

====================
CartAGen User Manual
====================

:Authors: Guillaume Touya, Justin Berli, Azelle Courtial
:Copyright:
  This work is licensed under a `European Union Public Licence v1.2`__.

.. __: https://eupl.eu/

.. _intro:

Introduction
============

Please note that this user guide is still under construction and some functions of the library are still not enough documented or not documented at all. 
Feel free to contact us if you need documentation for a specific function.

Generalisation operations
-------------------------

Lines simplification
^^^^^^^^^^^^^^^^^^^^

Several algorithms for line simplification are available, including Douglas-Peucker(1973), Visvalingam-Whyatt and Raposo hexagon-based simplification.

.. plot:: code/manual/line_simplification.py

Polygons
^^^^^^^^

cartagen4py contains algorithms that process any type of polygons, and others specific to some types of map polygons, such as buildings. Only the algorithms that process one polygon at a time are documented in this section.

.. method:: building_simplification(building, edge_threshold, parallel_limit = 20 * pi / 180, orthogonal_limit = 20 * pi / 180)

    Returns a simplified version of the building polygon using the algorithm from Anne Ruas (1988). The algorithm was later used to simplify buildings in the AGENT project.
    The 'edge_threshold' is the minimum length of an edge between two vertices of the building to be removed. The 'parallel_limit' and 'orthogonal_limit' parameters define to what extent two edges are considered parallel or orthogonal.

.. code-block:: pycon

  >>> building = Polygon([(0, 0), (0, 10), (2, 10), (2, 9), (10, 9), (10, 0), (0, 0)])
  >>> building_simplification(building, 2.5)
  <POLYGON ((0 0, 0 9.5, 10 9.5, 10 0, 0 0))>

.. plot:: code/building_simplification.py

Figure 4. Four buildings simplified with the Ruas algorithm.

.. method:: square_polygon(polygons, max_iteration=1000, norm_tolerance=0.05, right_tolerance=10, flat_tolerance=10, fixed_weight=5, right_weight=100, flat_weight=50)

    Squares the angles of a polygon using the algorithm from `(Lokhat & Touya, 2016) <https://josis.org/index.php/josis/article/view/72>`_. The algorithm is based on a least squares adjustment, where angles that are almost right or almost flat are adjusted to be exactly right, or exactly flat.

.. code-block:: pycon

  >>> building = Polygon([(0, 0), (0, 10), (9.8, 9.8), (10, 0), (0, 0)])
  >>> square_polygon([building])
  [POLYGON((-0.00002213 -0.00002213, 0.00000159  9.90002291, 9.89999763  9.89999763, 9.90002291  0.00000159, -0.00002213 -0.00002213))]

.. plot:: code/building_squaring.py

Figure 5. Two buildings squared with the algorithm from (Lokhat & Touya, 2016).

Groups of objects
^^^^^^^^^^^^^^^^^

.. method:: morphological_amalgamation(buildings, buffer_size, edge_length)

    Amalgamates a group of building polygons using morphological operators, with the algorithm presented in `(Damen et al., 2008) <https://www.semanticscholar.org/paper/High-Quality-Building-Generalization-by-Extending-Damen-Kreveld/b64618584b3ae3725da7eeb5a545d1580e5f2113>`_. 
    The algorithm chains morphological opening and closing to amalgamate close buildings into a larger building area.
    The 'buffer_size' is parameter used for the opening and closing operations. The 'edge_length' gives the length of edges that are later simplified in the amalgamated polygons.

.. code-block:: pycon

  >>> buildings = [Polygon([(1, 0), (9, 0), (9, 6), (1, 6), (1, 0)]),Polygon([(10, 0), (17, 0), (17, 6), (10, 6), (10, 0)])]
  >>> morphological_amalgamation(buildings, 1.0, 1.0)
  <POLYGON ((1.207 1.983, 2.547 5.885, 16.768 4.282, 15.42 0.148, 1.207 1.983))>

.. plot:: code/building_amalgamation.py

Figure 6. Buildings amalgamated using the algorithm from Damen et al. (2008).


.. class:: RandomDisplacement(max_trials=25, max_displacement=10, network_partitioning=False, verbose=False)

    This algorithm displaces buildings that overlap with each other and/or other features. The algorithm was never published but was available in CartAGen. It is an iterative process that selects the building with the most overlaps, and then pushes slightly the building in a random direction. If the overlaps are reduced, the displacement is commited and a new iteration starts. If the overlaps are worse, the displacement is backtracked, and another one is tried.
    The 'max_trials' parameter gives the maximum number of random displacements tried on one building. The 'max_displacement' parameter is the maximum distance a building is allowed to move. For large datasets, the algorithm can work on smaller partitions, using the 'network_partitioning' parameter.
    The name of the class mentions buildings but other objects can be similarly displaced, as long as GeoDataframe with polygons is provided.

.. method:: displace(self, buildings, roads, rivers, *networks)

    This method displaces the buildings with roads and rivers acting as obstacles for the buildings, i.e. the algorithm minimises the overlaps between buildings and with the geometries contained in those two collections.
    'buildings', 'roads', and 'rivers' are geopandas GeoDataframes, not arrays of geometries. If you want to avoid overlaps with road and river symbols, you need to provide polygons as the main geometry of these GeoDataframes, i.e. buffer the road and river lines.
    'networks' contained the lines that are used to partition space in case of a large dataset. The lines may be the same as the ones used as obstacles, or not.;
    The algorithm returns a geopandas GeoDataframe.

.. code-block:: pycon

  displacement = RandomDisplacement(network_partitioning=False)
  displaced_gdf = displacement.displace(building_gdf, road_gdf, rivers_gdf)

.. plot:: code/random_displacement.py

Figure 7. A block with buildings displaced because of the width of the road symbol, using the Random Displacement algorithm.

.. method:: kmeans_point_set_reduction(points, shrink_ratio, centroid_option = False)

    This algorithm reduces a set of points to a smaller set of points that is representative of the initial set. The algorithm uses a K-Means clustering to reduce the set to a number of clusters that corresponds to the shrinking ratio parameter.
    The 'shrink_ratio' parameter can vary between 0 (all points are removed) and 1 (all points are kept).
    Two options are possible: either keeping one of the initial points to replace a cluster (default option) or replace the cluster by its centroid.

.. code-block:: pycon

  >>> points = [Point(1,1), Point(1,2), Point(0,1), Point(2,1), Point(2,2), Point(5,5), Point(8,10), Point(10,10), Point(10,8), 
              Point(16,10), Point(16,9), Point(14,11)]
  >>> kmeans_point_set_reduction(points, 0.25)
  [<POINT (2.0 2.0)>, <POINT (10.0 10.0)>, <POINT (16.0 10.0)>]

.. plot:: code/kmeans_reduction.py

Figure 8. A set of points reduced to 25% of its initial amount, with the K-Means reduction algorithm.

.. method:: quadtree_point_set_reduction(points, depth, mode='simplification', attribute = "")

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
  >>> quadtree_point_set_reduction(points, 0.25)
  [<POINT (1 2)>, <POINT (10 8)>, <POINT (10 10)>, <POINT (14 11)>]

.. plot:: code/quadtree_reduction.py

Figure 9. A set of points reduced to depth 2 of the quadtree, with the selection mode. The selected points are displayed in red.


Road network
^^^^^^^^^^^^

Those functions are used to generalized specific features inside a road network. Those tools are used in conjonction with the
data enrichment tools.

.. method:: collapse_roundabouts(roads, roundabouts, crossroads=None, maximum_diameter=None)

    This function collapse detected roundabouts inside a road network.
    **It is recommended to detect both roundabouts and branching crossroads before collapsing them, this approach yields better results.**
    Returns the new road network as a geopandas GeoDataFrame.
    
    :param roads: The road network where roundabouts will be collapsed.
    :type roads: geopandas GeoDataFrame of LineStrings
    :param roundabouts: The polygons representing the faces of the network detected as roundabouts.
    :type roundabouts: geopandas GeoDataFrame of Polygons
    :param crossroads: The polygons representing the faces of the network detected as branching crossroads. This allows incoming branching crossroads on roundabouts to be collapsed as well. 
    :type crossroads: geopandas GeoDataFrame of Polygons, optional
    :param maximum_diameter: The diameter below which roundabouts are not collapsed.
    :type maximum_diameter: float, optional
    
.. code-block:: pycon

    # Detect roundabouts using default parameters
    roundabouts = detect_roundabouts(roads)

    # Detect branching crossroads using default parameters
    branching = detect_branching_crossroads(roads)

    # Collapse roundabouts with default parameters
    collapsed = collapse_roundabouts(roads, roundabouts, crossroads=branching)

.. plot:: code/collapse_roundabouts.py

.. method:: collapse_branching_crossroads(roads, crossroads, roundabouts=None, maximum_area=None)

    This function collapse detected branching crossroads inside a road network.
    **It is recommended to detect both roundabouts and branching crossroads before collapsing them, this approach yields better results.
    Then, the collapsing of branching crossroads connected to a roundabout is conducted using the roundabout collapsing algorithm.**
    Returns the new road network as a geopandas GeoDataFrame.
    
    :param roads: The road network where branching crossroads will be collapsed.
    :type roads: geopandas GeoDataFrame of LineStrings
    :param crossroads: The polygons representing the faces of the network detected as branching crossroads. 
    :type crossroads: geopandas GeoDataFrame of Polygons
    :param roundabouts: The polygons representing the faces of the network detected as roundabouts. This allows roundabouts to be collapsed at the same time.
    :type roundabouts: geopandas GeoDataFrame of Polygons, optional
    :param maximum_area: The area, in square meter, below which branching crossroads are collapsed. Default value is set to None. 
    :type maximum_area: float, optional
    
.. code-block:: pycon

    # Detect roundabouts using default parameters
    roundabouts = detect_roundabouts(roads)

    # Detect branching crossroads using default parameters
    branching = detect_branching_crossroads(roads)

    # Collapse branching crossroads with default parameters
    collapsed = collapse_branching_crossroads(roads, branching, roundabouts=roundabouts)

.. plot:: code/collapse_branching_crossroads.py

.. method:: collapse_dual_carriageways(roads, carriageways, sigma=None, propagate_attributes=None)

    This function collapse detected dual carriageways inside a road network.
    Returns the new road network as a geopandas GeoDataFrame.
    
    :param roads: The road network where dual carriageways will be collapsed.
    :type roads: geopandas GeoDataFrame of LineStrings
    :param carriageways: The polygons representing the faces of the network detected as dual carriageways.
    :type carriageways: geopandas GeoDataFrame of Polygons
    :param sigma: If not None, apply a gaussian smoothing to the collapsed dual carriageways to avoid jagged lines that can be created during the TIN skeleton creation.
    :type sigma: float, optional
    :param propagate_attributes: Propagate the provided list of column name to the resulting network. The propagated attribute is the one from the longest line.
    :type propagate_attributes: list of str, optional
    
.. code-block:: pycon

    # Detect branching crossroads using default parameters
    carriageways = detect_dual_carriageways(roads)

    # Collapse branching crossroads with default parameters
    collapsed = collapse_dual_carriageways(roads, carriageways, sigma=2)

.. plot:: code/collapse_dual_carriageways.py

.. method:: eliminate_dead_ends(roads, deadends, length, keep_longest=True)

    This function eliminates dead ends inside a road network if the length of their main component is below a given threshold.
    If the dead end is simple (i.e. just one road), the main component is the road.
    If the dead end contains multiple ramification of roads, the main component represents the road between the entry and the longest ramification.
    If the dead end contains inner network faces (i.e. enclosed roads), the main component represents the longest of the shortest paths between the entry and all the nodes of the dead ends.
    Returns the roads network without the unwanted dead ends as a GeoDataFrame.

    :param roads: The GeoDataFrame containing the dead ends as LineString geometries.
    :type roads: geopandas GeoDataFrame
    :param deadends: The LineString representing the roads of the network detected as dead ends.
    :type deadends: geopandas GeoDataFrame of Polygons
    :param length: The length below which dead ends are eliminated.
    :type length: float
    :param keep_longest: If set to true, in case of complex dead end, keep the main component (c.f. description) if above the provided length.
    :type keep_longest: boolean, optional

.. code-block:: pycon

    # Detect dead ends using default parameters
    deadends = detect_dead_ends(network)

    # Eliminate dead ends using a length threshold of 250
    eliminated = eliminate_dead_ends(network, deadends, 250)

.. plot:: code/collapse_dead_ends.py

Mountain roads
^^^^^^^^^^^^^^

Those functions are used to generalized mountain roads with a lot of bends.

.. method:: detect_pastiness(line, tolerance, cap_style='flat', quad_segs=8)

    Detect pastiness of a line object.
    Returns a list of dictionary as { "paste": **paste**, "geometry": **geometry** } where **paste** represents the number of conflicts (0 when no
    conflicts are detected, 1 when a conflict exists on one side only, two when conflicts are on both side of the line) and **geometry**
    is the shapely geometry of the line.
    This algorithm subdivide the provided line into multiple chunks, thus modifying the geometry,
    it is not a data enrichment function stricto sensu.
    
    :param line: The line to detect the pastiness from.
    :type line: shapely LineString
    :param tolerance: The Tolerance of the offset used to detect the pastiness.
    :type tolerance: float
    :param cap_style: The type of caps at the start and end of the provided line. Possible values are 'round' or 'flat'. Default to 'flat'.
    :type cap_style: str, optional
    :param quad_segs: The number of point allowed per circle quadrant when interpolating points using round method. Default to 8.
    :type quad_segs: int, optional
    
.. code-block:: pycon

    # Detect pastiness using a tolerance of 60 metres and default parameters
    pastiness = detect_pastiness(line, 60)

.. plot:: code/mountain_pastiness.py

Detection of the pastiness of a line (the width represent the number of conflict as described in the method description)


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

Stroke computation (in general)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Strokes are network segments that follow the perceptual grouping principle of Good Continuity (Gestalt).

.. class:: StrokeNetwork(lines, attributeNames)

    This Class contains methods allowing the computation of strokes in a line network representing geographic entities (e.g., roads). 
    
    :param lines: The geopanda dataframe from which the network must be initialized. It must contain an 'id' column with a unique id (the name is case sensitive). Geometry must be simple LineString (no MultiLineString). The geometry can have a Z value but inconsistencies in Z value may make the stroke research fails.  
    :type lines: GeoDataFrame
    :param attributeNames: List of attribute names to be used as a criteria for continuity.
    :type attributeNames: list[str]
    
    The initialization of this class is required prior to computing strokes, it includes the precomputing of neighbouring relations between edges of the network.

.. method:: buildStrokes(self, attributeNames,deviatAngle, deviatSum)

    This method computes the strokes in a Strokenetwork using a loop on network features, and updates its strokes attribute.
    
    :param self: The network in which we expect to compute strokes
    :type self: StrokeNetwork
    :param attributeNames: List of attribute names to be used as a criteria for continuity.
    :type attributeNames: list[str]
    :param deviatAngle: Thresholds for the maximum angle between two segments at the junction of two sections belonging to the same stroke.
    :type deviatAngle: float
    :param deviatSum: Thresholds for the maximum angle between two sections belonging to the same stroke.
    :type deviatAngle, deviatSum: float
    
     For each feature that does not already belong to a stroke, it creates a new object of class Stroke and applies the method one side stroke on both sides to find sections that belong to the same stroke as the current section.

.. code-block:: pycon

	from shapely.geometry import LineString, Point
	import geopandas as gpd
	from cartagen4py.enrichment import StrokeNetwork
	import matplotlib.pyplot as plt

	data={'geometry':
        [LineString([Point(0, 0),Point(1, 1)]),
        LineString([Point(1, 1),Point(1, 0)]),
        LineString([Point(1, 1),Point(2, 2.2)]),
        LineString([Point(1, 1),Point(2.2, 2)]),
        LineString([Point(2.2, 2),Point(3, 3)])],
        'name':["rue A",None,None,"rue A","rue A"],
        'id':[1,2,3,4,5]}
	lines =gpd.GeoDataFrame(data, crs="EPSG:4326")

	sn=StrokeNetwork(lines,['name'])

	sn.buildStrokes(['name'], 45,30)
	array=sn.reconstruct_strokes()
	gdf = gpd.GeoDataFrame(array,  columns = ['id','geom',"section"],crs="epsg:2154",geometry="geom")   
	gdf.plot('id')
	plt.show()


.. plot:: code/stroke.py

Figure 11. A set of lines with colour depicting the stroke it belongs to using the general algorithm for stroke computation algorithm, with parameters "name", 45 and 30 respectively for attributeNames, deviatAngle and deviatSum.

.. method:: save_strokes_shp(path)

    This method save the computed stroke in a shapefile. 
    
    :param path: The access path to the file in which the stroke must be recorded
    :type path: str
    
    The saved shapefile is made with segment belonging to a unique stroke merged in a geometries  the attributes of each geometries are an id (generated as a serial) and the comma-separated list of IDs of initial sections used to construct the stroke.


Stroke computation (for river networks)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: RiverStrokeNetwork(lines, attributeNames)

    This Class contains methods allowing the computation of the strokes in a river network. 
    
    :param lines: The geopanda dataframe from which the network must be initialized. It must contain an 'id' column with a unique id (the name is case sensitive). Geometry must be simple LineString (no MultiLineString). The geometry can have a Z value but inconsistencies in Z value may make the stroke research fails.  
    :type lines: GeoDataFrame
    :param attributeNames: List of attribute names to be used as a criteria for continuity.
    :type attributeNames: list[str]
    The initialization of this class is required prior to computing strokes, it includes the precomputing of neighbouring relations between edges of the network.


.. method:: buildRiverStrokes(self, attributeNames,deviatAngle, deviatSum)

    This method computes strokes in a RiverStrokeNetwork, and updates its strokes attributes. It can find strokes in complex braided networks.

    :param self: The RiverNetwork in which we expect to compute strokes
    :type self: RiverStrokeNetwork
    :param attributeNames: List of attribute names to be used as a criteria for continuity.
    :type attributeNames: list[str]
    :param deviatAngle: Thresholds for the maximum angle between two segments at the junction of two sections belonging to the same stroke.
    :type deviatAngle: float
    :param deviatSum: Thresholds for the maximum angle between two sections belonging to the same stroke.
    :type deviatAngle, deviatSum: float
    Stroke are computed from sources to sink while computing Strahler order.
    First, it identifies each source as a departure for a stroke adds them to the downstream section list and sets its Strahler order to 1.
    Then the main loop runs through the downstream section list, pops the current element and adds the next section in its stroke (if exists).


.. code-block:: pycon
    from shapely.geometry import LineString, Point
    import geopandas as gpd
    from cartagen4py.enrichment import RiverStrokeNetwork
    import matplotlib.pyplot as plt

    data={'geometry':
        [LineString([Point(1,4),Point(1, 3)]),
         LineString([Point(1.5,3.5),Point(1, 3)]),
         LineString([Point(1, 3),Point(1, 2.4)]),
         LineString([Point(1, 2.4),Point(0.8, 1.8),Point(0.9, 1.5)]),
         LineString([Point(1, 2.4),Point(1.2, 2.1)]),
         LineString([Point(1.2, 2.1),Point(0.9, 1.5)]),
         LineString([Point(0.9, 1.5),Point( 1.2,0.6)]),
         LineString([Point(1.2, 2.1),Point( 1.2,0.6)]),
         LineString([Point( 1.2,0.6),Point(1.1, 0.3)]),
         LineString([Point(1.1, 0.3),Point(1, 0)]),
         LineString([Point(0.5, 2),Point(1.1, 0.3)])],
        'id':[1,2,3,4,5,6,8,9,10,11,12]}
    lines =gpd.GeoDataFrame(data, crs="EPSG:4326")

    sn=RiverStrokeNetwork(lines,None)

    sn.buildRiverStrokes([], 45,30)
    array=sn.reconstruct_strokes()
    gdf = gpd.GeoDataFrame(array,  columns = ['id', 'geom',"strahler"],crs="epsg:4326",geometry="geom")

    a=gdf.plot('id')
    plt.show()

    b=gdf.plot('strahler')
    plt.show()

.. plot:: code/riverstroke.py

Figure 12. A river network with color depicting the stroke. 
Figure 13. A river network with color depicting the Horton order : purple =1; yellow=2.

.. method:: save_strokes_shp(path)

    This method save the computed stroke in a shapefile. 
    
    :param path: The access path to the file in which the stroke must be recorded
    :type path: str
    
    The saved shapefile is made with segment belonging to a unique stroke merged in a geometries  the attributes of each geometries are an id (generated as a serial) and the comma-separated list of IDs of initial sections used to construct the stroke.


Road network enrichment
^^^^^^^^^^^^^^^^^^^^^^^

Those functions characterize different specificities of a road network by interpreting its layout and shape.

.. method:: detect_roundabouts(network, area_threshold=40000, miller_index=0.95)

    This function detects roundabouts inside a road network.
    Returns the polygons representing the roundabouts extent as a geopandas GeoDataFrame.
    Returns None if no roundabouts where found.
    
    :param network: The GeoDataFrame containing the road network as LineString geometries.
    :type network: geopandas GeoDataFrame
    :param area_threshold: The area (in square meters) above which the object is not considered a roundabout. The default value is set to 40000.
    :type area_threshold: int, optional
    :param miller_index: Index of compactess that determines if the shape is round or not. The default value is set to 0.97.
    :type miller_index: float, optional
    
.. code-block:: pycon

    # Detect roundabouts using default parameters
    detect_roundabouts(network)

.. plot:: code/detect_roundabouts.py

.. method:: detect_branching_crossroads(roads, area_threshold=2500, maximum_distance_area=0.5, roundabouts=None, allow_middle_node=True, middle_angle_tolerance=10.0, allow_single_4degree_node=True)

    This function detects branching crossroads inside a road network.
    **Although the roundabouts parameter is optional, it is recommended to detect roundabouts before branching crossroads to help their characterization.**
    Returns a GeoDataFrame of polygons representing their extents.
    
    :param network: The GeoDataFrame containing the road network as LineString geometries.
    :type network: geopandas GeoDataFrame
    :param area_threshold: The area (in square meters) above which the object is not considered a branching crossroads. The default value is set to 2500.
    :type area_threshold: int, optional
    :param area_threshold: The maximum distance area between the actual polygon and the triangle formed by the 3 nodes connecting the junction to the rest of the network. The default value is set to 0.5.
    :type maximum_distance_area: float, optional
    :param roundabouts: The polygons representing the network faces considered as roundabouts. If provided, it offers a better detection of branching crossroads. The default value is set to None.
    :type roundabouts: geopandas GeoDataFrame of Polygons, optional
    :param allow_middle_node: If set to True, allow 4 nodes to form the crossroads, but each must have a degree of 3 and the 'middle' node must have an angle of 180°. Default value set to False.
    :type allow_middle_node: boolean, optional
    :param middle_angle_tolerance: If allow_middle_node is set to True, indicate an angle tolerance in degree for the fourth node of the crossroad to be considered the middle node. Default value is set to 10.0.
    :type middle_angle_tolerance: float, optional
    :param allow_single_4degree_node: If set to True, allow one and only one node to have a degree of 4. Default value set to False.
    :type allow_single_4degree_node: float, optional
    
.. code-block:: pycon

    # Detect branching crossroads using default parameters
    detect_branching_crossroads(network)

.. plot:: code/detect_branching_crossroads.py

.. method:: detect_dual_carriageways(roads, importance=None, value=None, concavity=0.85, elongation=6.0, compactness=0.12, area=60000.0, width=20.0, huber=16)

    Detects dual carriageways within a road network. Dual carriageways are derived from the network faces.
    Return road separators as GeoDataFrame polygons.
    
    :param roads: The GeoDataFrame containing the road network as LineString geometries.
    :type roads: geopandas GeoDataFrame
    :param importance: The attribute name of the data on which road importance is based. Default value is set to None which means every road is taken for the network face calculation.
    :type importance: str, optional
    :param value: Maximum value of the importance attribute. Roads with an importance higher than this value will not be taken. Default value is set to None.
    :type value: int, optional
    :param concavity: Minimum concavity above which the face is a dual carriageway. It represents the factor between the polygon surface and its convex hull surface. Default value is set to 0.85.
    :type concavity: float, optional
    :param elongation: Minimum elongation above which the face is a dual carriageway. It represents the ratio between the length and the width of the minimum rotated rectangle containing the polygon. Default value is set to 6.0.
    :type elongation: float, optional
    :param compactness: Maximum compactness below which the face is a dual carriageway. (4*pi*area/perimeter^2)Default value is set to 0.12.
    :type compactness: float, optional
    :param area: Area factor to detect very long motorways. Do not change if you don't know what you are doing. Default value is set to 60000.0.
    :type area: float, optional
    :param width: Minimum width above which the face is a dual carriageway. It represents the width of the minimum rotated rectangle containing the polygon. Default value is set to 20.0.
    :type width: float, optional
    :param huber: Huber width for long motorways. Do not change. Default value is set to 16.
    :type huber: int, optional

.. code-block:: pycon

    # Detect dual carriageways using default parameters
    detect_dual_carriageways(network)

.. plot:: code/detect_dual_carriageways.py

.. method:: detect_dead_ends(roads, outside_faces=True)

    This function detects dead ends inside a road network.
    Returns the roads detected as dead-ends. Return None if none were found.

    :param roads: The GeoDataFrame containing the road network as LineString geometries.
    :type roads: geopandas GeoDataFrame
    :param outside_faces: If set to true, detects dead-ends on the network faces located on the border.
    :type outside_faces: boolean, optional
    
Five attributes are added:

* **face**: the id of the network face it belongs to.
* **deid**: the id of the dead end group inside a given face.
* **connected**: set to 'true' if the dead end group is connected to the network.
* **root**: set to true if the road section is the root of the dead end group, i.e. the section connecting the dead end group to the road network.
* **hole**: set to true if the road section touches a hole inside the dead end group.

.. code-block:: pycon

    # Detect dead ends using default parameters
    detect_dead_ends(network)

.. plot:: code/detect_dead_ends.py


Apply map generalisation complex processes
------------------------------------------

AGENT model
^^^^^^^^^^^^^^^^^^^^^^^^
This user guide is not meant to fully explain the principles of the AGENT model, and how it works. If you are not familiar with the AGENT model, please read the scientific papers describing this model:
- `<http://icaci.org/files/documents/ICC_proceedings/ICC2001/icc2001/file/f13041.pdf>`_
- `<http://dx.doi.org/10.1016/b978-008045374-3/50016-8>`_
- `<https://hal.inria.fr/IFSTTAR/hal-01682131v1>`_
Though it is a tutorial for the JAVA CartAGen platform, this `webpage <https://ignf.github.io/CartAGen/docs/tuto_agents.html>`_ also contains complementary information on how to trigger agent-based map generalisation.

Micro agents
=============
You may want to use micro agents only, i.e. one cartographic feature such as a building that generalises itself without consideration for the other cartographic features. To generalise micro agents, you have to follow these steps:
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

Least squares adjustment
^^^^^^^^^^^^^^^^^^^^^^^^