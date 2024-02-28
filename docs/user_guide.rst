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


.. method:: gaussian_smoothing(line, sigma, threshold)

    Compute the gaussian smoothing of a set of a LineString. With a smoothing, the shape of the line is simpler but the number of vertices remains the same: the vertices are slightly moved to obtain a smoother shape.
    Sigma is the gaussian filter parameter (the bigger sigma is, the smoother the shape), and threshold is the subsampling parameter, i.e. the step in meters to densify the line prior to the smoothing.
    This code is a port from the GaussianFilter class in the GeOxygene Java library. See p. 119-120 of the book "Algorithmic Foundation of Multi-Scale Spatial Representation" by Z. Li.


.. code-block:: pycon

  line = LineString([(2, 0), (2, 4), (3, 4), (3, 5), (5, 7)])
  smoothed_line = gaussian_smoothing(line, 3.0, 2.0)

.. plot:: code/gaussian_smoothing.py

Figure 3. A polyline smoothed with the Gaussian smoothing algorithm.

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

Figure 4. Four buildings simplified with the Ruas algorithm.

.. method:: square_polygon(polygons, max_iteration=1000, norm_tolerance=0.05, right_tolerance=10, flat_tolerance=10, fixed_weight=5, right_weight=100, flat_weight=50)

    Squares the angles of a polygon using the algorithm from `(Lokhat & Touya, 2016) <https://josis.org/index.php/josis/article/view/72>`_. The algorithm is based on a least squares adjustment, where angles that are almost right or almost flat are adjusted to be exactly right, or exactly flat.

.. code-block:: pycon

  >>> building = Polygon([(0, 0), (0, 10), (9.8, 9.8), (10, 0), (0, 0)])
  >>> square_polygon([building])
  [POLYGON((-0.00002213 -0.00002213, 0.00000159  9.90002291, 9.89999763  9.89999763, 9.90002291  0.00000159, -0.00002213 -0.00002213))]

.. plot:: code/building_squaring.py

Figure 5. Two buildings squared with the algorithm from (Lokhat & Touya, 2016).

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

Figure 6. Buildings amalgamated using the algorithm from Damen et al. (2008).


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
^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: StrokeNetwork(lines, attributeNames)

    This Class contains methods allowing the computation of strokes in a line network representing geographic entities (e.g., roads). Strokes are network segments that follow the perceptual grouping principle of Good Continuity (Gestalt).
    It implements the general method based on both geometry criteria and attributes.
    The initialization of this class is required prior to computing strokes, it includes the precomputing of neighbouring relations between edges of the network.
    The 'lines' parameter gives the geopanda data frame containing the line network geometries and attributes. The 'attributesNames' parameter is the list of attribute names (as str) to be used as a criterion for continuity. 

.. method:: buildStrokes(self, attributeNames,deviatAngle, deviatSum)

    This method computes the stroke in the network using a loop on network features. For each feature that does not already belong to a stroke, it creates a new object of class Stroke and applies the method one side stroke on both sides to find sections that belong to the same stroke as the current section.
    'attributeNames', is the list of attribute names (as str) to be used as a criteria for continuity.
    'deviatAngle', and 'deviatSum' are the thresholds for geometric continuity that respectively represent the maximum angle between two segments at the junction of two sections belonging to the same stroke and the maximum angle between two sections belonging to the same stroke.
    The algorithm updates the attributes strokes of class NetworkStroke that contain the list of strokes in the network. 

.. code-block:: pycon

    data={
      'geometry':[LineString([Point(0, 0),Point(1, 1)]), LineString([Point(1, 1),Point(1, 0)]), LineString([Point(1, 1),Point(2, 2.2)]), LineString([Point(1, 1),Point(2.2, 2)]), LineString([Point(2.2, 2),Point(3, 3)]),],
      'name':["rue de la maison blanche",None,"rue de la maison blanche",None,None],
      'id':[1,2,3,4,5]}
    lines =gpd.GeoDataFrame(data, crs="EPSG:4326")
    sn=StrokeNetwork(lines,['name'])
    sn.buildStrokes(['name'], 45,30)

.. plot:: code/stroke.py

Figure 11. A set of lines with colour depicting the stroke it belongs to using the general algorithm for stroke computation algorithm, with parameters "name", 45 and 30 respectively for attributeNames, deviatAngle and deviatSum.

.. method:: save_strokes_shp(path)

    This algorithm allows to save the computed stroke in a shapefile. The algorithm merges all segments belonging to a stroke in a new entity that has as attribute an id generated as a serial and the comma-separated list of IDs of initial sections used to construct the stroke.
    The 'path' parameter allows us to specify where the output shapefile must be saved.


Stroke computation (for river networks)
^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: RiverStrokeNetwork(lines, attributeNames)

    This Class contains methods allowing the computation of the strokes in a river network. 
    
    :param lines: The geopanda dataframe from which the network must be initialized. It must contain an 'id' column with a unique id (the name is case sensitive). Geometry must be simple LineString (no MultiLineString). The geometry can have a Z value but inconsistencies in Z value may make the stroke research fails.  
    :type lines: GeoDataFrame
    :param attributeNames: List of attribute names to be used as a criteria for continuity.
    :type attributeNames: list[str]
    Strokes are network segments that follow the perceptual grouping principle of Good Continuity (Gestalt). The initialization of this class is required prior to computing strokes, it includes the precomputing of neighbouring relations between edges of the network.


.. method:: buildRiverStrokes(self, attributeNames,deviatAngle, deviatSum)

    This method computes strokes in a RiverStrokeNetwork, add updates its strokes attributes. It can find strokes in complex braided networks.

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
from cartagen4py.data_enrichment import RiverStrokeNetwork
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