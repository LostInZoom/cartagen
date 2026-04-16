.. _changelog:

Changelog
#########

1.2
===

- **New features**:
   #. New line smoothing algorithm :func:`smooth_wma <cartagen.smooth_wma>`.
   #. New line simplification algorithm :func:`simplify_angular <cartagen.simplify_angular>`.
   #. New line simplification algorithm :func:`simplify_topographic <cartagen.simplify_topographic>`.
   #. New area patches generalisation algorithm :func:`generalise_area_patches <cartagen.generalise_area_patches>`.

1.1
===

- **New features**:
   #. New building regularization algorithm :func:`regularize_building_fer <cartagen.regularize_building_fer>`.
   #. New building typification algorithm :func:`typify_buildings_matching <cartagen.typify_buildings_matching>`.
   #. New line smoothing algorithm :func:`smooth_catmull_rom <cartagen.smooth_catmull_rom>`.
   #. New line smoothing algorithm :func:`smooth_chaikin <cartagen.smooth_chaikin>`.
   #. New line smoothing algorithm :func:`smooth_taubin <cartagen.smooth_taubin>`.
   #. New line smoothing algorithm :func:`smooth_platre <cartagen.smooth_platre>`.
   #. New line simplification algorithm :func:`simplify_whirlpool <cartagen.simplify_whirlpool>`.
   #. New line simplification algorithm :func:`simplify_lang <cartagen.simplify_lang>`.
   #. New line simplification algorithm :func:`simplify_reumann_witkam <cartagen.simplify_reumann_witkam>`.
   #. New line displacement algorithm :func:`beams_displacement <cartagen.beams_displacement>`.
   #. New line displacement algorithm :func:`propagation_crow_flies <cartagen.propagation_crow_flies>`.
   #. New buildings algorithm :func:`building_amalgamation <cartagen.building_amalgamation>` to merge two buildings.
   #. New line simplification algorithms designed to work on polygon boundaries and keep their topology:

      - :func:`boundaries_douglas_peucker <cartagen.boundaries_douglas_peucker>`
      - :func:`boundaries_visvalingam_whyatt <cartagen.boundaries_visvalingam_whyatt>`
      - :func:`boundaries_raposo <cartagen.boundaries_raposo>`
      - :func:`boundaries_li_openshaw <cartagen.boundaries_li_openshaw>`

- **Improvements**:
   #. Renamed :func:`gaussian_smoothing <cartagen.smooth_gaussian>` to :func:`smooth_gaussian <cartagen.smooth_gaussian>`.
   #. Optimization of the :func:`coalescence_splitting <cartagen.coalescence_splitting>` algorithm
      along with all its related algorithms using numpy arrays.
   #. Handle MultiLineString input on line simplification algorithms.
   #. Optimize tesselation algorithm, now the spatial index is only calculated once.
   #. Enhanced :func:`visvalingam_whyatt <cartagen.simplify_visvalingam_whyatt>` using an external code made by
      `Elliot Hallmark. <https://github.com/Permafacture/Py-Visvalingam-Whyatt>`_
   #. Update :func:`strokes_rivers <cartagen.strokes_rivers>` to facilitate their generalization.
   #. Optimized :func:`accordion <cartagen.accordion>` algorithm.
   #. Optimized :func:`schematization <cartagen.schematization>` algorithm.
   #. Optimized :func:`inflexion_points <cartagen.inflexion_points>` algorithm.

- **Bug Fixes**:
   #. Fixed the :func:`collapse_dual_carriageways <cartagen.collapse_dual_carriageways>` algorithm to treat
      polygons connected by their short side.
   #. Fixed the :func:`LeastSquaresMethod <cartagen.LeastSquaresMethod>`, some matrices were wrongly built.

1.0.2
=====

- **New features**:

  #. Added a new :func:`make_planar <cartagen.make_planar>` function that creates a planar
     network from a provided network and keeps a link between input and output.

  #. Added a new :func:`strokes_rivers <cartagen.strokes_rivers>` function that generate strokes
     on river networks.

- **Improvements**:

  #. :func:`network_faces <cartagen.network_faces>` now returns a list of geometry instead
     of a GeometrySequence.

  #. Fixed some inconsistencies within the :func:`dilate_line <cartagen.dilate_line>` function
     that caused random issues with the precision of the vertexes coordinates.
     Now, the function never relies on distances, this is much quicker. Furthermore, wrapped
     line extremities are now accounted for.

  #. The point reduction algorithms now returns the full set of provided points
     with a boolean attribute to know if it was selected. This, of course, is not the case
     in aggregation mode. Furthermore, point reduction algorithms have been divided into
     multiple algorithms in order to be more obvious, easy to use. Read more about that
     on their respective reference page:

     - :func:`kmeans_selection <cartagen.kmeans_selection>`
     - :func:`kmeans_simplification <cartagen.kmeans_simplification>`
     - :func:`kmeans_aggregation <cartagen.kmeans_aggregation>`
     - :func:`labelgrid_selection <cartagen.labelgrid_selection>`
     - :func:`labelgrid_simplification <cartagen.labelgrid_simplification>`
     - :func:`labelgrid_aggregation <cartagen.labelgrid_aggregation>`
     - :func:`quadtree_selection <cartagen.quadtree_selection>`
     - :func:`quadtree_simplification <cartagen.quadtree_simplification>`
     - :func:`quadtree_aggregation <cartagen.quadtree_aggregation>`

- **Bug fixes**:

  #. Fixed an issue in :func:`li_openshaw <cartagen.simplify_li_openshaw>`. :pr:`14` (:user:`gowestmen`)

  #. Fixed an issue in :func:`collapse_branching_crossroads <cartagen.collapse_branching_crossroads>`
     when splitting a line by a point which causes problems of floating points in Shapely.

1.0.1
=====

- **New features**:

  #. Added new kernel methods for the :func:`heatmap <cartagen.heatmap>` function:
     epanechnikov, gaussian, uniform and triangular. :pr:`12` (:user:`gowestmen`)

- **Improvements**:

  #. Renamed :func:`detect_pastiness <cartagen.coalescence_splitting>` to
     :func:`coalescence_splitting <cartagen.coalescence_splitting>` to better
     reflect the usage and litterature.

- **Bug fixes**:

  #. Fixed an issue in the :func:`collapse_branching_crossroads <cartagen.collapse_branching_crossroads>`
     function where :func:`linemerge() <shapely.ops.linemerge>` could input a LineString instead of
     a MultiLineString.

1.0.0
=====

- **New features**:

  #. Exposed three new functions used inside the network enrichment functions:

     - :func:`is_roundabout <cartagen.is_roundabout>`
     - :func:`is_branching_crossroad <cartagen.is_branching_crossroad>`
     - :func:`is_dual_carriageway <cartagen.is_dual_carriageway>`

- **Improvements**:

  #. All undocumented/unwanted functions are now hidden.

  #. Removed all circular import inside the library.

  #. Renamed :class:`Constraint <cartagen.LeastSquaresMethod>` to :class:`LeastSquaresMethod <cartagen.LeastSquaresMethod>`
     to better reflect if usage and enhanced its documentation.

1.0rc2
======

- **New features**:

  #. Added :func:`li_openshaw <cartagen.simplify_li_openshaw>` to simplify lines. (:user:`jberli`)

  #. Added :func:`square_polygon_naive <cartagen.square_polygon_naive>` to square polygons. (:user:`jberli`)

  #. Added :func:`heatmap <cartagen.heatmap>` creation. :pr:`8` (:user:`gowestmen`)

  #. New covering algorithms to create convex and concave hull:
    
     - :func:`hull_delaunay <cartagen.hull_delaunay>` (:user:`gtouya`)
     - :func:`hull_swinging_arm <cartagen.hull_swinging_arm>` :pr:`4` (:user:`Vpech77`)

  #. Added :func:`reduce_labelgrid <cartagen.reduce_labelgrid>` function as a new point reduction method.
     :pr:`3` :pr:`6` (:user:`Vpech77`) :pr:`9` (:user:`gowestmen`)
   
  #. Added :func:`tessellate <cartagen.tessellate>` to create a tesselation of a given shape. This method
     has been taken from :pr:`3` :pr:`6` (:user:`Vpech77`) and wrapped inside a new function.

  #. Added :func:`partition_grid <cartagen.partition_grid>` to partition objects using the new tessellations.

  #. Added :func:`strokes_roads <cartagen.strokes_roads>` (:user:`ACourtial`) function to detect strokes inside a road network.
     This allows an easy strokes calculation by reducing the number of steps.

- **Improvements**:

  #. Renamed point reduction functions:

     - :func:`reduce_points_kmeans <cartagen.reduce_kmeans>` to :func:`reduce_kmeans <cartagen.reduce_kmeans>`.
     - :func:`reduce_points_quadtree <cartagen.reduce_quadtree>` to :func:`reduce_quadtree <cartagen.reduce_quadtree>`.
  
  #. Every point reduction method, namely :func:`reduce_kmeans <cartagen.reduce_kmeans>`,
     :func:`reduce_quadtree <cartagen.reduce_quadtree>` and :func:`reduce_labelgrid <cartagen.reduce_labelgrid>`
     now takes GeoDataFrame as input and have the same modes available (selection, simplification and aggregation).

  #. :func:`gaussian_smoothing <cartagen.gaussian_smoothing>` can now treat polygons.

  #. AGENT rectangle transformation now depends on the minimum rotated rectangle that shares at least
     one edge with the original rectangle. This allows the resulting rectangle to be more aligned
     with the original building.

- **Bug fixes**:

  #. Fixed the :func:`morphological_amalgamation <cartagen.amalgamate_buildings_morphological>` issues function caused by:

     - The ``__edge_removal`` function. The function was reworked.
     - The ``straight_line_intersection`` method of the ``Segment`` class crashed
       because of the use of the deprecated numpy array method ``itemset``.
     - The ``Vector2D.from_segment`` method which was fixed.

  #. Fixed bugs in the network enrichment functions:
     
     - :func:`detect_roundabouts <cartagen.detect_roundabouts>`
     - :func:`detect_branching_crossroads <cartagen.detect_branching_crossroads>`
     - :func:`detect_dead_ends <cartagen.detect_dead_ends>`
     - :func:`detect_dual_carriageways <cartagen.detect_dual_carriageways>`
     - :func:`rural_traffic <cartagen.rural_traffic>`
     
     They now return an empty GeoDataFrame if no entity was detected.

  #. Fixed a bug in :class:`PointSetQuadTree` where negative coordinates could cause problems.

1.0rc1
======

The first official beta pre-release of CartAGen.