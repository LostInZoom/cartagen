.. _changelog:

Changelog
#########

1.0rc2 (unreleased)
===================

- **New features**:

  #. New function for heatmap creation: :func:`heatmap <cartagen.heatmap>` :pr:`8` (:user:`gowestmen`)

  #. New covering algorithms to create convex and concave hull:
    
     - :func:`hull_delaunay <cartagen.hull_delaunay>` (:user:`gtouya`)
     - :func:`hull_swinging_arm <cartagen.hull_swinging_arm>` :pr:`4` (:user:`Vpech77`)

  #. New :func:`reduce_labelgrid <cartagen.reduce_labelgrid>` function as a new point reduction method.
     :pr:`3` :pr:`6` (:user:`Vpech77`)

  #. New :func:`strokes_roads <cartagen.strokes_roads>` (:user:`ACourtial`) function to detect strokes inside a road network.
     This allows an easy strokes calculation by reducing the number of steps.

- **Improvements**:

  #. Renamed point reduction functions:

     - :func:`reduce_points_kmeans <cartagen.reduce_kmeans>` to :func:`reduce_kmeans <cartagen.reduce_kmeans>`.
     - :func:`reduce_points_quadtree <cartagen.reduce_quadtree>` to :func:`reduce_quadtree <cartagen.reduce_quadtree>`.
  
  #. Every point reduction method, namely :func:`reduce_kmeans <cartagen.reduce_kmeans>`,
     :func:`reduce_quadtree <cartagen.reduce_quadtree>` and :func:`reduce_labelgrid <cartagen.reduce_labelgrid>`
     now takes GeoDataFrame as input and have the same modes available (selection, simplification and aggregation).

  #. :func:`gaussian_smoothing <cartagen.gaussian_smoothing>` can now treat polygons.

- **Bug fixes**:

  #. Fixed the :func:`morphological_amalgamation <cartagen.morphological_amalgamation>` issues function caused by:

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