.. _changelog:

Changelog
#########

1.0rc2 (unreleased)
===================

- **New features**:

  #. Exposed covering algorithms to create convex and concave hull:
    
     - :func:`hull_delaunay`
     - :func:`hull_swinging_arm`

  #. Exposed :func:`reduce_labelgrid` as a new point reduction method.

  #. Exposed :func:`strokes_roads` to detect strokes inside a road network.
     This allows an easy strokes calculation by reducing the number of steps.

- **Improvements**:

  #. Renamed point reduction functions:

     - :func:`reduce_points_kmeans` to :func:`reduce_kmeans`.
     - :func:`reduce_points_quadtree` to :func:`reduce_quadtree`.
  
  #. :func:`gaussian_smoothing` can now treat polygons.

- **Bug fixes**:

  #. Fixed the :func:`morphological_amalgamation` issues function caused by:

     - The ``__edge_removal`` function. The function was reworked.
     - The ``straight_line_intersection`` method of the ``Segment`` class crashed
       because of the use of the deprecated numpy array method ``itemset``.
     - The ``Vector2D.from_segment`` method which was fixed.

  #. Fixed bugs in the network enrichment functions: :func:`detect_roundabouts`,
     :func:`detect_branching_crossroads`, :func:`detect_dead_ends`, :func:`detect_dual_carriageways`,
     :func:`rural_traffic`. They now return an empty GeoDataFrame if no entity was detected.

1.0rc1
======

The first official beta pre-release of Cartagen.