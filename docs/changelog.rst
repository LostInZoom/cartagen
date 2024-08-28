.. _changelog:

Changelog
#########

1.0rc2 (unreleased)
===================

Features:

- Exposed covering algorithms to create convex and concave hull:
    
  - :func:`hull_delaunay`
  - :func:`hull_swinging_arm`

- Exposed :func:`reduce_labelgrid` as a new point reduction method.

Improvements:

- Renamed point reduction functions:

  - reduce_points_kmeans to :func:`reduce_kmeans`.
  - reduce_points_quadtree to :func:`reduce_quadtree`.

Bug fixes:

- Fixed the :func:`morphological_amalgamation` function caused by:

  - The ``__edge_removal`` function. The function was reworked.
  - The ``straight_line_intersection`` method of the ``Segment`` class crashed
    because of the use of the deprecated numpy array method ``itemset``.
  - The ``Vector2D.from_segment`` method which was fixed.
    


1.0rc1
======

The first official beta pre-release of Cartagen.