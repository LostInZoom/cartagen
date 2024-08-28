.. _changelog:

Changelog
#########

1.0rc2 (unreleased)
===================

Features:

- Exposed covering algorithms to create convex and concave hull:
    
  - :func:`hull_delaunay <cartagen.hull_delaunay>`
  - :func:`hull_swinging_arm <cartagen.hull_swinging_arm>`

- Exposed :func:`reduce_labelgrid <cartagen.reduce_labelgrid>` as a new point reduction method.

Improvements:

- Renamed point reduction functions:

  - reduce_points_kmeans to :func:`reduce_kmeans <cartagen.reduce_kmeans>`
  - reduce_points_quadtree to :func:`reduce_quadtree <cartagen.reduce_quadtree>`

1.0rc1: Beta Pre-release
========================

The first official beta pre-release of CartAGen.