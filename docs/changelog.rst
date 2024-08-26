.. _changelog:

Changelog
#########

Version 1.x
===========

1.0.0 (unreleased)
------------------

Coming soon.

Alpha version 0.x
=================

0.3.7 (unreleased)
------------------

Bug fixes:

- Partitioning a polygon by network now assign the polygon using shapely
  `point_on_surface <https://shapely.readthedocs.io/en/stable/reference/shapely.Polygon.html#shapely.Polygon.point_on_surface>`_
  method to keep the point inside the polygon.
- Functions :func:`detect_roundabouts <cartagen4py.detect_roundabouts>`,
  :func:`detect_branching_crossroads <cartagen4py.detect_branching_crossroads>` and
  :func:`detect_dual_carriageways <cartagen4py.detect_dual_carriageways>`
  now return an empty GeoDataFrame if no objects are detected instead of None.
- Functions :func:`collapse_roundabouts <cartagen4py.collapse_roundabouts>`,
  :func:`collapse_branching_crossroads <cartagen4py.collapse_branching_crossroads>`,
  :func:`collapse_dual_carriageways <cartagen4py.collapse_dual_carriageways>` and
  :func:`eliminate_dead_ends <cartagen4py.eliminate_dead_ends>`
  return the input GeoDataFrame if an empty GeoDataFrame is provided as the objects to collapse/eliminate.