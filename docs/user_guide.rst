==========================
User Guide for cartagen4py
==========================

Apply map generalisation operations
-----------------------------------

Operations for lines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: visvalingam_whyatt(line, area_tolerance)

    Returns a simplified version of the line using the Visvalingam-Whyatt algorithm '(Visvalingam & Whyatt, 1993) <https://www.tandfonline.com/doi/abs/10.1179/000870493786962263?journalCode=ycaj20>'_.
    The area_tolerance is the minimum area of the triangle formed by three consecutive vertices, to keep the middle vertex in the simplified line.

.. code-block:: pycon

  >>> line = LineString([(2, 0), (2, 4), (3, 4), (3, 5), (5, 7)])
  >>> visvalingam_whyatt(line, 1)
  <LINESTRING (2 0, 2 4, 3 5, 5 7)>

.. plot:: code/visvalingam.py

Figure 1. Two polylines simplified with the Visvalingam-Whyatt algorithm.

Operations for polygons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Operations for groups of objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Enrich your data prior to map generalisation
--------------------------------------------

Extracting implicit geographic structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Measures on map features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Apply map generalisation complex processes
------------------------------------------
