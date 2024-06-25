.. _algorithms:

==========
Algorithms
==========

Lines
-----

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    douglas_peucker
    visvalingam_whyatt
    raposo
    gaussian_smoothing
    offset_curve

Polygons and groups of polygons
-------------------------------

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    building_simplification_ruas
    square_polygons
    morphological_amalgamation
    random_displacement

Networks
--------

Algorithms specific for the generalisation of networks.

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    detect_roundabouts
    collapse_roundabouts
    detect_branching_crossroads
    collapse_branching_crossroads
    detect_dual_carriageways
    collapse_dual_carriageways
    detect_dead_ends
    eliminate_dead_ends

Mountain roads
--------------

Algorithms specific for the generalisation of mountain roads.

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    detect_pastiness
    max_break
    min_break
    accordion
    schematization

Dilation
--------

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    offset_curve
    offset_points
    circle_interpolation