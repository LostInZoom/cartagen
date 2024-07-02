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

Polygons and groups of polygons
-------------------------------

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    building_simplification
    square_polygons
    random_displacement
    morphological_amalgamation
    boffet_areas

Networks
--------

Detection
^^^^^^^^^

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    detect_roundabouts
    detect_branching_crossroads
    detect_dual_carriageways
    detect_dead_ends

Collapse
^^^^^^^^

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    collapse_roundabouts
    collapse_branching_crossroads
    collapse_dual_carriageways
    eliminate_dead_ends

Miscellaneous
^^^^^^^^^^^^^

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    network_faces
    skeletonize_natural
    skeletonize_artificial
    skeletonize_network
    spinalize_polygon
    spinalize_polygons


Bend and bend series
--------------------

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    detect_pastiness
    max_break
    min_break
    accordion
    schematization

Tools
-----

.. currentmodule:: cartagen4py

.. autosummary::
    :toctree: reference/

    dilate_line
    offset_line
    circle_interpolation
    resample_line
    inflexion_points