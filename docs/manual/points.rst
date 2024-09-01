Points
======

Reduction
~~~~~~~~~

The reduction of points is a useful method to declutter the map.
It can be used for several types of objects such as toponyms and symbols
but can also be used to generate thematic maps from a really dense
set of points.

Here are the 3 points reduction algorithms available in CartAGen:

- :func:`K-Means reduction <cartagen.reduce_kmeans>`
- :func:`Quad Tree reduction <cartagen.reduce_quadtree>`
- :func:`Label Grid reduction <cartagen.reduce_labelgrid>`

In this section, we will demonstrate those three methods on the Australian cities
from `Natural Earth <https://www.naturalearthdata.com/>`_ because their geographical
repartition is really localised on the southeastern coast of the country.

.. plot:: code/manual/points_reduction_original.py

    The Australian cities from Natural Earth designed for the 1:10m scale

The :func:`K-Means reduction <cartagen.reduce_kmeans>` is used to
reduce the amount of points by forming clusters based on the spatial
proximity of the cities. In selection mode, like the following image shows,
only the city with the largest population is kept.

.. plot:: code/manual/points_reduction_kmeans.py

    The Australian cities reduced by the K-Means method with a shrink ratio of 0.1
    (10% of the cities are kept) and keeping the city with the largest population for each cluster.

On the other hand, the :func:`Quad Tree reduction <cartagen.reduce_quadtree>` method
in selection mode can be used For the same purpose. This method can be interesting
to adress the same issue and is quicker as the amount of points is recursively
divided as the depth of the quad tree raises.

.. plot:: code/manual/points_reduction_quadtree.py

    World's cities reduced using the Quad Tree method with a depth of 3
    by selection on the largest population

Finally, the :func:`Label Grid reduction <cartagen.reduce_labelgrid>` can also be used
in selection mode, like the quad tree, but uses a regular grid of a given size and shape.

.. plot:: code/manual/points_reduction_labelgrid.py

    World's cities reduced using the Label Grid method with an hexagonal grid of 500,000x500,000 meters
    by selection on the largest population

Covering
~~~~~~~~