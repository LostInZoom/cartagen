Points
======

Reduction
~~~~~~~~~

The reduction of points is a useful method to declutter the map.
It can be used for several types of objects such as toponyms and symbols
but can also be used to generate thematic maps from a really dense
set of points.

K-Means, Quad Tree and Label Grid are 3 different algorithms to generate
clusters of points. Aggregation, selection and simplification are 3 different
methods for using those clusters and are aimed at different use cases.

K-Means-based algorithms :

- :func:`K-Means aggregation <cartagen.kmeans_aggregation>`
- :func:`K-Means selection <cartagen.kmeans_selection>`
- :func:`K-Means simplification <cartagen.kmeans_simplification>`

Quad Tree-based algorithms :footcite:p:`bereuter:2012` :

- :func:`Quad Tree aggregation <cartagen.quadtree_aggregation>`
- :func:`Quad Tree selection <cartagen.quadtree_selection>`
- :func:`Quad Tree simplification <cartagen.quadtree_simplification>`

Label Grid-based algorithms :footcite:p:`grobe:2021` :

- :func:`Label Grid aggregation <cartagen.labelgrid_aggregation>`
- :func:`Label Grid selection <cartagen.labelgrid_selection>`
- :func:`Label Grid simplification <cartagen.labelgrid_simplification>`

In this section, we will show those three algorithms on the Australian cities
from `Natural Earth <https://www.naturalearthdata.com/>`_ because of their repartition
on the island. Indeed, most Australian cities are localised on the southeastern coast of the country,
which offers a good way to test the different point reduction algorithms to see how they keep
the heterogeneity of the data.

.. plot:: code/manual/points_reduction_original.py

    The Australian cities from Natural Earth designed for the 1:10m scale

The K-Means clustering algorithm can be used used to
reduce the amount of points by forming clusters based on the spatial
proximity of the cities. In selection mode, like the following image shows,
only the city with the largest population are kept.

.. plot:: code/manual/points_reduction_kmeans.py

    The Australian cities reduced by the K-Means method with a shrink ratio of 0.1
    (10% of the cities are kept) and keeping the city with the largest population for each cluster.

On the other hand, the Quad Tree based algorithm
in selection mode can be used for the same purpose. This method can be interesting
to adress the same issue and is quicker as the amount of points is recursively
divided while the quad tree depth rises.

.. plot:: code/manual/points_reduction_quadtree.py

    World's cities reduced using the Quad Tree method with a depth of 3
    by selection on the largest population

Finally, the Label Grid algorithm can also be used
in selection mode, like the Quad Tree, but uses a regular grid of a given size and shape.

.. plot:: code/manual/points_reduction_labelgrid.py

    World's cities reduced using the Label Grid method with an hexagonal grid of 500,000x500,000 meters
    by selection on the largest population

Covering
~~~~~~~~

Covering algorithms are used to transform points into surfaces. For now, two methods are
available, which are:

- :func:`Delaunay convex and concave hull <cartagen.hull_delaunay>` :footcite:p:`duckham:2008`
- :func:`Swinging arm algorithm <cartagen.hull_swinging_arm>` :footcite:p:`galton:2006`

The Delaunay method can only create a single polygon, thus, if you want to use
this method, you should first generate clusters using the right library, such as
`scikit-learn clustering, <https://scikit-learn.org/stable/modules/clustering.html>`_
and then create the hulls for each cluster. On the other hand, the swinging arm algorithm
can create multiple polygons as can be seen on the following image.

.. plot:: code/manual/points_covering.py

    Comparison of the two covering methods

Heatmap
~~~~~~~

CartAGen also proposes a :func:`Heatmap <cartagen.heatmap>`
:footcite:p:`wilkinson:2009` function to visualize the density of points within a dataset.
For example, we can display the density of Australian cities weighted by the population.

.. plot:: code/manual/points_heatmap.py

    The density of cities in Australia weighted by the population