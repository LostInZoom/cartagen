.. _contribution:

Contribution
############

In the spirit of open-source development, you are encouraged to enhance this Python
library with you own contributions by forking the `GitHub repository. <https://github.com/LostInZoom/cartagen>`_
You can contribute in different ways, either by implementing new algorithms, resolving known bugs,
proposing enhancements for already-present algorithms, *etc.*

From your fork of the repository, you can simply propose a pull request and it will be
approved if no conflicts or code breakings bugs are detected.

Development guidelines
======================

Here are some guidelines to help you implement your generalisation algorithms in the Python library.
As a general rule of thumb, we ask you to:

#. Use GeoPandas GeoDataFrame as input and output if your algorithm relies on attributes
   or need multiple geometries as context for generalisation. But please never relies on
   column names or specific CRS.

#. Use Shapely Geometry types as input and output if your algorithm relies on single geometries, try to
   anticipate input geometry types to avoid breakage (MultiGeometry for example).

#. Try to use :class:`STRTree <shapely.STRTree>` when applicable, this can greatly improve your algorithm speed
   when dealing with a lot of geometry. 

#. If your algorithm can be decomposed into multiple functions which have their place as standalone algorithms,
   implement both of them separately. For example, :func:`partition_grid <cartagen.partition_grid>` relies on
   :func:`tessellate <cartagen.tessellate>`, thus, those two functions have been implemented separately.

As you might have seen if you've taken a look at the repository, algorithms are not
developped homogeneously throughout the library, this is because different people
with different coding habits have contributed. We don't expect you to have the same
coding habits as us, so we only ask that you try to debug as thoroughly as possible
your algorithms before submitting them. We also ask you to create docstrings using the
following template:

.. code-block:: Python
   
   def algorithm1(param1, param2, param3=False, *args, **kwargs):
      """
      A short description of the algorithm.

      This algorithm was proposed by/is described in :footcite:p:`name:date`
      (Don't hesitated to add references inside the docs/bibliography.bib file).
      A long description of the algorithm where you can describe
      how it works.

      Parameters
      ----------
      param1 : GeoDataFrame of Polygon
         Description of the first parameter which is a GeoDataFrame with Polygon geometries.
      param2 : float
         Description of the second parameter which is a float.
      param3 : bool, optional
         Description of the third parameter which is a boolean.

      Returns
      -------
      result : GeoDataFrame of Polygon
         If needed, a description of the result (added columns, modified geometries, etc.)

      Warning
      -------
      Description of a warning message if needed.
      
      See Also
      --------
      algorithm2 :
         The short description of a related algorithm.

      Notes
      -----
      Notes to user using your algorithm if needed.

      # Add the references section if you have :footcite:p:`name:date` inside the description.

      References
      ----------
      .. footbibliography::
      """

Contributors
============

.. contributors:: LostInZoom/cartagen
    :avatars:
    :order: DESC