.. image:: https://img.shields.io/readthedocs/cartagen4py?color=green
   :alt: Read the Docs
   :target: https://cartagen4py.readthedocs.io/en/latest/

.. image:: https://img.shields.io/pypi/v/cartagen4py?color=green
   :alt: PyPI - Version
   :target: https://pypi.org/project/cartagen4py/

.. image:: https://img.shields.io/github/last-commit/LostInZoom/cartagen4py?color=blue
   :alt: GitHub last commit
   :target: https://github.com/LostInZoom/cartagen4py

.. image:: https://img.shields.io/github/contributors/LostInZoom/cartagen4py?color=blue
   :alt: GitHub contributors
   :target: https://github.com/LostInZoom/cartagen4py/graphs/contributors

**CartAGen**: A cartographic generalisation Python library
==========================================================

.. image:: docs/img/logo.svg
   :height: 76px
   :alt: CartAGen logo
   :align: left
   :target: https://github.com/LostInZoom/cartagen4py

CartAGen is an open source Python library dedicated to cartogaphic generalisation, published under
the European Union Public `EUPL-1.2. <https://github.com/IGNF/CartAGen>`_ license.
It is a port of the `Java application, <https://github.com/IGNF/CartAGen>`_
originally developed at IGN France.

It aims at providing you with a set of tools to generalise your spatial data.
Those tools constitutes the foundation on which you have to construct your own
generalisation process. That being said, if you want to learn more about
cartographic generalisation or simply want to familiarize yourself with
the algorithms provided by CartAGen, some Jupyter notebooks are available
`here. <https://github.com/LostInZoom/cartagen-notebooks>`_

CartAGen relies on the usage of the geometry formats of `Shapely <https://github.com/shapely/shapely>`_
and the dataset formats of `GeoPandas. <https://github.com/geopandas/geopanda>`_
This approach is based on the idea those libraries are the most commonly used among the
community and provide advantages as powerful spatial operations, measures, indexes, *etc*.