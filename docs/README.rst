.. meta::
    :author: jberli
    :title: CartAGen - Cartographic generalisation for Python
    :description: Generalise cartographic objects using advanced algorithms

.. only:: pypi

    CartAGen - Generalise cartographic objects using advanced algorithms
    --------------------------------------------------------------------

    .. image:: img/logo.svg
        :alt: CartAGen - Generalise cartographic objects using advanced algorithms
        :align: center
        :width: 150px

.. raw:: html

    <div align="center">

.. only:: readme

   .. image:: img/logo.svg
      :alt: CartAGen - Generalise cartographic objects using advanced algorithms
      :align: center
      :width: 150px

.. only:: html

   .. image:: img/logo.svg
      :alt: CartAGen - Generalise cartographic objects using advanced algorithms
      :align: center
      :width: 150px

.. only:: readme

   .. raw:: html

        <h1>CartAGen</h1>

.. only:: html

    .. raw:: html
        
        <br>
        <span class="h1">CartAGen</span>


A Python library for cartogaphic generalisation using Shapely and GeoPandas

.. only:: readme or pypi

   .. |RTD| replace:: **Documentation**
   .. _RTD: https://cartagen4py.readthedocs.io/en/latest/

   |RTD|_

.. only:: pypi

   |

.. image:: https://img.shields.io/readthedocs/cartagen4py?color=306998&style=flat-square
   :alt: Read the Docs
   :target: https://cartagen4py.readthedocs.io/en/latest/

.. image:: https://img.shields.io/pypi/v/cartagen4py?color=306998&style=flat-square
   :alt: PyPI - Version
   :target: https://pypi.org/project/cartagen4py/

.. image:: https://img.shields.io/github/last-commit/LostInZoom/cartagen4py?color=ffd43b&style=flat-square
   :alt: GitHub last commit
   :target: https://github.com/LostInZoom/cartagen4py/commits/main/

.. image:: https://img.shields.io/github/contributors/LostInZoom/cartagen4py?color=ffd43b&style=flat-square
   :alt: GitHub contributors
   :target: https://github.com/LostInZoom/cartagen4py/graphs/contributors

.. raw:: html

    <br>
    <br>

.. image:: img/github.svg
   :alt: Repo GitHub
   :target: https://github.com/LostInZoom/cartagen4py
   :height: 40px

.. raw:: html

   </div>

|

**CartAGen** is an open source Python library dedicated to cartogaphic generalisation, published under
the `EUPL-1.2. <https://github.com/IGNF/CartAGen>`_ (European Union Public License).
It is a port of the `Java application, <https://github.com/IGNF/CartAGen>`_
originally developed at IGN France.

It aims at providing a set of tools to generalise spatial data.
Those tools constitutes the foundation on which you have to construct your own
generalisation process. That being said, if you want to learn more about
cartographic generalisation or simply want to familiarize yourself with
the algorithms provided by CartAGen, some Jupyter notebooks are available
`here. <https://github.com/LostInZoom/cartagen-notebooks>`_

CartAGen relies on the usage of the geometry formats of `Shapely <https://github.com/shapely/shapely>`_
and the dataset formats of `GeoPandas. <https://github.com/geopandas/geopanda>`_
This approach is based on the idea those libraries are the most commonly used among the
community and provide advantages as powerful spatial operations, measures, indexes, *etc*.