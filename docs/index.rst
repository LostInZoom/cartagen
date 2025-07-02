.. _index:

.. raw:: html

   <div align="center">

.. image:: img/logo.svg
   :alt: CartAGen - Generalise cartographic objects using advanced algorithms
   :align: center
   :width: 150px

.. raw:: html

   <div style="height: 0; visibility: hidden;">

CartAGen
========

.. raw:: html

   </div>

.. raw:: html

   <br>
   <h1>CartAGen</h1>

A Python library for cartogaphic generalisation using Shapely and GeoPandas

.. image:: https://img.shields.io/readthedocs/cartagen?color=306998&style=flat-square
   :alt: Read the Docs
   :target: https://cartagen.readthedocs.io/en/latest/

.. image:: https://img.shields.io/pypi/v/cartagen?color=306998&style=flat-square
   :alt: PyPI - Version
   :target: https://pypi.org/project/cartagen/

.. image:: https://img.shields.io/github/last-commit/LostInZoom/cartagen?color=ffd43b&style=flat-square
   :alt: GitHub last commit
   :target: https://github.com/LostInZoom/cartagen/commits/main/

.. image:: https://img.shields.io/github/contributors/LostInZoom/cartagen?color=ffd43b&style=flat-square
   :alt: GitHub contributors
   :target: https://github.com/LostInZoom/cartagen/graphs/contributors

.. raw:: html

   <br>
   <br>

   <div style="display: flex; justify-content: center;">
   <div style="margin: 0 10px 0 10px;">

.. image:: img/github.svg
   :alt: Repo GitHub
   :target: https://github.com/LostInZoom/cartagen
   :height: 40px

.. raw:: html

   </div>
   <div style="margin: 0 10px 0 10px;">

.. image:: img/qgis.svg
   :alt: Repo GitHub
   :target: https://github.com/LostInZoom/cartagen-qgis
   :height: 40px

.. raw:: html

   </div>
   <div style="margin: 0 10px 0 10px;">

.. image:: img/jupyter.svg
   :alt: Repo GitHub
   :target: https://github.com/LostInZoom/cartagen-notebooks
   :height: 40px

.. raw:: html

   </div>
   </div>
   </div>

|

**CartAGen** is an open source Python library dedicated to cartogaphic generalisation, published under
the `EUPL-1.2 <https://github.com/IGNF/CartAGen>`_ (European Union Public License).
It is a port of the `Java application, <https://github.com/IGNF/CartAGen>`_
originally developed at IGN France.

It aims at providing a set of tools to generalise spatial data.
Those tools constitutes the foundation on which you have to construct your own
generalisation process. That being said, if you want to learn more about
cartographic generalisation or simply want to familiarize yourself with
the algorithms provided by CartAGen, some Jupyter notebooks are available
`here. <https://github.com/LostInZoom/cartagen-notebooks>`_

CartAGen relies on the usage of the `Shapely <https://github.com/shapely/shapely>`_ geometry objects
and `GeoPandas <https://github.com/geopandas/geopanda>`_ dataset objects.
This approach is based on the idea those libraries are the most commonly used among the
community and provide advantages as powerful spatial operations, measures, indexes, *etc*.
It is recommended for users to have an understanding of those libraries as CartAGen heavily
relies on them.

.. raw:: html

   <div align="center">

Apply map generalisation operations
-----------------------------------

.. plot:: code/index/operation_simplify.py
   :show-source-link: False

.. plot:: code/index/operation_recursive.py
   :show-source-link: False

.. plot:: code/index/operation_blocks.py
   :show-source-link: False

Orchestrate multiple algorithms
-------------------------------

.. plot:: code/index/orchestration_network.py
   :show-source-link: False

.. plot:: code/index/orchestration_agent.py
   :show-source-link: False

Enrich your cartographic data
-----------------------------

.. plot:: code/index/enrichment_boffet.py
   :show-source-link: False

.. plot:: code/index/enrichment_roads.py
   :show-source-link: False

.. plot:: code/index/enrichment_rivers.py
   :show-source-link: False

.. raw:: html

   </div>

.. toctree::
   :caption: User Guide
   :hidden:

   installation
   manual
   changelog
   roadmap
   qgis
   notebooks
   contribution
   faq
   bibliography

.. toctree::
   :caption: API Reference
   :hidden:

   points
   lines
   polygons
   networks
   bends
   processes
   tools