.. _qgis:

===========
QGIS Plugin
===========

.. .. image:: https://img.shields.io/github/v/release/LostInZoom/cartagen-qgis?color=306998&style=flat-square
..    :alt: QGIS plugin
..    :target: https://github.com/LostInZoom/cartagen-qgis

CartAGen is also available as a QGIS processing plugin |logosub|. This
adds a new toolbox containing the different algorithms of the
Python library. The version of the plugin is tied to a specific
version of CartAGen which is shipped inside the plugin.

.. |logosub| image:: https://img.shields.io/github/v/release/LostInZoom/cartagen-qgis?color=306998&style=flat-square&label=%20
   :target: https://github.com/LostInZoom/cartagen-qgis

Installation
============

The plugin is currently available in the `official QGIS plugin repository <https://plugins.qgis.org/plugins/cartagen4qgis/>`_.
But, as the plugin relies on cartagen, you might need to install cartagen's Python dependencies, namely Numpy, Shapely, GeoPandas,
Matplotlib and NetworkX.

Debian based distributions
--------------------------

Depending on your linux distribution, the installation of CartAGen system-wide can be different.
Please keep in mind that installing system-wide pip packages using this solution will conflict
with the apt packages of your system. Continue at your own risks.

One way to install the CartAGen Python package for QGIS is to use this command outside of a python environment::

    $ pip install numpy shapely geopandas matplotlib networkx

If you are running Debian 12 or above, you might get an error from the system because you are
trying to install the package outside of a virtual environment.
You can bypass this error by using the ``--break-system-package`` flag::

    $ pip install --break-system-package numpy shapely geopandas matplotlib networkx

Windows
-------

To install a Python package for QGIS in Windows (from
`this blog post <https://landscapearchaeology.org/2018/installing-python-packages-in-qgis-3-for-windows/>`_):

#. Open OSGeo4W shell, it should be available in your start menu
#. Type ``py3_env`` in the console (This should print paths of your QGIS Python installation)
#. Use pip to install CartAGen::
    
    python -m pip install numpy shapely geopandas matplotlib networkx