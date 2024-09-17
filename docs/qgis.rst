.. _qgis:

===========
QGIS Plugin
===========

CartAGen is also available as a QGIS processing plugin. This
adds a new toolbox containing the different algorithms of the
Python library. The version of the plugin is tied to a specific
version of CartAGen which is shipped inside the plugin.

Installation
============

For now, the plugin is not available in the official QGIS repository,
so you can download it from the other `GitHub repository. <https://github.com/LostInZoom/cartagen-qgis>`_
In the :guilabel:`release` folder, download the latest `zip` file, decompress it
and copy/paste the :guilabel:`cartagen4qgis` folder in:

* On linux (QGIS 3.32)::

    ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/

* On Windows (QGIS 3.28)::

    C:\Users\<username>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\


Now, inside QGIS, open the :guilabel:`Plugins --> Manage and Install Plugins` menu item.
Activate the :guilabel:`CartAGen4QGIS` plugin. That's it!
A new toolbox containing CartAGen algorithms should appear
inside your :guilabel:`Processing Toolbox`.