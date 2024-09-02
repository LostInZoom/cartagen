.. _installation-qgis:

============
Installation
============

Install CartAGen package for QGIS
---------------------------------

Debian based distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on your linux distribution, the installation of CartAGen system-wide can be different.
Please keep in mind that installing system-wide pip packages using this solution will conflict with the apt packages of your system. Continue at your own risks.

One way to install the CartAGen Python package for QGIS is to use this command outside of a python environment::

    $ pip install cartagen

If you are running Debian 12 or above, you might get an error from the system because you are trying to install the package outside of a virtual environment.
You can bypass this error by using the ``--break-system-package`` flag::

    $ pip install --break-system-package cartagen

Windows
^^^^^^^

To install a Python package for QGIS in Windows (from `this blog post <https://landscapearchaeology.org/2018/installing-python-packages-in-qgis-3-for-windows/>`_):

#. Open OSGeo4W shell, it should be available in your start menu
#. Type ``py3_env`` in the console (This should print paths of your QGIS Python installation)
#. Use pip to install CartAGen::
    
    python -m pip install cartagen


Install CartAGen4QGIS
---------------------

Now that CartAGen is installed and available for QGIS, you can install the QGIS plugin simply by copying the folder
``cartagen4qgis`` from the `GitHub repository <https://github.com/LostInZoom/cartagen>`_ and paste it inside:

* On linux (QGIS 3.32)::

    ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/

* On Windows (QGIS 3.28)::

    C:\Users\<username>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\


Now, inside QGIS, open the :guilabel:`Plugins --> Manage and Install Plugins` menu item. Activate the :guilabel:`CartAGen4QGIS` plugin. That's it!
A new toolbox containing CartAGen algorithms should appear inside your :guilabel:`Processing Toolbox`.