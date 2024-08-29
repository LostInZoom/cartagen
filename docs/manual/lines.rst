Simplification
~~~~~~~~~~~~~~

Multiple algorithms for line simplification are available, including:

- :func:`Douglas-Peucker <cartagen.douglas_peucker>`
- :func:`Visvalingam-Whyatt <cartagen.visvalingam_whyatt>`
- :func:`Raposo <cartagen.raposo>`

Those line simplification algorithm are used for different purposes and their computational
time differs. For example, Raposo :footcite:p:`raposo:2013` is mainly used to simplify
natural lines such as rivers. You can read more about them in the API Reference section.

.. plot:: code/manual/lines_simplification.py
  
  Line simplification algorithms

Smoothing
~~~~~~~~~

For now, only one algorithm is available for line smoothing:

- :func:`Gaussian smoothing <cartagen.gaussian_smoothing>`

.. plot:: code/manual/lines_smoothing.py

    Tweaking of the gaussian smoothing parameters