Lines
=====

Simplification
~~~~~~~~~~~~~~

Multiple algorithms for line simplification are available, including:

- :func:`Douglas-Peucker <cartagen.douglas_peucker>` :footcite:p:`douglas:1973`
- :func:`Visvalingam-Whyatt <cartagen.visvalingam_whyatt>` :footcite:p:`visvalingam:1993`
- :func:`Raposo <cartagen.raposo>`  :footcite:p:`raposo:2013`

Those line simplification algorithm are used for different purposes and their computational
time differs. For example, Douglas-Peucker and Visvalingam-Whyatt are often used to
simplify roads or surfaces borders such as forests.
On the other hand, Raposo and Li-Openshaw are usually applied to simplify
natural lines like the hydrographic network.
You can read more about them in the API Reference section.

.. plot:: code/manual/lines_simplification.py
  
  Line simplification algorithms

Smoothing
~~~~~~~~~

For now, only one algorithm is available for line smoothing, the
:func:`Gaussian smoothing <cartagen.gaussian_smoothing>` :footcite:p:`babaud:1986` :footcite:p:`plazanet:1996`
but this method is very powerful and can be tweaked to best suit your needs.

.. plot:: code/manual/lines_smoothing.py

    Tweaking of the gaussian smoothing parameters