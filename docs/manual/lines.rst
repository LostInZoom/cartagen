Lines
=====

Simplification
~~~~~~~~~~~~~~

Multiple algorithms for line simplification are available, including:

- :func:`Douglas-Peucker <cartagen.simplify_douglas_peucker>` :footcite:p:`douglas:1973`
- :func:`Lang <cartagen.simplify_lang>`  :footcite:p:`lang:1969`
- :func:`Li-Openshaw <cartagen.simplify_li_openshaw>`  :footcite:p:`li:1993`
- :func:`Raposo <cartagen.simplify_raposo>`  :footcite:p:`raposo:2013`
- :func:`Reumann-Witkam <cartagen.simplify_reumann_witkam>`  :footcite:p:`reumann:1974`
- :func:`Visvalingam-Whyatt <cartagen.simplify_visvalingam_whyatt>` :footcite:p:`visvalingam:1993`
- :func:`Whirlpool <cartagen.simplify_whirlpool>`  :footcite:p:`dougenik:1979`

Those line simplification algorithm are used for different purposes and their computational
time differs. For example, Douglas-Peucker and Visvalingam-Whyatt are often used to
simplify roads or surfaces borders such as forests.
On the other hand, Raposo, Li-Openshaw and Whirlpool are usually applied to simplify
natural lines like the hydrographic network.
You can read more about them in the API Reference section.

.. plot:: code/manual/lines_simplification.py
  
  Line simplification algorithms

Smoothing
~~~~~~~~~

Multiple algorithms for polyline smoothing are available, including:

- :func:`Gaussian smoothing <cartagen.smooth_gaussian>` :footcite:p:`babaud:1986` :footcite:p:`plazanet:1996`
- :func:`PLATRE smoothing <cartagen.smooth_platre>` :footcite:p:`fritsch:1998`
- :func:`Taubin smoothing <cartagen.smooth_taubin>` :footcite:p:`taubin:1995`

The gaussian smoothing is the most well-known smoothing algorithm for polylines, it attenuates
inflexions in the line but isn't specialized in anything. On the other hand, PLATRE is designed
to preserve sharp turns and attenuate minor bends while Taubin is used to prevent shrinkage of bends
that is characteristic of the gaussian smoothing.

.. plot:: code/manual/lines_smoothing.py

    Polyline smoothing algorithms

CartAGen also contains two smoothing algorithms that relies on the iterative creation of vertexes:

- :func:`Catmull-Rom smoothing <cartagen.smooth_catmull_rom>` :footcite:p:`catmull:1974` :footcite:p:`barry:1988`
- :func:`Chaikin smoothing <cartagen.smooth_chaikin>` :footcite:p:`chaikin:1974` :footcite:p:`wu:2004`

Catmull-Rom smoothing preserves the original vertexes of the line and link them with spline curves while
Chaikin cut the corners of the original line.

.. plot:: code/manual/lines_smoothing_interpolation.py

    Comparison between Catmull-Rom and Chaikin