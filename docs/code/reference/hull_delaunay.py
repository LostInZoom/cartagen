from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import numpy
import geopandas as gpd
from geopandas import GeoDataFrame
import shapely
from shapely.wkt import loads
import cartagen as c4
from shapely import LineString, Polygon, Point

points = [
    loads('POINT (300868.686842791 6176312.581639586)'), 
    loads('POINT (300542.55552210764 6176237.320565582)'), 
    loads('POINT (300599.0013276105 6175973.906806569)'), 
    loads('POINT (300887.502111292 6176067.983149074)'), 
    loads('POINT (301257.53572514426 6176143.244223078)'), 
    loads('POINT (301182.2746511404 6175835.928170895)'), 
    loads('POINT (300956.49142912886 6175735.580072223)'), 
    loads('POINT (300486.10971660475 6175616.4167050505)'), 
    loads('POINT (300009.45624791365 6175873.558707897)'), 
    loads('POINT (299927.92341774283 6176212.233540914)'), 
    loads('POINT (300291.6852754281 6176569.723642433)'), 
    loads('POINT (300178.79366442235 6176136.97246691)'), 
    loads('POINT (301840.8090486741 6175284.0136282)'), 
    loads('POINT (302656.1373503826 6175265.198359699)'), 
    loads('POINT (302442.897640705 6175014.3281130195)'), 
    loads('POINT (302160.66861319053 6175189.937285695)'), 
    loads('POINT (301941.15714734595 6174882.621233513)'), 
    loads('POINT (301991.33119668183 6174750.914354006)'), 
    loads('POINT (302455.441153039 6174644.294499167)'), 
    loads('POINT (302292.3754926973 6174895.164745847)'), 
    loads('POINT (302668.6808627166 6174863.805965012)'), 
    loads('POINT (302762.7572052214 6174794.816647175)'), 
    loads('POINT (303145.3343314076 6174907.708258181)'), 
    loads('POINT (302950.90989023104 6174506.315863494)'), 
    loads('POINT (301226.1769443093 6172292.385936547)'), 
    loads('POINT (300900.0456236259 6172461.723353055)'), 
    loads('POINT (300548.82727827464 6172317.472961214)'), 
    loads('POINT (300736.9799632843 6172066.602714535)'), 
    loads('POINT (300987.8502099638 6172223.39661871)'), 
    loads('POINT (301113.28533330356 6172060.330958368)'), 
    loads('POINT (301445.6884101539 6172085.4179830365)'), 
    loads('POINT (301690.2869006664 6172298.657692714)'), 
    loads('POINT (301602.4823143286 6172499.353890058)'), 
    loads('POINT (301307.70977448014 6172405.277547552)'), 
    loads('POINT (301044.2960154667 6172424.092816054)'), 
    loads('POINT (301069.3830401346 6172204.581350209)'), 
    loads('POINT (301244.99221281026 6172054.059202202)'), 
    loads('POINT (301439.4166539869 6172317.472961214)'), 
    loads('POINT (301408.05787315196 6172179.494325541)'), 
    loads('POINT (301765.5479746703 6172129.320276205)'), 
    loads('POINT (304669.3710799857 6172881.931016244)'), 
    loads('POINT (304587.8382498149 6173208.062336927)'), 
    loads('POINT (305001.7741568361 6173264.50814243)'), 
    loads('POINT (305070.7634746729 6173044.996676586)'), 
    loads('POINT (304456.1313703081 6173113.985994422)'), 
    loads('POINT (304224.0763921296 6172718.865355902)'), 
    loads('POINT (304725.81688548863 6172587.1584763955)'), 
    loads('POINT (305051.948206172 6172838.028723075)'), 
    loads('POINT (305585.04748036596 6173019.909651917)'), 
    loads('POINT (305296.5466966845 6173358.584484935)'), 
    loads('POINT (305020.58942533704 6173571.824194612)'), 
    loads('POINT (304625.4687868168 6173590.639463114)'), 
    loads('POINT (304380.8702963043 6173308.410435599)'), 
    loads('POINT (304311.8809784674 6172894.474528578)'), 
    loads('POINT (304155.0870742927 6172888.202772411)'), 
    loads('POINT (304380.8702963043 6172957.192090248)'), 
    loads('POINT (304719.54512932163 6172775.3111614045)'), 
    loads('POINT (305509.7864063621 6172681.2348189)'), 
    loads('POINT (304236.6199044635 6172605.973744896)'), 
    loads('POINT (304500.03366347705 6172399.005791386)'), 
    loads('POINT (305277.73142818356 6172787.854673739)'), 
]

h1 = c4.hull_delaunay(points, 2000)
h2 = c4.hull_delaunay(points, 3000)
h3 = c4.hull_delaunay(points, 4000)

fig = plt.figure(1, (12, 4))

sub1 = fig.add_subplot(131)
sub1.set_title("a) length=2000", pad=10, family='sans-serif')
sub1.axes.get_xaxis().set_visible(False)
sub1.axes.get_yaxis().set_visible(False)

sub2 = fig.add_subplot(132)
sub2.set_title("b) length=3000", pad=10, family='sans-serif')
sub2.axes.get_xaxis().set_visible(False)
sub2.axes.get_yaxis().set_visible(False)

sub3 = fig.add_subplot(133)
sub3.set_title("c) length=4000", pad=10, family='sans-serif')
sub3.axes.get_xaxis().set_visible(False)
sub3.axes.get_yaxis().set_visible(False)

poly1 = Path.make_compound_path(Path(numpy.asarray(h1.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in h1.interiors])
sub1.add_patch(PathPatch(poly1, facecolor="red", edgecolor='none'))
poly2 = Path.make_compound_path(Path(numpy.asarray(h2.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in h2.interiors])
sub2.add_patch(PathPatch(poly2, facecolor="red", edgecolor='none'))
poly3 = Path.make_compound_path(Path(numpy.asarray(h3.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in h3.interiors])
sub3.add_patch(PathPatch(poly3, facecolor="red", edgecolor='none'))

for p in points:
    c = p.coords[0]
    sub1.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=3)
    sub2.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=3)
    sub3.plot(c[0], c[1], linestyle="", marker='o', color="black", markersize=3)

sub1.autoscale_view()
sub2.autoscale_view()
sub3.autoscale_view()

plt.show()