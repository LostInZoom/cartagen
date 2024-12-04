import shapely
import geopandas as gpd

from shapely import hausdorff_distance

from cartagen.algorithms.lines.coalescence import coalescence_splitting
from cartagen.algorithms.lines.bends import accordion, schematization

from cartagen.utils.geometry.dilation import dilate_line 
from cartagen.utils.debug import plot_debug, geojson_debug

def galbe(line, width, hausdorff=100):
    """
    Apply GALBE to the provided line.
    """
    split = coalescence_splitting(line, width)

    # plot_debug(*[s['geometry'] for s in split])

    # First treating coalescence on two sides
    for two in [ s for s in split if s['coalescence'] == 2 ]:
        geom1 = two['geometry']
        geom2 = accordion(geom1, width)
        h = hausdorff_distance(geom1, geom2)

        if h > hausdorff:
            geom2 = schematization(geom1)

        print(hausdorff)

    # Second, treating coalescence on one side only
    for one in [ s for s in split if s['coalescence'] == 2 ]:
        pass
    
    # left, right = dilate_line(line, width)

    # splitted = [ s['geometry'] for s in split ]

    # plot_debug(left, right, *splitted)

    # return split