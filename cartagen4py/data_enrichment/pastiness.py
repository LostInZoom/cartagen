import shapely
import geopandas as gpd

def detect_pastiness(line):
    """
    Detect pastiness of a line object. The result is a list of sections of the original line.
    """

    