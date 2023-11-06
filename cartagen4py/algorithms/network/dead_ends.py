import geopandas as gpd

from cartagen4py.utils.graph import *

def eliminate_dead_ends(roads, deadends, length):
    """
    Eliminate dead ends from a road network if they are smaller than the given length.
    """

    # Retrieve crs for output
    crs = roads.crs

    # Convert geodataframe to list of dicts
    roads = roads.to_dict('records')
    deadends = deadends.to_dict('records')

    # Create a list of all the roads geometry of the network
    network = []
    for n in roads:
        network.append(n['geometry'])