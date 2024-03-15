
def measure_block_density(block, buildings, roads):
    """
    Measures the density of the block, i.e. the sum of the areas of the symbols inside the block
    (including the road symbols and the building symbols) divided by the area of the block. The value is usually between 0 and 1 
    but can be bigger than 1 if the block is very dense and scale is small (increasing road symbol size). 
    Parameters
    ----------
    block : a shapely Polygon 
        The polygon represents the measured block, i.e. the face of the road graph.
    buildings : list of shapely Polygon features.
        The geometries of the buildings inside the block.
    roads : list of shapely Polygon features.
        The polygonal symbol of the roads surrounding the block.
        Usually, the symbol width depends on scale.
    """

    symbol_area = 0

    # first add the road symbols
    for road in roads:
        inside = road.intersection(block)
        symbol_area += inside.area
    
    # then add the buildings
    for building in buildings:
        symbol_area += building.area

    return symbol_area / block.area

def measure_block_simulated_density(block, buildings, building_min_size, roads):
    """
    Measures the density of the block, i.e. the sum of the areas of the symbols inside the block
    (including the road symbols and the building symbols) divided by the area of the block. The value is usually between 0 and 1 
    but can be bigger than 1 if the block is very dense and scale is small (increasing road symbol size). 
    Parameters
    ----------
    block : a shapely Polygon 
        The polygon represents the measured block, i.e. the face of the road graph.
    buildings : list of shapely Polygon features.
        The geometries of the buildings inside the block.
    buildings : float.
        The minimal area in m² of the buildings at the target scale of generalisation (depends on scale).
    roads : list of shapely Polygon features.
        The polygonal symbol of the roads surrounding the block.
        Usually, the symbol width depends on scale.
    """

    symbol_area = 0

    # first add the road symbols
    for road in roads:
        inside = road.intersection(block)
        symbol_area += inside.area
    
    # then add the buildings
    for building in buildings:
        symbol_area += max(building.area, building_min_size)

    return symbol_area / block.area