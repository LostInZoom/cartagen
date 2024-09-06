from shapely.ops import unary_union

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
    building_min_size : float.
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

def mean_building_overlap_rate(buildings, min_sep, roads, road_sizes):
    """
    Measures the mean of the overlapping rates for the buildings of the block. The overlapping rate of a building is the ratio between
    the portions of the building overlapping another symbol (road or building), and the are of the building. The value is between 0 and 1 
    where 0 means that no building overlap another building or the road symbol. 
    Parameters
    ----------
    buildings : list of shapely Polygon features.
        The geometries of the buildings inside the block.
    min_sep : float.
        The minimal distance separating two buildings (in meters).
    roads : list of shapely Polygon features.
        The polygonal symbol of the roads surrounding the block.
        Usually, the symbol width depends on scale.
    road_sizes : list of floats
        a lists of symbol widths for all the sections around the block. 
        The list uses the same order as the sections around the block. 
    """
    mean = 0.0
    nb = 0
    for building in buildings:
        nb+=1
        overlap_geoms = []
        for other in buildings:
            if other == building:
                continue
            if other.distance(building) > min_sep:
                continue
            overlap = other.buffer(min_sep).intersection(building)
            if overlap.is_empty:
                continue
            overlap_geoms.append(overlap)
        
        road_id = 0
        for road in roads:
            symbol_width = road_sizes[road_id]
            if road.distance(building) > (min_sep+symbol_width):
                continue
            overlap = road.buffer(symbol_width).intersection(building)
            if overlap.is_empty:
                continue
            overlap_geoms.append(overlap)
            road_id += 1
        overlap_geom = unary_union(overlap_geoms)
        if overlap_geom.is_empty:
            continue
        mean += overlap_geom.area
    
    if nb > 0:
        return mean / nb
    else:
        return 0.0

