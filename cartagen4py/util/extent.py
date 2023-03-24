import shapely, numpy

# Returns an extent as a numpy array of all provided objects -> [xmin,ymin xmax ymax]
def get_extent_from_objects(*objects):
    # Function to replace bounds coordinates if needeed
    def replace(a, b):
        # Loop through the list b of coordinates
        for i, coord in enumerate(b):
            # If index is 0 or 2 (ie. xmin or ymin) and b coordinates below a coordinates, a value is replaced
            if i in [0, 2] and coord < a[i]:
                a[i] = coord
            # If index is 1 or 3 (ie. xmax or ymax) and b coordinates above a coordinates, a value is replaced
            elif i in [1, 3] and coord > a[i]:
                a[i] = coord
        return a

    # If no objects are passed, throws an error
    if len(objects) < 1:
        raise Exception('Provide at least one object to calculate the extent.')
    
    # Setting the result extent as None
    extent = None

    # Looping in objects
    for object in objects:
        # Calculating the extent of the object
        object_extent = shapely.total_bounds(object).tolist()
        # If the extent varibale is not None, launching the replace function, else, setting the variable as the current extent
        extent = object_extent if extent is None else replace(extent, object_extent)

    return numpy.asarray(extent)