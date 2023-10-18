def attributes_from_longest(lines, attributes):
    """
    Returns the attribute value which is has the longest representation from the provided list of entries.
    Returns a list of x entry (same as attributes entry) where each entry the value with the longest line representation.
    """

    # Dict to store the attributes with the associated value
    result = {}

    # Loop through the list of wanted attributes
    for a in attributes:
        # Dict to store the attribute value along with the total length of lines having it
        lengths = {}
        # Loop through lines
        for l in lines:
            # Check to make sure the attribute exists
            if a not in l:
                raise Exception("The provided attribute {0} does not exist.".format(a))

            # Retrieve the attribute value and the length of the line
            value = l[a]
            length = l["geometry"].length

            # If the key value already exists, add the length to the total
            if value in lengths:
                lengths[value] += length
            # Else, create the dict entry
            else:
                lengths[value] = length
        
        # Select the max length from the attributes value dict and create a new result dict entry
        result[a] = max(lengths, key=lengths.get)

    return result