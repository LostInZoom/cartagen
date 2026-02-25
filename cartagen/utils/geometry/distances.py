import shapely

def distance_area(geom1, geom2):
    """
    Calculcate the distance area between two geometries.
    distance area = area(intersection) / area(union)
    """
    if geom1.geom_type != 'Polygon' or geom2.geom_type != 'Polygon':
        raise Exception('Distance area calculation can be calculated only on Polygon geometries.')

    intersection = shapely.intersection(geom1, geom2)
    union = shapely.union(geom1, geom2)

    if intersection.is_empty == False and union.is_empty == False:
        return 1 - (intersection.area / union.area)
    else:
        return None

# def group_intersecting(geoms):
#     """
#     Group intersecting geometries in nested lists.
#     """
#     done = []
#     groups = []
#     groupsgeom = []
#     for i, p1 in enumerate(geoms):
#         if i not in done:
#             done.append(i)
#             group = [ i ]
#             groupgeom = p1
#             for j, p2 in enumerate(geoms):
#                 if j not in done:
#                     if shapely.intersects(groupgeom, p2):
#                         done.append(j)
#                         group.append(j)
#                         groupgeom = shapely.union(groupgeom, p2)
#             groups.append(group)
#             groupsgeom.append(groupgeom)
#     return groups

def group_intersecting(polygons):
    """
    Group intersecting geometries in nested lists.
    """
    n = len(polygons)
    parent = list(range(n))
    
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    
    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y:
            parent[root_x] = root_y
    
    # Create a spatial index
    tree = shapely.STRtree(polygons)
    
    # For each polygon, find candidates for intersection
    for i, poly in enumerate(polygons):
        # Spatial request
        candidats_idx = tree.query(poly)
        
        for j in candidats_idx:
            # Union only if current polygon index is below candidate
            if i < j and polygons[i].intersects(polygons[j]):
                union(i, j)
    
    # Regroup indexes
    groups = {}
    for i in range(n):
        root = find(i)
        if root not in groups:
            groups[root] = []
        groups[root].append(i)
    
    return list(groups.values())