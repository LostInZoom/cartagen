import shapely, numpy
import geopandas as gpd
import pandas as pd

from cartagen.utils.attributes import attributes_from_longest
from cartagen.utils.network.roads import Crossroad
from cartagen.utils.network.faces import NetworkFace
from cartagen.utils.geometry.skeletonization import SkeletonTIN
from cartagen.utils.geometry.spinalization import spinalize_polygon
from cartagen.utils.geometry.line import extend_line_with_point, merge_linestrings, merge_linestrings

def collapse_dual_carriageways(roads, carriageways, propagate_attributes=None, sigma=None, blend_smoothing=False):
    """
    Collapse dual carriageways using a TIN skeleton.

    This algorithm proposed by Thom :footcite:p:`thom:2005`
    collapses the network faces considered as dual carriageways
    using a skeleton calculated from a Delaunay triangulation.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        The road network.
    carriageways : GeoDataFrame of Polygon
        The polygons representing the faces of the network detected as dual carriageways.
    propagate_attributes : list of str, optional
        Propagate the provided list of column name to the resulting network.
        The propagated attribute is the one from the longest line.
    sigma : float, optional
        If not None, apply a gaussian smoothing to the collapsed dual carriageways to
        avoid jagged lines that can be created during skeletonization.
    blend_smoothing : bool, optional
        If True and a sigma is defined, apply the gaussian smoothing to the skeleton once it has
        been blended within the provided network. This will modify the network outside the extent
        of the detected dual carriageway polygon. If false and a sigma is defined, apply the
        gaussian smoothing only to the skeleton before blending with the network.

    See Also
    --------
    detect_dual_carriageways : 
        Detect dual carriageways inside a road network.
    skeletonize_network :
        Blends a TIN skeleton inside a road network.

    References
    ----------
    .. footbibliography::
    """

    def merge_overlapping_lines(gdf):
        """
        Union overlapping lines. Attributes of the first line are kept.
        Force the result to be a simple linestring.
        
        Parameters
        ----------
        gdf : GeoDataFrame
            GeoDataFrame containing lines
        
        Returns
        -------
        GeoDataFrame
            GeoDataFrame containing unioned lines where applicable
        """
        gdf = gdf.copy().reset_index(drop=True)
        
        # Spatial index
        tree = shapely.STRtree(gdf.geometry)
        
        to_remove = set()
        merged_geoms = {}
        
        for i in range(len(gdf)):
            if i in to_remove:
                continue
                
            geom_i = gdf.geometry.iloc[i]
            
            # Find candidates
            candidates = tree.query(geom_i, predicate='intersects')
            
            lines_to_merge = [geom_i]
            indices_to_merge = [i]
            
            for j in candidates:
                if j <= i or j in to_remove:
                    continue
                
                geom_j = gdf.geometry.iloc[j]
                
                # Check if lines overlaps
                if geom_i.overlaps(geom_j):
                    lines_to_merge.append(geom_j)
                    indices_to_merge.append(j)
                    to_remove.add(j)
            
            # If they do, union
            if len(lines_to_merge) > 1:
                # Create progressive union
                merged = lines_to_merge[0]
                for line in lines_to_merge[1:]:
                    merged = merged.union(line)
                
                # Force simple linestring
                if merged.geom_type == 'MultiLineString':
                    # Try linemerge
                    merged = shapely.ops.linemerge(merged)
                    
                    # If still MultiLineString, takes longest component
                    if merged.geom_type == 'MultiLineString':
                        merged = max(merged.geoms, key=lambda x: x.length)
                
                merged_geoms[i] = merged
        
        # Apply unioned geometries
        for idx, geom in merged_geoms.items():
            gdf.loc[idx, 'geometry'] = geom
        
        # Drop old geometries
        return gdf[~gdf.index.isin(to_remove)].reset_index(drop=True)

    def group_carriageways(polygons):
        """
        Group intersecting geometries in nested lists.
        Separate intersections that are only points (0-dimensional).
        Excludes polygons that are connected along their length.
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
        
        # Track polygons to exclude (connected by length or by a single point)
        exclude_indices = set()
        
        # For each polygon, find candidates for intersection
        for i, poly in enumerate(polygons):
            # Spatial request
            candidats_idx = tree.query(poly)
            for j in candidats_idx:
                # Union only if current polygon index is below candidate
                if i < j and polygons[i].intersects(polygons[j]):
                    # Check if intersection is more than just a point
                    intersection = polygons[i].intersection(polygons[j])

                    if intersection.geom_type == 'MultiLineString':
                        line = shapely.ops.linemerge(intersection.geoms)
                    elif intersection.geom_type == 'LineString':
                        line = intersection
                    else:
                        # If it's not a line, exclude and continue the loop
                        # It's probably a point, but it also exclude polygons if there are topological issues
                        exclude_indices.add(i)
                        exclude_indices.add(j)
                        continue

                    face_i = NetworkFace(polygons[i])
                    face_j = NetworkFace(polygons[j])
                    
                    # Check if either polygon is connected by its length
                    if abs(line.length - face_i.length) < abs(line.length - face_i.width):
                        # Shared edge is the long side of polygon i
                        exclude_indices.add(i)
                        exclude_indices.add(j)
                        continue
                    
                    if abs(line.length - face_j.length) < abs(line.length - face_j.width):
                        # Shared edge is the long side of polygon j
                        exclude_indices.add(i)
                        exclude_indices.add(j)
                        continue
                
                    union(i, j)
        
        # Regroup indexes
        groups = {}
        for i in range(n):
            if i not in exclude_indices:
                root = find(i)
                if root not in groups:
                    groups[root] = []
                groups[root].append(i)

        return list(groups.values())
    
    # Retrieve crs for output
    crs = roads.crs

    if len(carriageways) == 0:
        return roads

    carriageways = carriageways.to_dict('records')

    polygons = [ c['geometry'] for c in carriageways ]
    # Create groups of intersecting carriageways
    groups = group_carriageways(polygons)

    # Transform potential multi geometries into simple geometries
    roads = roads.explode(index_parts=False).reset_index(drop=True)

    faces = []
    remove = []
    for group in groups:
        # Create the unioned polygon using the whole group
        polygon = shapely.ops.unary_union([ polygons[i] for i in group ])
        faces.append(polygon)
        # Look for road contained within the polygon, they need to be removed
        contains = roads.sindex.query(polygon, predicate='contains').tolist()
        remove.extend(contains)

    removed = roads.loc[remove]

    # This list will store the indexes of roads to throw away
    skeletons = []
    roads = roads.drop(remove)
    roads = roads.reset_index(drop=True)
    network = roads.geometry.tolist()

    removed.to_file('cartagen/data/removedinside.geojson')
    roads.to_file('cartagen/data/noinside.geojson')

    remove = []
    crossroads = []
    for face in faces:
        # Create the crossroad object without the removed roads
        crossroad = Crossroad(network, roads.sindex, face)
        crossroads.append(crossroad)

        # Stores internal roads for removal
        remove.extend(crossroad.internals)

        # If there are more than one external road to the crossroad
        if len(crossroad.externals) > 1:
            # Retrieve the id of the external roads that have not been changed by the conversion to a crossroad object
            unchanged = crossroad.get_unchanged_roads('externals')

            # Check if attribute propagation is wanted
            attributes = None
            if propagate_attributes is not None:
                # Retrieve roads entry to calculate longest attributes
                proads = []
                boundary = face.boundary
                # Loop through the original roads id of the crossroad
                for rid in crossroad.original:
                    # Check that the road is contained by the boundary of the carriageway or overlaps it
                    if boundary.contains(network[rid]) or boundary.overlaps(network[rid]):
                        # If so, append it
                        proads.append(roads.iloc[rid])
                
                # Retrieve the longest attributes value
                attributes = attributes_from_longest(proads, propagate_attributes)
            
            # Retrieve incoming roads, i.e. the external network of the crossroad
            incoming = []
            # Looping through external roads
            for ext in crossroad.externals:
                original = None
                # Retrieve the geometry of the line
                egeom = crossroad.network[ext]
                # Looping through unchanged network
                for u in unchanged:
                    # Retrieve geometry
                    ugeom = network[u]
                    # If the line equals an unchanged line
                    if shapely.equals(egeom, ugeom):
                        # Set original to be the road object with its attributes
                        original = roads.iloc[u].to_dict()
                        remove.append(u)
                # If an unchanged line match has been found, add the object to the list
                if original is not None:
                    incoming.append(original)
                # Else, create a new object without attributes
                else:
                    # Check if attribute propagation is wanted
                    if attributes is not None:
                        incoming.append({ **attributes, **{"geometry": egeom} })
                    else:
                        incoming.append({ "geometry": egeom })
            
            # Calculate the skeleton
            skeleton = SkeletonTIN(face)
            skeleton.add_incoming_lines(incoming)
            skeleton.create_network()
            skeleton.blend(attributes, sigma=sigma, blend_smoothing=blend_smoothing)
            skeletons.append(gpd.GeoDataFrame(skeleton.blended, crs=crs))

            # elif method == 'spinalization':
            #     entries = crossroad.get_entry_points()
            #     spines = spinalize_polygon(face, densify=densify, sigma=sigma, entries=entries)

            #     print(entries)
            #     print(incoming)

            #     skeleton = []
            #     remove_inc = []
            #     for spine in spines:
            #         modified = False
            #         for i, inc in enumerate(incoming):
            #             igeom = inc['geometry']
            #             merged = None
            #             # Try to merge incoming with current spine
            #             try:
            #                 merged = merge_linestrings(spine, igeom)
            #             except:
            #                 continue

            #             inc['geometry'] = merged
            #             skeleton.append(inc)
            #             # remove this incoming line
            #             remove_inc.append(i)
            #             modified = True

            #         if not modified:
            #             skeleton.append({"geometry": spine})
                
            #     if len(skeleton) > 0:
            #         skeletons.append(gpd.GeoDataFrame(skeleton, crs=crs))

            # Storing the original geometries of the crossroad
            remove.extend(crossroad.original)

    removed = roads.loc[remove]
    roads = roads.drop(remove)

    cleaned = pd.concat([roads, *skeletons], ignore_index=True)
    cleaned = cleaned.drop_duplicates(subset='geometry')
    cleaned = merge_overlapping_lines(cleaned)

    incoming = []
    for c in crossroads:
        for g in c.externals:
            incoming.append({'geometry': c.network[g]})
    inc = gpd.GeoDataFrame(incoming, crs=3857)
    inc.to_file('cartagen/data/incoming.geojson')
            
    cleaned.to_file('cartagen/data/cleaned.geojson')
    removed.to_file('cartagen/data/removedoustide.geojson')

    return cleaned



    # # This list will store the indexes of roads to throw away
    # originals = []
    # # This will store the new collapsed geometries
    # collapsed = []

    # # Stores carriageways that touches either on their long side or short side
    # longside = []
    # shortside, shortside_list = [], []
    # # Stores carriageways that touches on a single point
    # pointside, pointside_list = [], []

    # # This will store future skeletons
    # skeletons = []

    # # Here, calculate the tin skeleton if applicable
    # # --------------------------------------------------------
    # for cid, carriageway in enumerate(carriageways):
    #     # Get the geometry of the face
    #     polygon = carriageway['geometry']

    #     # Calculate the crossroad object
    #     crossroad = Crossroad(network, tree, polygon)

    #     # If there are more than one external road to the crossroad
    #     if len(crossroad.externals) > 1 and cid not in longside:
    #         # Retrieve the id of the external roads that have not been changed by the conversion to a crossroad object
    #         unchanged = crossroad.get_unchanged_roads('externals')

    #         # Check if attribute propagation is wanted
    #         attributes = None
    #         if propagate_attributes is not None:
    #             # Retrieve roads entry to calculate longest attributes
    #             proads = []
    #             boundary = polygon.boundary
    #             # Loop through the original roads id of the crossroad
    #             for rid in crossroad.original:
    #                 # Get the geometry of the road
    #                 rgeom = roads[rid]["geometry"]
    #                 # Check that the road is contained by the boundary of the carriageway or overlaps it
    #                 if boundary.contains(rgeom) or boundary.overlaps(rgeom):
    #                     # If so, append it
    #                     proads.append(roads[rid])
                
    #             # Retrieve the longest attributes value
    #             attributes = attributes_from_longest(proads, propagate_attributes)

    #         # Retrieve incoming roads, i.e. the external network of the crossroad
    #         incoming = []
    #         # Looping through external roads
    #         for ext in crossroad.externals:
    #             original = None
    #             # Retrieve the geometry of the line
    #             egeom = crossroad.network[ext]
    #             # Looping through unchanged network
    #             for u in unchanged:
    #                 # Retrieve geometry
    #                 ugeom = network[u]
    #                 # If the line equals an unchanged line
    #                 if shapely.equals(egeom, ugeom):
    #                     # Set original to be the road object with its attributes
    #                     original = roads[u]
    #             # If an unchanged line match has been found, add the object to the list
    #             if original is not None:
    #                 incoming.append(original)
    #             # Else, create a new object without attributes
    #             else:
    #                 # Check if attribute propagation is wanted
    #                 if attributes is not None:
    #                     incoming.append({ **attributes, **{"geometry": egeom} })
    #                 else:
    #                     incoming.append({ "geometry": egeom })

    #         # Calculate the skeleton
    #         skeleton = SkeletonTIN(polygon)
    #         skeleton.add_incoming_lines(incoming)
    #         skeleton.create_network()
    #         skeleton.blend(attributes, sigma=sigma)
    #         skeletons.append(skeleton)

    #         # Storing the original geometries of the crossroad
    #         originals.extend(crossroad.original)

    #     else:
    #         # Add None to keep the order of the indexes
    #         skeletons.append(None)

    # # Here, carriageways connected by their short sides are treated
    # # --------------------------------------------------------------

    # # Get the bone junction between the two provided entry points
    # def get_junction(points, skeleton):
    #     # Get the first junction starting from an entry point
    #     def get_next_joint(point, network):
    #         for ni, n in enumerate(network):
    #             # Retrieve start and end point of bone
    #             start, end = n.coords[0], n.coords[-1]
    #             if start == point:
    #                 return ni, end
    #             elif end == point:
    #                 return ni, start
    #         return None, None

    #     ni1, j1 = get_next_joint(points[0], skeleton.network)
    #     ni2, j2 = get_next_joint(points[1], skeleton.network)

    #     # If that junction is the same point, return the point
    #     if j1 == j2:
    #         return [ni1, ni2], j1
    #     # Else, return None
    #     else:
    #         return None, None

    # # Remove 'fake' incoming lines from both skeletons object
    # # i.e. incoming lines from one skeleton objects that is contained or overlap the other polygon geometry
    # def remove_fake_incoming(skeleton1, skeleton2):
    #     def keep_loop(incoming, boundary):
    #         keep = []
    #         for i in incoming:
    #             if boundary.contains(i['geometry']) or boundary.overlaps(i['geometry']):
    #                 continue
    #             keep.append(i)
    #         return keep

    #     b1 = skeleton1.polygon.boundary
    #     b2 = skeleton2.polygon.boundary

    #     return keep_loop(skeleton1.incoming, b2), keep_loop(skeleton2.incoming, b1)

    # # Retrieve incoming lines
    # def retrieve_incoming_lines(incoming, entries, boundary):
    #     iresult = []
    #     for i, inc in enumerate(incoming):
    #         add = False
    #         # Add only if the line intersects an entry point
    #         for e in entries:
    #             if shapely.intersects(inc['geometry'], shapely.Point(e)):
    #                 add = True
    #         if add:
    #             iresult.append([i, inc])
    #     return iresult
    
    # # Get the 'middle' line of a skeleton network
    # def get_middle_line(junction, shortentries, skeleton):
    #     for ni, n in enumerate(skeleton.network):
    #         start, end = n.coords[0], n.coords[-1]
    #         if junction == start:
    #             if end not in shortentries:
    #                 return ni, n, 'start'
    #         elif junction == end:
    #             if start not in shortentries:
    #                 return ni, n, 'end'
    #     return None, None, None

    # # Treating short side connections between carriageways
    # for shorts in shortside:
    #     # Retrieve carriageways index and the shortside geometry
    #     cid1, cid2, shortline = shorts[0], shorts[1], shorts[2]

    #     if cid1 in longside or cid2 in longside:
    #         continue
        
    #     # Create a list containing entry points intersecting the short side
    #     shortentries = list(filter(lambda x: shapely.intersects(x, shortline), skeletons[cid1].entries))

    #     # Remove those entries from skeletons list of entries
    #     skeletons[cid1].entries = [x for x in skeletons[cid1].entries if x not in shortentries]
    #     skeletons[cid2].entries = [x for x in skeletons[cid2].entries if x not in shortentries]

    #     # Converts this list to tuples of coordinates
    #     shortentries = [x.coords[0] for x in shortentries]

    #     # Retrieve the junction between both entries inside both skeletons
    #     diamond1, j1 = get_junction(shortentries, skeletons[cid1])
    #     diamond2, j2 = get_junction(shortentries, skeletons[cid2])

    #     # Continue only if both junctions were found
    #     if j1 is not None and j2 is not None:
    #         # Remove 'fake' incoming lines
    #         skeletons[cid1].incoming, skeletons[cid2].incoming = remove_fake_incoming(skeletons[cid1], skeletons[cid2])

    #         # Retrieve incoming lines
    #         incoming = retrieve_incoming_lines(skeletons[cid1].incoming, shortentries, carriageways[cid2]['geometry'].boundary)

    #         # Get the middle line of both skeletons along with the direction of the line
    #         index1, line1, pos1 = get_middle_line(j1, shortentries, skeletons[cid1])
    #         index2, line2, pos2 = get_middle_line(j2, shortentries, skeletons[cid2])
    #         # Create the line between both junctions
    #         linejunction = shapely.LineString([j1, j2])

    #         # If there are incoming lines intersecting the short side entries
    #         if len(incoming) > 0:
    #             # Stores for distances and positions
    #             distances = []
    #             positions = []

    #             # Loop through incoming lines
    #             for i in incoming:
    #                 # Stores geometry
    #                 index = i[0]
    #                 geom = i[1]['geometry']

    #                 # Get start and end point
    #                 start, end = shapely.Point(geom.coords[0]), shapely.Point(geom.coords[-1])
    #                 # Calculate distance between start point and the line between junctions
    #                 startdist = shapely.distance(start, linejunction)
    #                 # Same for end point
    #                 enddist = shapely.distance(end, linejunction)

    #                 # Here, project the start or end point of the incmoming line on the junction line
    #                 # The stored distance represents the distance between the start of the junction line and the projected point on this same line
    #                 if startdist < enddist:
    #                     distances.append(linejunction.project(start))
    #                     positions.append('start')
    #                 else:
    #                     distances.append(linejunction.project(end))
    #                     positions.append('end')

    #             # Calculate the mean of the list of distances
    #             meandist = numpy.mean(distances)
    #             # This value is used to create the new intersection point on the junction line
    #             point = shapely.Point(linejunction.interpolate(meandist).coords)

    #             # Extend both interior lines with the new intersection point
    #             line1 = extend_line_with_point(line1, point, position=pos1)
    #             line2 = extend_line_with_point(line2, point, position=pos2)

    #             # Update each incoming line geometry with a new extended line with the intersection point
    #             for i, inc in enumerate(incoming):
    #                 inc[1]['geometry'] = extend_line_with_point(inc[1]['geometry'], point, position=positions[i])

    #             # Update the geometry of the 'middle' lines in their skeleton
    #             skeletons[cid1].network[index1] = line1
    #             skeletons[cid2].network[index2] = line2

    #         # Here, there is no incoming line, which means both interior lines can be merged along with the junction line
    #         else:
    #             # TODO: handle this situation by updating geometries somewhere...
    #             # This line is one long line which is both skeletons middle lines merged together.
    #             # This situation should not appear.
    #             merged = merge_linestrings(merge_linestrings(line1, linejunction), line2)

    #         # Remove the diamond shaped skeletons part
    #         skeletons[cid1].network = [x for i, x in enumerate(skeletons[cid1].network) if i not in diamond1]
    #         skeletons[cid2].network = [x for i, x in enumerate(skeletons[cid2].network) if i not in diamond2]
    
    # shortdone, pointdone = [], []

    # for skeleton in skeletons:
    #     if skeleton is not None:
    #         # Storing the blended skeleton
    #         collapsed.extend(skeleton.blended)

    # result = []
    # remove = []
    # for i, c in enumerate(collapsed):
    #     cgeom = c['geometry']
    #     add = True
    #     for o in originals:
    #         if shapely.equals(cgeom, network[o]):
    #             remove.append(o)
    #             add = False
    #     if add:
    #         result.append(c)

    # removeroad = []
    # for o in originals:
    #     if o not in remove:
    #         removeroad.append(o)

    # for rid, road in enumerate(roads):
    #     if rid not in removeroad:
    #         result.append(road)

    # return gpd.GeoDataFrame(result, crs=crs)


# import numpy as np
# import geopandas as gpd
# from shapely.geometry import LineString, Point, Polygon, MultiPoint
# from shapely.ops import linemerge, unary_union
# from scipy.spatial import Delaunay
# import warnings

# def simplify_dual_carriageways(road_network, separators, importance_threshold=3):
#     """
#     Simplifie les dual carriageways en créant une route centrale unique.
    
#     Parameters:
#     -----------
#     road_network : GeoDataFrame
#         Le réseau routier avec au minimum une colonne 'geometry' et 'importance'
#     separators : GeoDataFrame
#         Les polygones représentant les séparateurs entre les chaussées
#     importance_threshold : int
#         Seuil d'importance des routes à simplifier (défaut: 3)
    
#     Returns:
#     --------
#     GeoDataFrame
#         Le réseau routier simplifié
#     """
#     def collapse(road_network, separators, importance_threshold):
#         """
#         Effondre les dual carriageways en utilisant une triangulation.
        
#         Returns:
#         --------
#         tuple : (GeoDataFrame des nouvelles routes, liste des indices à supprimer)
#         """
#         collapsed_roads = []
#         roads_to_remove = []
        
#         for idx, separator in separators.iterrows():
#             face_geom = separator.geometry
            
#             # Récupération des routes qui intersectent ce séparateur
#             intersecting_roads = road_network[
#                 road_network.geometry.intersects(face_geom.buffer(0.1))
#             ]
            
#             # Triangulation du contour du séparateur
#             centerline = create_centerline_from_triangulation(face_geom, separators)
            
#             if centerline is None or centerline.is_empty:
#                 continue
            
#             # Trouver l'importance locale et une route modèle pour copier les attributs
#             local_importance = importance_threshold
#             template_road = None
            
#             for road_idx, road in intersecting_roads.iterrows():
#                 if not face_geom.buffer(0.1).intersects(road.geometry):
#                     # Cette route intersecte le séparateur mais n'est pas dedans
#                     local_importance = road.get('importance', importance_threshold)
#                     roads_to_remove.append(road_idx)
#                     if template_road is None:
#                         template_road = road
            
#             # Création de la nouvelle route centrale
#             new_road = create_collapsed_road(
#                 centerline, local_importance, template_road
#             )
#             collapsed_roads.append(new_road)
        
#         # Conversion en GeoDataFrame
#         if collapsed_roads:
#             collapsed_gdf = gpd.GeoDataFrame(
#                 collapsed_roads, 
#                 crs=road_network.crs
#             )
#         else:
#             collapsed_gdf = gpd.GeoDataFrame(columns=road_network.columns, crs=road_network.crs)
        
#         return collapsed_gdf, roads_to_remove


#     # Copie du réseau pour ne pas modifier l'original
#     simplified_network = road_network.copy()
    
#     # Vérifier que les colonnes nécessaires existent
#     if 'importance' not in simplified_network.columns:
#         simplified_network['importance'] = importance_threshold
    
#     # 1. Effondrement des dual carriageways
#     collapsed_roads, roads_to_remove = collapse(
#         simplified_network, separators, importance_threshold
#     )
    
#     # 2. Suppression des anciennes routes
#     simplified_network = simplified_network[~simplified_network.index.isin(roads_to_remove)]
    
#     # 3. Ajout des nouvelles routes simplifiées
#     if len(collapsed_roads) > 0:
#         simplified_network = gpd.GeoDataFrame(
#             pd.concat([simplified_network, collapsed_roads], ignore_index=True),
#             crs=simplified_network.crs
#         )
    
#     # 4. Reconnexion des bretelles (slip roads)
#     simplified_network = reconnect_slip_roads(simplified_network, separators)
    
#     return simplified_network

# def create_centerline_from_triangulation(face_geom, all_separators):
#     """
#     Crée une ligne centrale à partir de la triangulation d'un polygone.
    
#     Parameters:
#     -----------
#     face_geom : Polygon
#         Le polygone du séparateur
#     all_separators : GeoDataFrame
#         Tous les séparateurs (pour détecter les intersections)
    
#     Returns:
#     --------
#     LineString : La ligne centrale créée
#     """
#     # Extraction des points du contour
#     coords = list(face_geom.exterior.coords[:-1])  # Enlever le dernier point (doublon)
    
#     if len(coords) < 3:
#         return None
    
#     points = np.array(coords)
    
#     # Triangulation de Delaunay
#     try:
#         tri = Delaunay(points)
#     except Exception as e:
#         print(f"Erreur de triangulation: {e}")
#         return None
    
#     # Création des segments à partir des triangles
#     segments = []
    
#     for simplex in tri.simplices:
#         triangle_points = points[simplex]
#         triangle_geom = Polygon(triangle_points)
        
#         # Vérifier que le triangle est contenu dans le polygone
#         if not face_geom.contains(triangle_geom.centroid):
#             continue
        
#         # Calculer les segments entre les centres des arêtes
#         segments_triangle = create_segments_from_triangle(
#             triangle_points, face_geom
#         )
#         segments.extend(segments_triangle)
    
#     # Cas spécial: triangulation avec seulement 2 triangles
#     if len(tri.simplices) == 2:
#         # Utiliser le centroïde directement
#         return LineString([face_geom.centroid.coords[0]] * 2)
    
#     # Connexion des segments adjacents
#     if len(segments) == 0:
#         return None
    
#     # Ajouter les segments de connexion avec les faces adjacentes
#     segments = add_connection_segments(face_geom, all_separators, segments)
    
#     # Ordonner les segments pour former une ligne continue
#     centerline = order_segments_to_linestring(segments)
    
#     return centerline


# def create_segments_from_triangle(triangle_points, face_geom):
#     """
#     Crée les segments appropriés à partir d'un triangle de la triangulation.
#     """
#     segments = []
    
#     # Créer les trois arêtes du triangle
#     edges = [
#         LineString([triangle_points[0], triangle_points[1]]),
#         LineString([triangle_points[1], triangle_points[2]]),
#         LineString([triangle_points[2], triangle_points[0]])
#     ]
    
#     # Compter combien d'arêtes sont sur le contour
#     edges_on_contour = sum(
#         1 for edge in edges 
#         if face_geom.exterior.buffer(0.1).contains(edge)
#     )
    
#     # Cas 1: Une arête sur le contour (triangle typique)
#     if edges_on_contour == 1:
#         # Relier les centres des deux autres arêtes
#         edge_centers = [
#             edge.centroid 
#             for edge in edges 
#             if not face_geom.exterior.buffer(0.1).contains(edge)
#         ]
#         if len(edge_centers) == 2:
#             segments.append(LineString([edge_centers[0], edge_centers[1]]))
    
#     # Cas 2: Aucune arête sur le contour (fork/fourche)
#     elif edges_on_contour == 0:
#         # Créer un point central et relier aux centres des arêtes
#         triangle_center = Point(triangle_points.mean(axis=0))
#         for edge in edges:
#             segments.append(LineString([triangle_center, edge.centroid]))
    
#     return segments


# def add_connection_segments(face_geom, all_separators, segments):
#     """
#     Ajoute les segments de connexion avec les faces adjacentes.
#     """
#     # Trouver les faces adjacentes
#     for other_face in all_separators.geometry:
#         if face_geom.equals(other_face):
#             continue
        
#         intersection = face_geom.buffer(0.1).intersection(other_face.buffer(0.1))
        
#         if not intersection.is_empty and intersection.area > 0:
#             # Centroïde de l'intersection
#             connection_point = intersection.centroid
            
#             # Trouver le point le plus proche dans nos segments
#             min_dist = float('inf')
#             closest_point = None
            
#             for segment in segments:
#                 for coord in segment.coords:
#                     pt = Point(coord)
#                     dist = pt.distance(connection_point)
#                     if dist < min_dist:
#                         min_dist = dist
#                         closest_point = pt
            
#             if closest_point and min_dist < 1000:
#                 segments.append(LineString([closest_point, connection_point]))
    
#     return segments


# def order_segments_to_linestring(segments):
#     """
#     Ordonne les segments pour former une LineString continue.
#     """
#     if len(segments) == 0:
#         return None
    
#     if len(segments) == 1:
#         return segments[0]
    
#     # Utiliser linemerge pour connecter les segments
#     try:
#         merged = linemerge(segments)
#         if isinstance(merged, LineString):
#             return merged
#         elif hasattr(merged, 'geoms'):
#             # Prendre la plus longue ligne si plusieurs
#             return max(merged.geoms, key=lambda x: x.length)
#     except:
#         pass
    
#     # Méthode alternative: ordonner manuellement
#     ordered_coords = list(segments[0].coords)
#     remaining = segments[1:]
    
#     max_iterations = 1000
#     iteration = 0
    
#     while remaining and iteration < max_iterations:
#         iteration += 1
#         found = False
        
#         for seg in remaining:
#             seg_coords = list(seg.coords)
            
#             # Vérifier les connexions possibles
#             if Point(ordered_coords[-1]).distance(Point(seg_coords[0])) < 0.01:
#                 ordered_coords.extend(seg_coords[1:])
#                 remaining.remove(seg)
#                 found = True
#                 break
#             elif Point(ordered_coords[-1]).distance(Point(seg_coords[-1])) < 0.01:
#                 ordered_coords.extend(reversed(seg_coords[:-1]))
#                 remaining.remove(seg)
#                 found = True
#                 break
#             elif Point(ordered_coords[0]).distance(Point(seg_coords[-1])) < 0.01:
#                 ordered_coords = list(reversed(seg_coords[:-1])) + ordered_coords
#                 remaining.remove(seg)
#                 found = True
#                 break
#             elif Point(ordered_coords[0]).distance(Point(seg_coords[0])) < 0.01:
#                 ordered_coords = list(reversed(seg_coords))[:-1] + ordered_coords
#                 remaining.remove(seg)
#                 found = True
#                 break
        
#         if not found:
#             break
    
#     if len(ordered_coords) < 2:
#         return None
    
#     return LineString(ordered_coords)


# def create_collapsed_road(centerline, importance, template_road):
#     """
#     Crée une nouvelle route effondrée avec les attributs appropriés.
#     """
#     road_data = {'geometry': centerline, 'importance': importance}
    
#     # Copier les autres attributs de la route modèle si disponible
#     if template_road is not None:
#         for col in template_road.index:
#             if col not in ['geometry', 'importance']:
#                 road_data[col] = template_road[col]
    
#     return road_data


# def reconnect_slip_roads(road_network, separators):
#     """
#     Reconnecte les bretelles (slip roads) sur les nouvelles routes simplifiées.
#     """
#     # Détection des dead-ends (nœuds de degré 1)
#     slip_road_nodes = detect_slip_roads(road_network, separators)
    
#     if len(slip_road_nodes) == 0:
#         return road_network
    
#     # Reconnexion des bretelles
#     reconnected_network = road_network.copy()
    
#     for node_info in slip_road_nodes:
#         node_point = node_info['point']
#         road_idx = node_info['road_idx']
        
#         # Trouver la route centrale la plus proche
#         central_road = find_nearest_central_road(
#             node_point, reconnected_network, separators
#         )
        
#         if central_road is None:
#             continue
        
#         # Projeter le nœud sur la route centrale
#         projected_point = central_road.interpolate(
#             central_road.project(node_point)
#         )
        
#         # Distance maximale pour la reconnexion
#         if node_point.distance(projected_point) > 200:
#             continue
        
#         # Modifier la géométrie de la bretelle
#         road_geom = reconnected_network.loc[road_idx, 'geometry']
#         coords = list(road_geom.coords)
        
#         # Déterminer quel bout prolonger
#         if Point(coords[0]).distance(node_point) < 1:
#             # Prolonger au début
#             new_coords = [projected_point.coords[0]] + coords
#         else:
#             # Prolonger à la fin
#             new_coords = coords + [projected_point.coords[0]]
        
#         reconnected_network.loc[road_idx, 'geometry'] = LineString(new_coords)
    
#     return reconnected_network


# def detect_slip_roads(road_network, separators):
#     """
#     Détecte les bretelles (slip roads) qui se terminent près d'un séparateur.
#     """
#     slip_roads = []
    
#     # Créer un buffer de tous les séparateurs
#     separators_union = unary_union(separators.geometry)
    
#     # Analyser chaque route
#     for idx, road in road_network.iterrows():
#         geom = road.geometry
        
#         # Vérifier les extrémités
#         start_point = Point(geom.coords[0])
#         end_point = Point(geom.coords[-1])
        
#         # Compter les connexions pour chaque extrémité
#         for point in [start_point, end_point]:
#             # Vérifier si c'est un dead-end
#             connections = sum(
#                 1 for _, other_road in road_network.iterrows()
#                 if other_road.geometry.distance(point) < 1 and idx != _
#             )
            
#             # Si dead-end et près d'un séparateur
#             if connections == 0 and point.distance(separators_union) < 10:
#                 slip_roads.append({
#                     'point': point,
#                     'road_idx': idx
#                 })
    
#     return slip_roads


# def find_nearest_central_road(point, road_network, separators):
#     """
#     Trouve la route centrale la plus proche d'un point.
#     """
#     separators_union = unary_union(separators.geometry)
    
#     min_dist = float('inf')
#     nearest_road = None
    
#     for idx, road in road_network.iterrows():
#         # Vérifier si la route est une route centrale (intersecte un séparateur)
#         if not road.geometry.intersects(separators_union.buffer(1)):
#             continue
        
#         dist = road.geometry.distance(point)
#         if dist < min_dist:
#             min_dist = dist
#             nearest_road = road.geometry
    
#     return nearest_road


# # Exemple d'utilisation
# if __name__ == "__main__":
#     import pandas as pd
    
#     # Exemple de données (à remplacer par vos données réelles)
#     # road_network = gpd.read_file("roads.shp")
#     # separators = gpd.read_file("separators.shp")
    
#     # simplified_network = simplify_dual_carriageways(road_network, separators, importance_threshold=3)
#     # simplified_network.to_file("simplified_roads.shp")
    
#     print("Module de simplification des dual carriageways chargé.")
#     print("Utilisez la fonction: simplify_dual_carriageways(road_network, separators, importance_threshold)")