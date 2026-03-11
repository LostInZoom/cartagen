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

    remove = []
    crossroads = []
    for face in faces:
        # Create the crossroad object without the removed roads
        crossroad = Crossroad(network, roads.sindex, face)
        crossroads.append(crossroad)

        # If there are more than one external road to the crossroad
        if len(crossroad.externals) > 1:
            # Stores internal roads for removal
            remove.extend([ crossroad.original[x] for x in crossroad.internals])

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

            # TODO: Include spinalization method
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

    return cleaned