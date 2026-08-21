import geopandas as gpd
from shapely import MultiPoint, MultiLineString, MultiPolygon

def multi_to_simple(gdf):
    """
    Convert the provided GeoDataFrame to simple geometries.
    """
    def __convert(geometry):
        match geometry:
            case MultiPoint() | MultiLineString() | MultiPolygon():
                l = []
                for simple in geometry.geoms:
                    l.append(simple)
                return l
            case _:
                return [geometry]

    if isinstance(gdf, gpd.GeoDataFrame):
        crs = gdf.crs
        records = gdf.to_dict('records')
        newlist = []
        for r in records:
            newlist.extend(__convert(r))
        return gpd.GeoDataFrame(newlist, crs=crs)
    else:
        raise TypeError('Unsupported data type.')
    

def convert_3d_to_2d(gdf):
    """
    Convert 3D geometries to 2D geometries.
    """
    def __convert(geometry):
        match geometry:
            case _ if hasattr(geometry, 'z'):
                return geometry.buffer(0)
            case _:
                return geometry

    if isinstance(gdf, gpd.GeoDataFrame):
        crs = gdf.crs
        records = gdf.to_dict('records')
        newlist = []
        for r in records:
            newlist.extend(__convert(r))
        return gpd.GeoDataFrame(newlist, crs=crs)
    else:
        raise TypeError('Unsupported data type.')