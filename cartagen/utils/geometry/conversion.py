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