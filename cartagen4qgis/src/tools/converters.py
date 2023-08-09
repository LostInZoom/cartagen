import geopandas as gpd
from shapely.wkt import loads
from qgis.core import QgsFeature, QgsGeometry

def qgis_source_to_geodataframe(source):
    """
    Converts a QGIS source to a geopandas GeoDataFrame.
    """
    crs = source.sourceCrs().authid()

    features = source.getFeatures()
    f = []
    for feature in features:
        entity = feature.__geo_interface__["properties"]
        entity['geometry'] = loads(feature.geometry().asWkt())
        f.append(entity)

    return gpd.GeoDataFrame(f, crs=crs)

def list_to_qgis_feature(dicts):
    """
    Converts a list of dicts with attributes and geometry to a list of qgis features.
    """
    # gdf = gpd.GeoDataFrame(dicts, crs)
    # return QgsVectorLayer(gdf.to_json())

    features = []

    for d in dicts:
        feature = QgsFeature()
        feature.setGeometry(QgsGeometry.fromWkt(d['geometry'].wkt))
        features.append(feature)

    return features