from shapely.geometry import LineString, Point
import geopandas as gpd
from cartagen.enrichment import StrokeNetwork
import matplotlib.pyplot as plt

data={'geometry':
        [LineString([Point(0, 0),Point(1, 1)]),
         LineString([Point(1, 1),Point(1, 0)]),
         LineString([Point(1, 1),Point(2, 2.2)]),
         LineString([Point(1, 1),Point(2.2, 2)]),
         LineString([Point(2.2, 2),Point(3, 3)]),
         ],
        'name':
            ["rue de la maison blanche",None,None,"rue de la maison blanche","rue de la maison blanche"],
        'id':
            [1,2,3,4,5]}
lines =gpd.GeoDataFrame(data, crs="EPSG:4326")


sn=StrokeNetwork(lines,['name'])

sn.buildStrokes(['name'], 45,30)
array=sn.reconstruct_strokes()
gdf = gpd.GeoDataFrame(array,  columns = ['id', 'geom',"section"],crs="epsg:2154",geometry="geom")   
gdf.plot('id')
plt.show()
