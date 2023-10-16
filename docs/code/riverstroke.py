from shapely.geometry import LineString, Point
import geopandas as gpd
from cartagen4py.data_enrichment import RiverStrokeNetwork

data={'geometry':
        [LineString([Point(1,4),Point(1, 3)]),
         LineString([Point(1.5,3.5),Point(1, 3)]),
         LineString([Point(1, 3),Point(1, 2.4)]),
         LineString([Point(1, 2.4),Point(0.8, 1.8),Point(0.9, 1.5)]),
         LineString([Point(1, 2.4),Point(1.2, 2.1)]),
         LineString([Point(1.2, 2.1),Point(0.9, 1.5)]),
         
         LineString([Point(0.9, 1.5),Point( 1.2,0.6)]),
         
         LineString([Point(1.2, 2.1),Point( 1.2,0.6)]),
         LineString([Point( 1.2,0.6),Point(1.1, 0.3)]),
         LineString([Point(1.1, 0.3),Point(1, 0)]),
         LineString([Point(0.5, 2),Point(1.1, 0.3)]),
         ],

        'id':
            [1,2,3,4,5,6,8,9,10,11,12]}
lines =gpd.GeoDataFrame(data, crs="EPSG:4326")

sn=RiverStrokeNetwork(lines,None)

sn.buildRiverStrokes([], 45,30)
array=sn.reconstruct_strokes()
gdf = gpd.GeoDataFrame(array,  columns = ['id', 'geom',"strahler"],crs="epsg:4326",geometry="geom")   
gdf.plot('id')
gdf.plot('strahler')