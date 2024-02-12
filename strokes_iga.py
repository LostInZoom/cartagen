import geopandas as gpd
from cartagen4py.data_enrichment import Stroke,StrokeNetwork
import shapely
import geopandas

def main():
    shapefile = gpd.read_file('C://Users//BLe-Mao//workspace_cartagen4py//cartagen4py//cartagen4py//data//data_correction_Iga//nysa2.shp')#je te recommande de le mettre dans data
    sn=StrokeNetwork(shapefile,['idMPHP'])#les attributs sont à passer en chaine de caractère
    sn.buildStrokes(['idMPHP'], 70, 90)#deviatA et deviatSum sont des nombre en degrés
    sn.save_strokes_shp('C://Users//BLe-Mao//Documents//these//Espace_de_travail//data_retravail//projet_preliminaire//strokes_Iga//nysa2.shp')#pour récupérer le SHP des résultats.
if __name__ == "__main__":
    main()