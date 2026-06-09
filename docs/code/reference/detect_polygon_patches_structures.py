from matplotlib import pyplot as plt
from matplotlib.path import Path

import geopandas as gpd
from shapely.geometry import Polygon
import cartagen as c4

def generate_test_islands():
    """
    Génère un GeoDataFrame de polygones simples (carrés/hexagones) 
    simulant des structures d'îles (une ligne d'îles, une grappe et des îles isolées).
    """
    polygons = []
    
    # --- GROUPE 1 : Un alignement d'îles (structure linéaire / chaîne d'îles) ---
    # 5 îles proches les unes des autres le long de l'axe X
    for i in range(5):
        x_center = i * 8  # Espacées de 8 unités
        y_center = 10
        # Création d'un petit carré de côté 6 (rayon ~3)
        poly = Polygon([
            (x_center - 3, y_center - 3),
            (x_center + 3, y_center - 3),
            (x_center + 3, y_center + 3),
            (x_center - 3, y_center + 3)
        ])
        polygons.append(poly)

    # --- GROUPE 2 : Une grappe dense d'îles (structure en cluster) ---
    # Un nœud central (graine potentielle) entouré de 4 îles très proches
    cluster_center_x = 30
    cluster_center_y = 50
    
    # L'île centrale du cluster (un peu plus grande)
    polygons.append(Polygon([
        (cluster_center_x - 5, cluster_center_y - 5),
        (cluster_center_x + 5, cluster_center_y - 5),
        (cluster_center_x + 5, cluster_center_y + 5),
        (cluster_center_x - 5, cluster_center_y + 5)
    ]))
    
    # 4 îles satellites autour du centre
    satellites_offsets = [(-12, 0), (12, 0), (0, -12), (0, 12)]
    for dx, dy in satellites_offsets:
        cx = cluster_center_x + dx
        cy = cluster_center_y + dy
        poly = Polygon([
            (cx - 3, cy - 3),
            (cx + 3, cy - 3),
            (cx + 3, cy + 3),
            (cx - 3, cy + 3)
        ])
        polygons.append(poly)

    # --- HORS GROUPE : Îles isolées ---
    # Île isolée 1 : Très loin en haut à droite
    polygons.append(Polygon([(80, 80), (86, 80), (86, 86), (80, 86)]))
    
    # Île isolée 2 : Entre la chaîne et le cluster mais trop éloignée pour être connectée
    polygons.append(Polygon([(65, 30), (69, 30), (69, 34), (65, 34)]))

    # Création du GeoDataFrame
    gdf = gpd.GeoDataFrame(geometry=polygons)
    return gdf

# 1. Générer le jeu de données
gdf_test = generate_test_islands()
    
# 2. Exécuter l'algorithme fourni précédemment 
# (Assurez-vous que la fonction `detect_meso_structures` est définie ou importée)
gdf_results = c4.detect_polygon_patches_structures(gdf_test, k=3, radius_method='area')
        
 # 3. Affichage des résultats
fig, ax = plt.subplots(figsize=(8, 8))
# Afficher les îles colorées par leur ID de groupe (-1 = isolée)
gdf_results.plot(column='group_id', cmap='Set1', categorical=True, legend=True, ax=ax, edgecolor='black')
        
# Ajouter les index pour faciliter la lecture
for idx, row in gdf_results.iterrows():
    ax.annotate(text=str(idx), xy=(row.geometry.centroid.x, row.geometry.centroid.y), color='white', weight='bold', ha='center', va='center')
            
plt.show()