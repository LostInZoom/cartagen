import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import polygonize, unary_union
from shapely.geometry import LineString

class StreetSelector:
    def __init__(self, gdf_routes, 
                 col_stroke_length="stroke_length", 
                 col_importance="stroke_degree",
                 col_centrality="centrality",
                 weights=None):
        """
        Initialise le généralisateur urbain configurable.
        
        Parameters:
        -----------
        gdf_routes : gpd.GeoDataFrame
            Le réseau de routes enrichi.
        col_stroke_length : str
            Nom de la colonne pour la longueur du stroke.
        col_importance : str
            Nom de la colonne pour l'importance/degré (trafic).
        col_centrality : str
            Nom de la colonne pour la centralité pré-calculée.
        weights : dict, optional
            Dictionnaire indiquant l'activation et le poids de chaque critère.
            Ex: {'compact': 1.0, 'stroke': 1.0, 'traffic': 0.0, 'size': 1.0, 'centrality': 1.0}
        """
        self.routes = gdf_routes.copy()
        self.col_stroke_length = col_stroke_length
        self.col_importance = col_importance
        self.col_centrality = col_centrality
        
        # Configuration par défaut des poids (tous activés à 1.0 si None)
        default_weights = {
            'compact': 1.0,
            'stroke': 1.0,
            'traffic': 1.0,
            'size': 1.0,
            'centrality': 1.0
        }
        self.weights = weights if weights is not None else default_weights
        
        # Statistiques globales pour la normalisation des attributs des routes
        self.mean_stroke_length = self.routes[self.col_stroke_length].mean() if self.col_stroke_length in self.routes.columns else 1.0
        self.mean_traffic = self.routes[self.col_importance].mean() if self.col_importance in self.routes.columns else 1.0
        self.mean_centrality = self.routes[self.col_centrality].mean() if self.col_centrality in self.routes.columns else 1.0

    def _calculer_compacite_miller(self, geom):
        """Calcule l'indice de compacité de Miller: 4 * pi * Aire / Perimètre^2"""
        if geom.is_empty or geom.length == 0:
            return 0.0
        return (4 * np.pi * geom.area) / (geom.length ** 2)

    def _construire_blocs(self):
        """Génère les blocs urbains à partir des lignes de routes."""
        lignes = self.routes.geometry.unary_union
        polygones = list(polygonize(lignes))
        
        gdf_blocs = gpd.GeoDataFrame(geometry=polygones, crs=self.routes.crs)
        gdf_blocs["block_id"] = gdf_blocs.index
        return gdf_blocs
    
    def __mettre_a_jour_gdf_blocs(self, gdf_blocs, bloc_a_id, bloc_b_id, taille_min_bloc):
        """
        Fusionne deux blocs spécifiques dans le GeoDataFrame des blocs
        et met à jour les attributs géométriques (aire, compacité).
        """
        # 1. Extraire les deux blocs à fusionner
        bloc_a = gdf_blocs.loc[gdf_blocs["block_id"] == bloc_a_id].iloc[0]
        bloc_b = gdf_blocs.loc[gdf_blocs["block_id"] == bloc_b_id].iloc[0]
        
        # 2. Fusionner les géométries avec l'opération topologique d'union
        nouvelle_geom = unary_union([bloc_a.geometry, bloc_b.geometry])
        
        # Nettoyer d'éventuels artefacts linéaires résiduels de la fusion (ex: SinglePolygon)
        if nouvelle_geom.geom_type == 'MultiPolygon':
            # En principe, ils partagent une frontière donc c'est un Polygon simple.
            # Si c'est un MultiPolygon, on prend le composant principal ou on garde le tout.
            pass

        # 3. Supprimer les deux anciens blocs du GeoDataFrame
        gdf_blocs = gdf_blocs[~gdf_blocs["block_id"].isin([bloc_a_id, bloc_b_id])].copy()
        
        # 4. Créer la ligne pour le nouveau bloc fusionné
        # On conserve l'ID du bloc A (ou on génère un nouvel ID unique)
        nouveau_bloc_id = bloc_a_id 
        
        nouveau_bloc = pd.DataFrame([{
            "geometry": nouvelle_geom,
            "block_id": nouveau_bloc_id,
            "_aire_temp": nouvelle_geom.area
        }])
        
        # Convertir en GeoDataFrame et fusionner avec les blocs restants
        gdf_nouveau = gpd.GeoDataFrame(nouveau_bloc, crs=gdf_blocs.crs)
        gdf_blocs = pd.concat([gdf_blocs, gdf_nouveau], ignore_index=True)
        
        return gdf_blocs

    def _evaluer_cout_aggregation(self, bloc_candidat, bloc_voisin, segment_separateur, taille_min_bloc):
        """
        Calcule le coût d'agrégation pondéré et configurable selon Touya (2010).
        """
        cout_total = 1.0
        
        # 1. Facteur de Compacité (Forme du bloc fusionné)
        if self.weights.get('compact', 0.0) > 0:
            bloc_fusionne_geom = unary_union([bloc_candidat.geometry, bloc_voisin.geometry])
            compacite = self._calculer_compacite_miller(bloc_fusionne_geom)
            factor_compact = (1.0 - compacite) ** 3
            # Application du poids multiplicatif ou exponentiel selon la sensibilité voulue
            cout_total *= (factor_compact * self.weights['compact'])
            
        # 2. Facteur de Taille (Size factor lié à la maille)
        if self.weights.get('size', 0.0) > 0:
            # Plus le bloc candidat est petit par rapport au seuil, plus il est critique à fusionner
            # (Le coût diminue pour privilégier l'élimination des micro-blocs résiduels)
            factor_size = (bloc_candidat.geometry.area + bloc_voisin.geometry.area) / taille_min_bloc
            cout_total *= (factor_size * self.weights['size'])

        # 3. Facteur de Longueur de Stroke
        if self.weights.get('stroke', 0.0) > 0:
            stroke_L = segment_separateur.get(self.col_stroke_length, self.mean_stroke_length)
            factor_stroke = max(0.1, stroke_L / self.mean_stroke_length)
            cout_total *= (factor_stroke * self.weights['stroke'])
            
        # 4. Facteur d'Importance de l'axe (Trafic / Degré)
        if self.weights.get('traffic', 0.0) > 0:
            traffic_val = segment_separateur.get(self.col_importance, self.mean_traffic)
            factor_traffic = max(0.1, traffic_val / self.mean_traffic)
            cout_total *= (factor_traffic * self.weights['traffic'])

        # 5. Facteur de Centralité de la route séparatrice
        if self.weights.get('centrality', 0.0) > 0:
            centrality_val = segment_separateur.get(self.col_centrality, self.mean_centrality)
            # Une forte centralité indique un axe structurel majeur pour le réseau complet: coût élevé
            factor_centrality = max(0.1, centrality_val / self.mean_centrality)
            cout_total *= (factor_centrality * self.weights['centrality'])

        return cout_total

    def generalise(self, taille_min_bloc, cout_max_aggregation):
        """
        Exécute la boucle d'agrégation en mettant à jour le GeoDataFrame 
        des blocs de manière incrémentale.
        """
        # Construction INITIALE unique des blocs
        gdf_blocs = self._construire_blocs()

        blocks_fixed = list()
        while True:
            if gdf_blocs.empty:
                break
                
            # Calcul/Mise à jour des aires
            gdf_blocs["_aire_temp"] = gdf_blocs.geometry.area
            
            # Sélection des blocs trop petits mais qui ne sont pas encore fixés
            petits_blocs = gdf_blocs[(gdf_blocs["_aire_temp"] < taille_min_bloc) & (~gdf_blocs["block_id"].isin(blocks_fixed))]
            if petits_blocs.empty:
                break  # Tous les blocs respectent la taille minimale
                
            # Bloc cible (le plus petit)
            bloc_cible = petits_blocs.nsmallest(1, "_aire_temp").iloc[0]
            frontiere_cible = bloc_cible.geometry.boundary
            voisins_potentiels = gdf_blocs[gdf_blocs.block_id != bloc_cible["block_id"]]
            
            meilleur_voisin_id = None
            cout_minimal = float("inf")
            route_a_supprimer_idx = None
            
            # Évaluation des coûts auprès des voisins
            for _, voisin in voisins_potentiels.iterrows():
                if voisin.geometry.intersects(frontiere_cible):
                    intersection = voisin.geometry.intersection(frontiere_cible)
                    if intersection.length > 1e-3:
                        
                        segments_intersectants = self.routes[self.routes.geometry.intersects(intersection)]
                        if segments_intersectants.empty:
                            continue
                            
                        idx_segment = segments_intersectants.geometry.intersection(intersection).length.idxmax()
                        segment_separateur = self.routes.loc[idx_segment]
                        
                        cout = self._evaluer_cout_aggregation(bloc_cible, voisin, segment_separateur, taille_min_bloc)
                        print(f"Évaluation coût fusion bloc {bloc_cible['block_id']} avec voisin {voisin['block_id']}: {cout:.4f}")
                        if cout < cout_minimal:
                            cout_minimal = cout
                            meilleur_voisin_id = voisin["block_id"]
                            route_a_supprimer_idx = idx_segment

            # Application de la fusion si le coût est valide
            if meilleur_voisin_id is not None and cout_minimal <= cout_max_aggregation:
                # 1. Retirer la route du réseau de sortie
                self.routes = self.routes.drop(index=route_a_supprimer_idx)
                
                # 2. METTRE À JOUR LE GEODATAFRAME DES BLOCS DIRECTEMENT
                gdf_blocs = self.__mettre_a_jour_gdf_blocs(
                    gdf_blocs, 
                    bloc_cible["block_id"], 
                    meilleur_voisin_id, 
                    taille_min_bloc
                )
            else:
                # on ajoute le bloc à la liste des blocs fixes pour éviter de le réévaluer
                blocks_fixed.append(bloc_cible["block_id"])

        return self.routes

def street_selection(gdf_roads, 
                     block_min_size, 
                     aggregation_max_cost,
                     col_stroke_length="stroke_length", 
                     col_importance="stroke_degree",
                     col_centrality="centrality", 
                     weights=None):
    """
    This algorithm proposed by Touya :footcite:p:`touya:2010`, as an extension of the work 
    from A. Ruas, selects and generalizes the urban road network (i.e. the streets in 
    cities) based on configurable criteria. This algorithm is based on the concept of urban blocks 
    and their aggregation to simplify the road network while preserving its essential structure.
    The criteria used to evaluate the cost of merging two blocks are: compactness, stroke length, 
    road importance, block size, and centrality.
    
    Parameters:
    -----------
    gdf_roads : gpd.GeoDataFrame
        A geodataframe containing the road network enriched with additional attributes (stroke length, centrality, importance).
    block_min_size : float
        The minimal size of a block to avoid its aggregation.
    aggregation_max_cost : float
        The maximum acceptable cost for merging two blocks.
    weights : dict, optional
        Dictionary indicating the activation and weight of each criterion.
        
    Returns:
    --------
    gpd.GeoDataFrame
        A geodataframe containing the road network after selection and generalization.
    
    See Also
    --------
        mesh_density_selection : 
            Select road segments based on mesh density.
        strokes_road : 
            Calculate strokes inside a road network.
        rural_traffic :
            Detect central roads inside a network using traffic simulation.
        
    References
    ----------
    .. footbibliography::
    """
    selector = StreetSelector(gdf_roads, 
                              col_stroke_length=col_stroke_length,
                              col_importance=col_importance,
                              col_centrality=col_centrality,
                              weights=weights)
    gdf_roads_generalisees = selector.generalise(block_min_size, aggregation_max_cost)
    return gdf_roads_generalisees


def __compute_mesh_density(gdf_routes):
    """Construit les mailles à partir des lignes de routes et calcule leur densité."""
    # Fusionner les géométries pour créer une topologie propre pour polygonize
    union_lignes = gdf_routes.geometry.unary_union
    polygones = list(polygonize(union_lignes))
    print(f"Nombre de mailles : {len(polygones)}")

    # Si aucune maille fermée n'est formée
    if not polygones:
        return gpd.GeoDataFrame(columns=["geometry", "densite"])

    # Créer le GeoDataFrame des mailles
    gdf_mailles = gpd.GeoDataFrame(geometry=polygones, crs=gdf_routes.crs)

    # Densité = Périmètre / Aire (Formule 6 de l'article)
    gdf_mailles["perimetre"] = gdf_mailles.geometry.length
    gdf_mailles["aire"] = gdf_mailles.geometry.area
    gdf_mailles["densite"] = gdf_mailles["perimetre"] / gdf_mailles["aire"]

    return gdf_mailles


def mesh_density_selection(gdf_roads, seuil_densite, col_importance=None):
    """
    Select the important road segments based on mesh density, following the algorithm by Chen et al. :footcite:p:`chen_selective_2009`
    The algorithm iteratively removes the least important road segments that are located on the boundary of the densest mesh, 
    until all meshes have a density below the specified threshold.
    
    Parameters:
    -----------
    gdf_roads : gpd.GeoDataFrame
        The GeoDataFrame containing the LineString of the road network.
    seuil_densite : float
        The density threshold above which a segment must be eliminated.
    col_importance : str, optional
        The name of the column indicating the importance (e.g., road class).
        If None, the length of the segment is used (the shortest is eliminated).
    
    Returns:
    --------
    gpd.GeoDataFrame
        A geodataframe containing the road network after selection and generalization.

    See Also
    --------
    street_selection : 
        Select urban road segments based on block aggregation.
    strokes_road : 
        Calculate strokes inside a road network.
    rural_traffic :
        Detect central roads inside a network using traffic simulation.
    
    References
    ----------
    .. footbibliography::
    """
    routes_courantes = gdf_roads.copy()

    while True:
        # 1. Recalculer les mailles et leurs densités à l'itération courante
        gdf_mailles = __compute_mesh_density(routes_courantes)

        if gdf_mailles.empty:
            break

        # Trouver la maille avec la densité maximale (Step 1 de l'article)
        maille_max = gdf_mailles.sort_values(by="densite", ascending=False).iloc[0]

        # Si la densité maximale respecte le seuil, l'algorithme a convergé
        if maille_max["densite"] <= seuil_densite:
            break

        polygon_cible = maille_max.geometry
        frontiere = polygon_cible.boundary

        # 2. Identifier les segments routiers situés sur la frontière de cette maille
        indices_frontiere = routes_courantes.index[
            routes_courantes.geometry.intersects(frontiere)
            & (routes_courantes.geometry.intersection(frontiere).length > 1e-3)
        ]

        if len(indices_frontiere) == 0:
            break

        # 3. FILTRAGE : Ne garder que les segments qui séparent deux mailles (segments internes)
        autres_mailles = gdf_mailles.loc[gdf_mailles.index != maille_max.name]
        
        indices_internes = []
        for idx in indices_frontiere:
            segment_geom = routes_courantes.loc[idx].geometry
            intersection_voisins = autres_mailles.geometry.intersects(segment_geom)
            
            if intersection_voisins.any():
                inter_geom = autres_mailles.geometry.intersection(segment_geom)
                if inter_geom.length.sum() > 1e-3:
                    indices_internes.append(idx)

        # Si la maille ne peut pas être fusionnée (uniquement entourée par l'extérieur), on arrête
        if not indices_internes:
            break

        segments_candidats = routes_courantes.loc[indices_internes].copy()

        # 4. TRI MULTI-CRITÈRES (Section 5.1 de l'article Chen 2009)
        # On stocke temporairement la longueur propre du segment routier
        segments_candidats["_longueur_temp"] = segments_candidats.geometry.length

        # Préparation des colonnes et des directions de tri (Ascending)
        colonnes_tri = []
        directions_tri = []

        # Criterium 1: Classe de la route (si fournie)
        if col_importance and col_importance in segments_candidats.columns:
            colonnes_tri.append(col_importance)
            directions_tri.append(True) # Les classes faibles (ex: 1) d'abord

        # Criterium 2: Degré du stroke (Ds)
        if "stroke_degree" in segments_candidats.columns:
            colonnes_tri.append("stroke_degree")
            directions_tri.append(True) # Un petit degré (ex: 1) est moins important

        # Criterium 3: Longueur du stroke (Ls)
        if "stroke_length" in segments_candidats.columns:
            colonnes_tri.append("stroke_length")
            directions_tri.append(False) # Un stroke LONG doit être préservé -> Tri décroissant

        # Criterium 4: Longueur propre du segment (L)
        colonnes_tri.append("_longueur_temp")
        directions_tri.append(True) # Le segment le plus court d'abord

        # Application du tri pour extraire le segment le moins important (index[0])
        segment_a_supprimer = segments_candidats.sort_values(
            by=colonnes_tri,
            ascending=directions_tri
        ).index[0]

        # 5. Élimination du segment
        routes_courantes = routes_courantes.drop(index=segment_a_supprimer)

    return routes_courantes