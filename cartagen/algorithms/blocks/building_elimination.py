from cartagen.utils.multicriteria.promethee import (
    PrometheeCriterion,
    Type1PreferenceFunction, Type2PreferenceFunction,
    Type3PreferenceFunction, Type4PreferenceFunction,
    Type5PreferenceFunction, Type6PreferenceFunction,
    make_prometheeII_ranking_decision, make_prometheeII_decision
)
from cartagen.enrichment.urban.building_measures import (
    triangulation_edges_around_building
)

class SizeCriterion(PrometheeCriterion):

    def value(self, parameters):
        polygon = parameters[0]
        simulated_area = parameters[2]
        size_threshold = parameters[3]
        area = polygon.area
        if abs(simulated_area-area) > 15:
            area = simulated_area
        if area <= size_threshold:
            return 1
        if area <= 2*size_threshold:
            a = -0.5 * 1.0 / size_threshold
            b = (2.0 * size_threshold - 1.0) / (2.0 * size_threshold)
            return a * area + b
        a = -1.0 / (4.0 * size_threshold)
        b = 1.0
        return a * area + b

class CornerCriterion(PrometheeCriterion):

    def value(self, parameters):
        building_id = parameters[1]
        corner_buildings = parameters[4]
        if building_id in corner_buildings:
            return 1
        return 0

class CongestionCriterion(PrometheeCriterion):

    def value(self, parameters):
        building_id = parameters[1]
        congestion_data = parameters[5]
        building_congestion = congestion_data[building_id]
        congestion = 0.0
        nb = 0
        for direction in building_congestion[0]:
            if direction > 0.0:
                nb += 1
                congestion += direction
        if nb > 0:
            congestion /= nb

        return congestion
    
class CongestionDirectionCriterion(PrometheeCriterion):

    def value(self, parameters):
        building_id = parameters[1]
        congestion_data = parameters[5]
        building_congestion = congestion_data[building_id]
        nb = 0
        for direction in building_congestion[0]:
            if direction > 0.08:
                nb += 1

        return nb

class NeighboursCriterion(PrometheeCriterion):

    def value(self, parameters):
        polygon = parameters[0]
        building_id = parameters[1]
        congestion_data = parameters[5]
        triangulation = parameters[6]
        building_congestion = congestion_data[building_id]
        edges = triangulation_edges_around_building(polygon, triangulation)
        if len(edges) == 0:
            return 0.0
        
        return building_congestion[1] / len(edges)

def __set_default_criteria():
    criteria = []
    weights = []

    criteria.append(SizeCriterion("size",Type5PreferenceFunction(0.1, 0.7)))
    weights.append(0.2)

    criteria.append(CongestionCriterion("congestion",Type4PreferenceFunction(0.2, 0.4)))
    weights.append(0.2)

    criteria.append(CongestionDirectionCriterion("congestion direction",Type4PreferenceFunction(0.2, 0.4)))
    weights.append(0.2)

    criteria.append(CornerCriterion("corner buildings",Type2PreferenceFunction(0.5)))
    weights.append(0.2)

    criteria.append(NeighboursCriterion("types of neighbours",Type5PreferenceFunction(0.15, 0.6)))
    weights.append(0.2)

    return [criteria, weights]

def building_elimination_best_in_block(buildings, triangulation, congestion, min_area_buildings, corner_buildings, building_elim_thresh=70.0):
    criteria, weights = __set_default_criteria()

    building_id = 0
    candidates = []
    for building in buildings:
        parameters = []
        parameters.append(building)
        parameters.append(building_id)
        parameters.append(max(min_area_buildings,building.area))
        parameters.append(building_elim_thresh)
        parameters.append(corner_buildings)
        parameters.append(congestion)
        parameters.append(triangulation)
        criteria_values = []
        for criterion in criteria:
            criteria_values.append(criterion.value(parameters))
        candidates.append([building_id,criteria_values])

        building_id += 1
    
    return make_prometheeII_decision(candidates, criteria, weights)

def building_elimination_ranking_in_block(buildings, triangulation, congestion, min_area_buildings, corner_buildings, building_elim_thresh=70.0):
    criteria, weights = __set_default_criteria()

    building_id = 0
    candidates = []
    for building in buildings:
        parameters = []
        parameters.append(building)
        parameters.append(building_id)
        parameters.append(max(min_area_buildings,building.area))
        parameters.append(building_elim_thresh)
        parameters.append(corner_buildings)
        parameters.append(congestion)
        parameters.append(triangulation)
        criteria_values = []
        for criterion in criteria:
            criteria_values.append(criterion.value(parameters))
        candidates.append([building_id,criteria_values])

        building_id += 1
    
    return make_prometheeII_ranking_decision(candidates, criteria, weights)