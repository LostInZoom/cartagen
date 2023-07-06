# This file contains algorithms to eliminate buildings in a block
from cartagen4py.utils.multicriteria.promethee import *

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