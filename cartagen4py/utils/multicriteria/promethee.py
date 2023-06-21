from math import exp

class PrometheeCriterion:
    def __init__(self, name, preferenceFunction):
        self.name = name
        self.preferenceFunction = preferenceFunction
    
    def value(self, parameters):
        pass

# Type 1 preference function from Brans & Mareschal 2005, it's a simple function giving 0 for negative 
# (or zero) deviations and 1 for positive deviations.
class Type1PreferenceFunction:
    def __init__(self):
        return
    
    def value(self, deviation):
        if deviation <= 0:
            return 0
        return 1

# Type 2 preference function from Brans & Mareschal 2005, it's a simple function giving 0 for deviations 
# below or equal to a parameter q and 1 for deviations over parameter q.
class Type2PreferenceFunction:
    def __init__(self, q):
        self.q = q
        return
    
    def value(self, deviation):
        if deviation <= self.q:
            return 0
        return 1
    
# Type 3 preference function from Brans & Mareschal 2005, it's a simple function giving 0 for negative (or zero) deviations, 
# linear from 0 to 1 for positive but below parameter p deviations, and giving 1 for deviations over parameter p.
class Type3PreferenceFunction:
    def __init__(self, p):
        self.p = p
        return
    
    def value(self, deviation):
        if deviation <= 0:
            return 0
        elif deviation <= self.p:
            return deviation/self.p
        return 1

# Type 4 preference function from Brans & Mareschal 2005, it's a stair-shaped function giving 0 for deviations below a parameter q, 
# 0.5 for deviations between parameter q and parameter p (p>q) and 1 for deviations over p.
class Type4PreferenceFunction:
    def __init__(self, p, q):
        self.p = p
        self.q = q
        return
    
    def value(self, deviation):
        if deviation <= self.q:
            return 0
        elif deviation <= self.p:
            return 0.5
        return 1

# Type 5 preference function from Brans & Mareschal 2005, it's a function giving 0 for deviations below parameter q, 
# linear between parameter q and parameter p (p>q), and 1 for deviations over p.
class Type5PreferenceFunction:
    def __init__(self, p, q):
        self.p = p
        self.q = q
        return
    
    def value(self, deviation):
        if deviation <= self.q:
            return 0
        elif deviation <= self.p and self.p != self.q:
            return (deviation - self.q) / (self.p - self.q)
        elif deviation <= self.p and self.p == self.q:
            return 0
        return 1

# Type 6 preference function from Brans & Mareschal 2005, it's a gaussian criterion, with a function giving 0 
# for negative (or zero) deviations and tends to 1 when deviation tends to infinity, a parameter s giving the
# inflexion point of the curve.
class Type6PreferenceFunction:
    def __init__(self, s):
        self.s = s
        return
    
    def value(self, deviation):
        if deviation <= 0:
            return 0
        expo = exp(-deviation * deviation / (2 * self.s * self.s))
        return 1-expo

# Picks the best candidate considering the criteria, using the Promethée II complete ranking
# candidates are arrays like [id, criteria_values]
def make_prometheeII_decision(candidates, criteria, weights):
    best = None
    maxNetOutranking = float('-inf')
    # loop on each candidate to compute its positive and negative outranking flow
    for a in candidates:
        positiveOutrankFlow = 0
        negativeOutrankFlow = 0
        for b in candidates:
            if(a[0] == b[0]):
                continue
            # compute the aggregated preference indices of a and b
            aggrPrefIndexAB = 0
            aggrPrefIndexBA = 0
            for i in range(0,len(criteria)):
                criterion = criteria[i]
                weight = weights[i]
                # first compute deviation between a and b
                devAtoB = a[1][i] - b[1][i]
                prefAtoB = criterion.preferenceFunction.value(devAtoB)
                aggrPrefIndexAB += prefAtoB * weight
                prefBtoA = criterion.preferenceFunction.value(-devAtoB)
                aggrPrefIndexBA += prefBtoA * weight
            
            positiveOutrankFlow += aggrPrefIndexAB
            negativeOutrankFlow += aggrPrefIndexBA
        
        positiveOutrankFlow = positiveOutrankFlow / (len(candidates) - 1)
        negativeOutrankFlow = negativeOutrankFlow / (len(candidates) - 1)
        netOutrankingFlow = positiveOutrankFlow - negativeOutrankFlow
        if (netOutrankingFlow > maxNetOutranking):
            best = a
            maxNetOutranking = netOutrankingFlow

    return best


# Rank the candidates considering the criteria, using the Promethée II complete ranking
# candidates are arrays like [id, criteria_values]
def make_prometheeII_ranking_decision(candidates, criteria, weights):
    list_candidates = []

    # loop on each candidate to compute its positive and negative outranking flow
    for a in candidates:
        positiveOutrankFlow = 0
        negativeOutrankFlow = 0
        for b in candidates:
            if(a[0] == b[0]):
                continue
            # compute the aggregated preference indices of a and b
            aggrPrefIndexAB = 0
            aggrPrefIndexBA = 0
            for i in range(0,len(criteria)):
                criterion = criteria[i]
                weight = weights[i]
                # first compute deviation between a and b
                devAtoB = a[1][i] - b[1][i]
                prefAtoB = criterion.preferenceFunction.value(devAtoB)
                aggrPrefIndexAB += prefAtoB * weight
                prefBtoA = criterion.preferenceFunction.value(-devAtoB)
                aggrPrefIndexBA += prefBtoA * weight
            
            positiveOutrankFlow += aggrPrefIndexAB
            negativeOutrankFlow += aggrPrefIndexBA
        
        positiveOutrankFlow = positiveOutrankFlow / (len(candidates) - 1)
        negativeOutrankFlow = negativeOutrankFlow / (len(candidates) - 1)
        netOutrankingFlow = positiveOutrankFlow - negativeOutrankFlow
        list_candidates.append((a,netOutrankingFlow))

    return sorted(list_candidates, key=lambda candidate: candidate[1])