import math
from itertools import combinations
from shapely.geometry import Point, MultiLineString
from shapely import ops
import geopandas as gpd
import networkx as nx

from cartagen.utils.geometry.angle import angle_3_pts, angle_between_2lines

def strokes_roads(roads, attributes, angle=45.0, angle_sum=30.0):
    """
    Calculate strokes inside a road network.

    This algorithm is based on the 'Good Continuation' principle
    as defined by Thomson & Richardson. :footcite:p:`thomson:2002`
    It defines strokes, which are segments that follow this perceptual
    grouping also known as Gestalt's principle.

    Parameters
    ----------
    roads : GeoDataFrame of LineString
        The roads to calculate the strokes from.
    attributes : list of str
        The list of attributes to help the derivation of continuity.
    angle : float, optional
        Thresholds for the maximum angle between two segments
        at the junction of two sections belonging to the same stroke.
    angle_sum : float, optional
        Thresholds for the maximum angle between
        two sections belonging to the same stroke.

    Returns
    -------
    GeoDataFrame of LineString

    References
    ----------
    .. footbibliography::
    """
    crs = roads.crs
    strokes = StrokeNetwork(roads, attributes)
    strokes.buildStrokes(attributes, angle, angle_sum)
    s = strokes.reconstruct_strokes()
    return gpd.GeoDataFrame(s,  columns = ['id','geometry','section'], crs=crs, geometry='geometry') 

def strokes_rivers(rivers, attributes, angle=45.0, angle_sum=30.0):
    """
    Calculate strokes inside a river network.

    Strokes are network segments that follow the perceptual
    grouping principle of Good Continuity (Gestalt).
    """
    crs = rivers.crs
    strokes = RiverStrokeNetwork(rivers, attributes)
    strokes.buildRiverStrokes(attributes, angle, angle_sum)
    s = strokes.reconstruct_strokes()
    return gpd.GeoDataFrame(s,  columns = ['id','geometry','strahler'], crs=crs, geometry='geometry') 

class Stroke:
    """
    
    """
    COUNTER = 0
    def __init__(self, network, root):
        self.network = network
        self.root=root
        self.features = []
        self.features.append(root)
        self.id = Stroke.COUNTER
        Stroke.COUNTER += 1
    
    def is_good_continuity(temp_geom, geom_foll, deviatAngle, deviatSum):
        angleThresh = deviatAngle / 180.0 * math.pi
        sumThresh = deviatAngle / 180.0 * math.pi
    
        coord_ini1 = temp_geom["geom"].coords[0]
        coord_fin1 = temp_geom["geom"].coords[-1]
        coord_ini2 = geom_foll["geom"].coords[0]
        coord_fin2 = geom_foll["geom"].coords[-1]
        coord_inter = None
        
        if coord_ini2==coord_ini1:
            coord_inter = coord_ini1
            inter_geom1 = True
            inter_geom2 = True
        if coord_fin2==coord_ini1:
            coord_inter = coord_ini1
            inter_geom1 = True
            inter_geom2 = False
        if coord_fin2==coord_fin1:
            coord_inter = coord_fin1
            inter_geom1 = False
            inter_geom2 = False
        if coord_ini2==coord_fin1:
            coord_inter = coord_ini2
            inter_geom1 = False
            inter_geom2 = True
        if coord_inter is None:
            return -1.0
        
        v1g1, v2g1, v1g2, v2g2 = None, None, None, None
        # define the node to compute the angle
        # get the first vertex on geometry 1
        if inter_geom1:
            v1g1 = Point(temp_geom["geom"].coords[1])
        else:
            v1g1 = Point(temp_geom["geom"].coords[-2])
            
        # get the first vertex on geometry 2
        if inter_geom2:
            v1g2 = Point(geom_foll["geom"].coords[1])
        else:
            v1g2 = Point(geom_foll["geom"].coords[-2])
            
        # if nbVert1 > 2, get the second vertex in geometry 1
        if len(temp_geom["geom"].coords) > 2:
            if inter_geom1:
                v2g1 = Point(temp_geom["geom"].coords[2])
            else:
                v2g1 = Point(temp_geom["geom"].coords[-3])
    
        # if nbVert2 > 2, get the second vertex in geometry 2
        if len(geom_foll["geom"].coords) > 2:
            if inter_geom2:
                v2g2 = Point(geom_foll["geom"].coords[2])
            else:
                v2g2 = Point(geom_foll["geom"].coords[-3])
    
        # now, compute interAngle between geom and geomFoll
        inter_angle = angle_3_pts(v1g1,  Point(coord_inter), v1g2)
        # put the angle between -pi and pi
        if inter_angle > math.pi:
            inter_angle -= 2 * math.pi    
        #case where both geometries have only 2 vertices
        if ((v2g1 is None) and (v2g2 is None)) :
          #then, there is good continuity if the angle is < 45°
          if ((inter_angle < (-angleThresh)) or (inter_angle > angleThresh)) :
            return math.pi - abs(inter_angle)
          return -1.0
        #case where geom has 2 vertices
        elif (v2g1 is None):
          angleTotalDiff = 0.0
          #on calcule angleGeom2
          angleGeom2 = angle_3_pts( Point(coord_inter), v1g2, v2g2)
         #on calcule l'écart entre les angles
          angleDiff = max(angleGeom2, inter_angle)- min(angleGeom2, inter_angle)
          if (angleDiff > math.pi) :
            angleTotalDiff = abs(angleDiff - 2 * math.pi)
          else:
            angleTotalDiff = angleDiff
          #il y a bonne continuité si l'angle est < 45° et la différence desangles < à 30°
          if (((inter_angle < (-angleThresh)) or (inter_angle > angleThresh)) and (angleTotalDiff < sumThresh)):
            return 2 * angleTotalDiff
          return -1.0
        #case where geomFoll has 2 vertices
        elif (v2g2 is None):
          angleTotalDiff = 0.0
          #on calcule angleGeom2
          angleGeom1 = angle_3_pts(v2g1, v1g1,  Point(coord_inter))
          #on calcule l'écart entre les angles
          angleDiff = max(angleGeom1, inter_angle)- min(angleGeom1, inter_angle)
          if (angleDiff > math.pi) :
            angleTotalDiff = abs(angleDiff - 2 * math.pi)
          else:
            angleTotalDiff = angleDiff
          #il y a bonne continuité si l'angle est < 45° et la différence des angles < à 30°
          if (((inter_angle < (-angleThresh)) or (inter_angle > angleThresh)) and (angleTotalDiff < sumThresh)) :
            return 2 * angleTotalDiff
          return -1.0
        #general case
        else:
            angleTotalDiff1 = angleTotalDiff2 = 0.0
            angleGeom1 = angle_3_pts(v2g1, v1g1,  Point(coord_inter))
            angleGeom2 = angle_3_pts( Point(coord_inter), v1g2, v2g2)
            #on calcule l'écart entre les angles 1 et inter
            angleDiff1 = max(angleGeom1, inter_angle)-min(angleGeom1, inter_angle)
            if (angleDiff1 > math.pi) :
                angleTotalDiff1 = abs(angleDiff1 - 2 * math.pi)
            else :
                angleTotalDiff1 = angleDiff1
            #on calcule l'écart entre les angles
            DiffAngles2 = max(angleGeom2, inter_angle)-min(angleGeom2, inter_angle)
            if (DiffAngles2 > math.pi):
                angleTotalDiff2 = abs(DiffAngles2 - 2 * math.pi)
            else :
                angleTotalDiff2 = DiffAngles2
            #il y a bonne continuité si l'angle est < 45° et les différences des angles < à 30°
            if (((inter_angle < (-angleThresh)) or (inter_angle > angleThresh)) and (angleTotalDiff1 < sumThresh) and (angleTotalDiff2 < sumThresh)):
                return angleTotalDiff1 + angleTotalDiff2
            return -1.0
    
    def attribute_filter(arc,followers,attributeNames):
        loopFoll = []
        for attribute in attributeNames:
            # get the attribute value for 'arc'
            value = arc[attribute]
            if len(followers) != 0:
                # loop on the followers to filter them
                loopFoll = []
                loopFoll+=followers
                for a in loopFoll:
                    # get the value of a for attribute
                    valueA = a[attribute]
                    if value != valueA:
                        # remove 'a' from the followers set
                        followers.remove(a)

    def filterFollowers(self, arc, followers):
        loopFoll = []
        #loop on the followers to filter them
        loopFoll+=followers
        for a in loopFoll:
            if (a==arc):
                #remove it from the set
                followers.remove(a)
                continue
             #check that a does not belong to another stroke
            if (a in self.network.groupedFeatures):
                  #remove it from the set
                  followers.remove(a)
                  continue
                  #check if it belongs to this stroke
            if (a in self.features):
             #remove it from the set
                 followers.remove(a)
                 continue
            #check if it belongs to the network
            if (not a in self.network.features) :
            # remove it from the set
                followers.remove(a)
      
                        
    def chooseNextSegment(self,arc, followers, attributeNames, deviatAngle, deviatSum):
        #first, if it's a node of degree two
        if (len(followers) == 1) :
            follower = next(iter(followers))
            if (not follower in self.features):#COR
                return follower
            return None
        # then, filter the followers
        self.filterFollowers(arc, followers)
        if (len(followers) == 0) :
            return None
        # then, filter the followers from the attributeNames
        Stroke.attribute_filter(arc, followers, attributeNames)
        if (len(followers) == 0) :
            return None
        continuity,bestSegment = True, None
        # Loop on the followers to choose the best continuity
        minDiff = math.pi
        for follower in followers:
            # get the continuity difference with this follower
            diffContinuity = Stroke.is_good_continuity(arc, follower, deviatAngle, deviatSum) 
            if (diffContinuity > -1.0) :
                continuity = True
                if (diffContinuity < minDiff) :
                    # this is the current best continuity update the difference
                    minDiff = diffContinuity
                    # change the bestSegment
                    bestSegment = follower
        
        # final verification: if we found a continuous follower, check if there is a better continuity between the followers themselves
        if (((continuity) and not bestSegment in self.network.groupedFeatures) or (minDiff < deviatAngle)) :
            for pair in combinations(followers, 2):
                diffContinuity = Stroke.is_good_continuity(pair[0], pair[1], deviatAngle, deviatSum)
                if (diffContinuity > -1.0 and diffContinuity < minDiff) :
                    # there is a pair of followers with a better continuity so we stop the stroke here and return None
                    return None
            return bestSegment
        return None

  
    def one_side_stroke(self, side, attributeNames, deviatAngle, deviatSum):
        # get the following network segments of the root of this stroke
        node =None        
        node = self.root["geom"].coords[side]
        if (node is None):
            return
        followers =[]
        followers+=self.network.dic_neighbours[self.root["geom"].coords[side]]
        if self.root in followers:
            followers.remove(self.root)
        next1 = self.root
        continuity = True
        self.network.groupedFeatures.append(next1)
        while (continuity) :
            # get the best candidate among the followers (the one with best continuity)
            best = self.chooseNextSegment(next1, followers, attributeNames,deviatAngle, deviatSum)
            if (best is None):
                break
            # add this to the 2 sets (the network one and the stroke one)
            if not side:
                self.features.insert(0, best)
            else:
                self.features.append(best)
            self.network.groupedFeatures.append(best) 
            # get the followers of 'best'
            followers=[]
            nextNode = best["geom"].coords[0]
            
            side2=0
            if (node==nextNode):
                nextNode = best["geom"].coords[-1]
                side2=-1 
                
            followers+=self.network.dic_neighbours[best["geom"].coords[side2]]#AC
            followers.remove(best)
            # if there is no follower, break
            if (len(followers)== 0):
                break
            # update the 'next' segment with 'best'
            next1 = best
            node = nextNode;         
    def __str__(self):
        liste=""
        for elem in self.features: 
            liste+=str(elem["id"])
            liste+=","
        return liste
    
class StrokeNetwork:
    """
    This Class contains methods allowing the computation of strokes
    in a line network representing geographic entities (e.g., roads).

    The initialization of this class is required prior to computing strokes,
    it includes the precomputing of neighbouring relations between edges of the network.
    """
    def __init__(self, features):
            #Initialisation from a list of shapely geometry with the correct attributes
          self.features = features
          self.groupedFeatures = []
          self.id = 0
          self.strokes = []
          
    def __init__(self, gdf, attributeNames):
        #Initialisation from a geopanda dataframe and the liste of desired attribute name
        features=[]

        entries = gdf.to_dict('records')

        for i, entry in enumerate(entries): 
            elem = {}
            elem["geom"] = entry['geometry']
            elem["id"] = i

            for attr in attributeNames:
                if attr in entry:
                    elem[attr] = entry[attr]
                else:
                    raise Exception("No attribute named '{0}' in the provided dataset.".format(attr))
                
            features.append(elem)

        self.features = features
        self.groupedFeatures = []
        self.id = 0
        self.strokes = []
        self.dic_neighbours=StrokeNetwork.compute_neighbours(features)
        
    def compute_neighbours(network):
        ntx=nx.DiGraph()
        dic_neighbours={}
        for edge in network:
            ntx.add_edge(edge["geom"].coords[0],edge["geom"].coords[-1],elem=edge)
            ntx.add_edge(edge["geom"].coords[-1],edge["geom"].coords[0],elem=edge)
        for node in ntx.nodes:
            dic_neighbours[node]=[]
            for succnode in ntx.successors(node):
                if ntx.has_edge(node,succnode) :
                    temp=ntx[node][succnode]["elem"]
                    dic_neighbours[node].append(temp)
                if ntx.has_edge(succnode,node) : 
                    temp=ntx[node][succnode]["elem"]
                    dic_neighbours[node].append(temp)
        return dic_neighbours


    def buildStrokes(self, attributeNames, deviatAngle, deviatSum) :
        """
        This method computes the strokes in a Strokenetwork
        using a loop on network features, and updates its strokes attribute.
        """
        #loop on the network features
        for obj in self.features :
            #test if the feature has already been treated
            if obj in self.groupedFeatures :
                continue
            #build a new stroke object
            stroke = Stroke(self, obj)
            #// build the stroke on the initial side
            stroke.one_side_stroke(0, attributeNames, deviatAngle, deviatSum)
            #// build the stroke on the final side
            stroke.one_side_stroke(-1, attributeNames, deviatAngle, deviatSum)
            #// add the stroke to the strokes set
            self.strokes.append(stroke)
    
    def reconstruct_strokes(self):
        strokes=self.strokes
        array = []
        for i, stroke in enumerate(strokes):
            listline=[]
            section=""
            for j, seg in enumerate(stroke.features):
                listline+=[seg["geom"]]
                section += str(seg["id"])
                if j < len(stroke.features):
                    section += ","
            multi_line = MultiLineString(listline)
            merged_line = ops.linemerge(multi_line)
            array += [[i,merged_line,section]]
        return array

    def save_strokes_shp(self,path):
        #Save a shapefile of stroke eometry in the desired folder
        strokes_list = self.strokes
        array = self.reconstruct_strokes(strokes_list)
        gdf = gpd.GeoDataFrame(array,  columns = ['id', 'geom',"section"],crs="epsg:2154",geometry="geom")   
        gdf.to_file(path, driver='ESRI Shapefile')

    def get_strokes(self):
        return self.reconstruct_strokes(self.strokes)

class RiverStroke:
    def __init__(self, network, section): 
        if type(section) is list:
            self.features = section 
        else:
            self.features = []
            self.features.append(section)
        self.network = network
        self.isBraided = None

    def setBraided(self,isBraided) :
        self.isBraided = isBraided;
        
    def getlength(self) :
        l=0
        for feat in self.features:
            l+=feat['geom'].length
        return l

class RiverStrokeNetwork:
    def __init__(self, lines, attributeName):
        #Initialisation from a geopanda dataframe and the liste of desired attribute name
        features=[]

        if attributeName is None: 
            self.NoAttribute=True
        else:
            self.NoAttribute=False
            self.attributeName=attributeName
        for idx in lines.index: 
            elem={}
            elem["geom"]=lines.geometry[idx]
            elem["id"]=idx
            if not self.NoAttribute:
                elem["name"]=getattr(lines,attributeName)[idx]
            features+=[elem]

        self.features = features
        self.sources = []
        self.sinks = []
        self.strahlerOrders={}  
        self.groupedFeatures = []
        self.strokes = []
        self.dic_neighbours,self.dicentrant=self.compute_neighbours(features)
        self.dtf=lines
        
    def compute_neighbours(self, network):
        ntx=nx.DiGraph()
        dic_neighbours={}
        dicentrant={}
        for edge in network:
            ntx.add_edge(edge["geom"].coords[0],edge["geom"].coords[-1],elem=edge)
        for node in ntx.nodes:
            dicentrant[node]=[]
            dic_neighbours[node]=[]
            la= list(ntx.successors(node))
            lb= list(ntx.predecessors(node))
            if (not len(la)==0) and  len(lb)==0: 
                self.sources.append(node)
            if (not len(lb)==0) and  len(la)==0: 
                self.sinks.append(node)    
            for prednode in ntx.predecessors(node):
                if ntx.has_edge(prednode,node) : 
                    temp=ntx[prednode][node]["elem"]
                    dicentrant[node].append(temp)
            for succnode in ntx.successors(node):
                if ntx.has_edge(node,succnode) :
                    temp=ntx[node][succnode]["elem"]
                    dic_neighbours[node].append(temp)  
        return dic_neighbours,dicentrant


    def allBelongToStroke(self, sections) :
        for section in sections:
            if (not section in self.groupedFeatures) :
                return False
        return True
    
    def getNonBraidedStroke(self, upstreamStrokes) :
        nb = 0
        unbraided = None
        for stroke in upstreamStrokes :
            if (stroke.isBraided) :
                continue
            unbraided = stroke;
            nb+=1
        if (nb == 1) :
            return unbraided
        return None
    
    def makeDecisionAtConfluence(self, node, downstreamSection, upstreamStrokes) :
        #if there is only one upstream river, it continues
        if (len(upstreamStrokes) == 1) :
            return upstreamStrokes[0]
        #first, make a decision on braided strokes
        unbraidedStroke = self.getNonBraidedStroke(upstreamStrokes)
        if (unbraidedStroke is not None): 
            return unbraidedStroke
        #first, make a decision on river name
        if not self.NoAttribute: 
            if not downstreamSection["name"] is None:        
                for stroke in upstreamStrokes:
                    if (downstreamSection["name"]==stroke.features[-1]["name"]):
                        return stroke
        #arrived here, decision is made on stroke length
        longest = None
        second = None
        maxLength = 0.0
        secondLength = 0.0
        for stroke in upstreamStrokes :
            length = stroke.getlength()
            if (length > maxLength):
                second = longest
                secondLength = maxLength
                longest = stroke
                maxLength = length
            elif (length > secondLength) :
                second = stroke
                secondLength = length
        if (longest is not None and  second is not None):#AC a verifier ici
            if (longest.getlength() / second.getlength() >= 3) :
                return longest
        # arrived here, make difference on angle difference
        best =None
        mini = 1.0
        for stroke in upstreamStrokes:
            angle = angle_between_2lines(stroke.features[-1]['geom'],downstreamSection['geom'])
            if (math.cos(angle) < mini) :
                mini = math.cos(angle)
                best = stroke
        return best;


    def getUpstreamStrokes(self, node) :
        upstreamStrokes = []
        pfeats=[]
        for stroke in self.strokes :
            feats = stroke.features
            pfeats=[]
            for feat in feats : 
                if feat in self.dicentrant[node]:
                    pfeats.append(feat)
            if (not len(pfeats) == 0):
                upstreamStrokes.append(stroke)
        return upstreamStrokes

    def computeStrahlerAtConfluence(self, orders) :
        if (len(orders) == 1):
            return orders[0]
        if (len(orders) == 0):
            return 1
        maxi = max(orders);
        nbMax =orders.count(maxi)
        if (nbMax > 1) :
            return maxi + 1
        return maxi
    
    def manageBraidedConfluence(self, node, upStream, downstreamNodes) :
        mainBranch = None
        mini = 1.0
        for branch in self.dic_neighbours[node]:
			# first, make a decision on river name
            if not self.NoAttribute: 
                if branch['name']==upStream["name"]:
                    mainBranch = branch
                    break
			#then decide on angles
            angle = angle_between_2lines(upStream["geom"], branch["geom"])
            if (math.cos(angle) < mini):
                mini = math.cos(angle)
                mainBranch = branch
        if mainBranch is None :
            return None
		#continue the upstream stroke with mainBranch
        upstreamStroke = self.getUpstreamStrokes(node)[0]
        upstreamStroke.features.append(mainBranch)
        self.groupedFeatures.append(mainBranch)
		#find the next node
        if (not mainBranch["geom"].coords[-1] in downstreamNodes):
            downstreamNodes.insert(0, mainBranch["geom"].coords[-1])
        self.strahlerOrders[mainBranch['id']]=self.strahlerOrders[upStream['id']]
        return mainBranch

    def getMainNonBraidedStroke(self, upstreamStrokes) :
        maxLength = 0
        unbraided = None
        for stroke in upstreamStrokes :
            if (stroke.isBraided):
                continue
            l=stroke.getlength()
            if (l > maxLength):
                unbraided = stroke
                maxLength = l
        return unbraided

    def buildRiverStrokes(self, attributeNames, deviatAngle, deviatSum) :
        downstreamNodes = []
        total=len(self.features)
        for source in self.sources:
			#first get the downstream section
            section = self.dic_neighbours[source][0]
			#build a new RiverStroke with this section
            stroke = RiverStroke(self, section)
            self.strokes.append(stroke)
            self.groupedFeatures.append(section)
            if (not  section['geom'].coords[-1] in downstreamNodes) :
                downstreamNodes.insert(0,section['geom'].coords[-1])
			# compute Strahler orders
            self.strahlerOrders[section['id']] = 1
		# ****************************
		# MAIN LOOP
		# ****************************
        counter = 0;
        preva=int(len(self.groupedFeatures)*100/total)
        
        while (len(downstreamNodes)>0) :
            
			#test to break the loop when it's stuck
            if (counter == len(downstreamNodes)) :
                break
            node = downstreamNodes.pop()
            if (node in self.sinks): 
                continue
			#check that all entering sections already belong to a stroke
            if (not self.allBelongToStroke(self.dicentrant[node])) :
                downstreamNodes.insert(0, node)
                counter+=1
                continue

			#arrived there, it has to be decided which stroke is stopped and which one continues.
            counter = 0
            if (len(self.dic_neighbours[node]) == 1) :
                downstreamSection = self.dic_neighbours[node][0]
                if downstreamSection in self.groupedFeatures:
                    continue
                #get the upstream strokes and find the one that continues
                continuing = self.makeDecisionAtConfluence(node, downstreamSection,self.getUpstreamStrokes(node))#ac[]
    
                #now extends the continuing stroke with downstreamSection
                continuing.features.append(downstreamSection)#AC features
                self.groupedFeatures.append(downstreamSection)
                
                # find the next node
                if (not downstreamSection['geom'].coords[-1] in downstreamNodes):
                    
                    downstreamNodes.insert(0,downstreamSection['geom'].coords[-1])
                #compute Strahler order
                orders =[]
                for arc in self.dicentrant[node] :
                    orders.append(self.strahlerOrders[arc['id']])
                if self.dicentrant[node]==[]:
                    self.strahlerOrders[downstreamSection['id']] = 1
                else:
                    self.strahlerOrders[downstreamSection['id']] = self.computeStrahlerAtConfluence(orders)
                
                
            elif (len(self.dic_neighbours[node]) > 1) :
                #braided stream case
                if (len(self.dicentrant[node]) == 1):
                    #normal case the main branch has to be found
                    upStream =  self.dicentrant[node][0]
                    mainBranch = self.manageBraidedConfluence(node, upStream, downstreamNodes)
                    #build new RiverStrokes with remaining branches
                    remainingBranches = self.dic_neighbours[node]
                    remainingBranches.remove(mainBranch)
                    for branch in remainingBranches :
                        stroke = RiverStroke(self, branch)
                        stroke.setBraided(True)
                        self.strokes.append(stroke)
                        self.groupedFeatures.append(branch)
                        if (not branch['geom'].coords[-1] in downstreamNodes) :
                            downstreamNodes.insert(0, branch['geom'].coords[-1])
                        # compute Strahler orders
                        self.strahlerOrders[branch['id']] = 1
                else :
                    #complex case with braids at a confluence point
                    upstreamStrokes = self.getUpstreamStrokes(node)
                    remainingBranches = self.dic_neighbours[node]
                    unbraided = self.getMainNonBraidedStroke(upstreamStrokes)
                    if (unbraided is not None) :
                        upstreamStrokes.remove(unbraided)
                        mainBranch = self.manageBraidedConfluence(node,unbraided.features[-1], downstreamNodes)
                        remainingBranches.remove(mainBranch)
                    for branch in remainingBranches:
                        stroke = RiverStroke(self, branch)
                        stroke.setBraided(True)
                        self.strokes.append(stroke)
                        self.groupedFeatures.append(branch)
                        if (not branch['geom'].coords[-1] in downstreamNodes):
                            
                            downstreamNodes.insert(0, branch['geom'].coords[-1])
                        #compute Strahler orders
                        self.strahlerOrders[branch['id']]= 1
            a=int(len(self.groupedFeatures)*100/total)
            preva=a

    def reconstruct_strokes(self):
        array = []
        for i, stroke in enumerate(self.strokes):
            listline=[]
            section=0
            for j, seg in enumerate(stroke.features):
                            listline+=[seg["geom"]]
                            section=max(section,self.strahlerOrders[seg['id']])
            multi_line = MultiLineString(listline)
            merged_line = ops.linemerge(multi_line)
            array += [[i,merged_line,section]]
        return array
    
    def save_strokes_shp(self,path):
        strokes_list = self.strokes
        array = self.__reconstruct_strokes(strokes_list)
        gdf = gpd.GeoDataFrame(array,  columns = ['id', 'geom',"order"],crs="epsg:2154",geometry="geom")   
        gdf.to_file(path, driver='ESRI Shapefile')
        
    def add_strahler_and_strokes(self,path):
        gdf = self.dtf 
        strahlerL=[]
        strokeL=[]
        for i in gdf.index: 
            strahlerL.append(self.strahlerOrders[gdf["id"][i]])
            strokeid=self.getstrokeid(gdf["id"][i])
            strokeL.append(strokeid)
        gdf['strahler'] = strahlerL
        gdf['stroke'] = strokeL
        gdf.to_file(path, driver='ESRI Shapefile')

    def get_strokes(self):
        return self.__reconstruct_strokes(self.strokes)
    
    def getstrokeid(self,ids):
        i=0
        for elem in self.strokes:
            for section in elem.features:
                if section['id']==ids: 
                    return i
            i+=1
        return -1