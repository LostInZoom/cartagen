# -*- coding: utf-8 -*-
import math

import numpy as np
from shapely.geometry import LineString, Point, MultiLineString
from cartagen4py.utils.geometry import angle_between_2lines
from shapely import ops
import geopandas as gpd
import networkx as nx

class RiverStroke:
    
    def __init__(self,network,section): 
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
            elem["id"]=lines.id[idx]
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


    def allBelongToStroke(self,sections) :
        for section in sections:
            if (not section in self.groupedFeatures) :
                return False
        return True
    
    def getNonBraidedStroke(self,upstreamStrokes) :
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
            if (not downstreamSection["name"] is None) :        
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
        #Save a shapefile of stroke geometry in the desired folder
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

                



        
