import math

import numpy as np
from shapely.geometry import LineString, Point, MultiLineString
from shapely import ops
import geopandas as gpd
import networkx as nx

from cartagen4py.utils.geometry.angle import angle_3_pts

class Stroke:
    COUNTER = 0
    def __init__(self, network, root):
        self.network = network
        self.root=root
        self.features = []
        self.features.append(root)
        self.id = Stroke.COUNTER
        Stroke.COUNTER += 1
    
    def is_good_continuity(temp_geom, geom_foll,deviatAngle, deviatSum):
        angleThresh = deviatAngle / 180.0 * math.pi;
        sumThresh = deviatAngle / 180.0 * math.pi;
    
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
          return -1.0;
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

    def filterFollowers(self,arc, followers):
        loopFoll = []
        #loop on the followers to filter them
        loopFoll+=followers
        for a in loopFoll:
            if (a==arc):
                #remove it from the set
                followers.remove(a);
                continue;
             #check that a does not belong to another stroke
            if (a in self.network.groupedFeatures):
                  #remove it from the set
                  followers.remove(a);
                  continue;
                  #check if it belongs to this stroke
            if (a in self.features):
             #remove it from the set
                 followers.remove(a);
                 continue;
            #check if it belongs to the network
            if (not a in self.network.features) :
            #// remove it from the set
                followers.remove(a);
      
                        
    def chooseNextSegment(self,arc, followers, attributeNames,deviatAngle, deviatSum):
        #first, if it's a node of degree two
        if (len(followers) == 1) :
            follower = next(iter(followers))
            if (not follower in self.features):#COR
                return follower;
            return None
        #then, filter the followers
        self.filterFollowers(arc, followers)
        if (len(followers) == 0) :
            return None
        #then, filter the followers from the attributeNames
        Stroke.attribute_filter(arc, followers, attributeNames)
        if (len(followers) == 0) :
            return None
        continuity,bestSegment = True, None
        #Loop on the followers to choose the best continuity
        minDiff = math.pi
        for follower in followers:
            #get the continuity difference with this follower
            diffContinuity = Stroke.is_good_continuity(arc, follower, deviatAngle, deviatSum) 
            if (diffContinuity > -1.0) :
                  continuity = True;
                  if (diffContinuity < minDiff) :
                      #this is the current best continuity update the difference
                      minDiff = diffContinuity;
                      #change the bestSegment
                      bestSegment = follower;
        if (((continuity) and not bestSegment in self.network.groupedFeatures) or (minDiff < deviatAngle)) :
            return bestSegment;
        return None

  
    def one_side_stroke(self, side, attributeNames,deviatAngle,deviatSum):
        #get the following network segments of the root of this stroke
        node =None        
        node = self.root["geom"].coords[side]
        if (node is None):
            return;
        followers =[]
        followers+=self.network.dic_neighbours[self.root["geom"].coords[side]]
        if self.root in followers:
            followers.remove(self.root)
        next1 = self.root
        continuity = True
        self.network.groupedFeatures.append(next1)
        while (continuity) :
                #get the best candidate among the followers (the one with best continuity)
                best = self.chooseNextSegment(next1, followers, attributeNames,deviatAngle, deviatSum)
                if (best is None):
                    break
                #      // add this to the 2 sets (the network one and the stroke one)
                if not side:
                    self.features.insert(0, best)
                else:
                    self.features.append(best)
                self.network.groupedFeatures.append(best) 
                #get the followers of 'best'
                followers=[]
                nextNode = best["geom"].coords[0];
                
                side2=0
                if (node==nextNode):
                    nextNode = best["geom"].coords[-1];
                    side2=-1 
                    
                followers+=self.network.dic_neighbours[best["geom"].coords[side2]]#AC
                followers.remove(best)
                #if there is no follower, break
                if (len(followers)== 0):
                    break;
                #update the 'next' segment with 'best'
                next1 = best;
                node = nextNode;         
    def __str__(self):
        liste=""
        for elem in self.features: 
            liste+=str(elem["id"])
            liste+=","
        return liste
    
class StrokeNetwork:
    def __init__(self,features):
            #Initialisation from a list of shapely geometry with the correct attributes
          self.features = features
          self.groupedFeatures = []
          self.id = 0;
          self.strokes = []
          
    def __init__(self, shapefile, attributeNames):
        #Initialisation from a geopanda dataframe and the liste of desired attribute name
        features=[]
        for idx in shapefile.index: 
            elem={}
            elem["geom"]=shapefile.geometry[idx]
            elem["id"]=shapefile.id[idx]

            for attr in attributeNames:
                elem[attr]=getattr(shapefile,attr)[idx]
                #setattr(elem,attr,)
            features+=[elem]
        self.features = features
        self.groupedFeatures = []
        self.id = 0;
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
    #loop on the network features
        for obj in self.features :
                #test if the feature has already been treated
                if obj in self.groupedFeatures :
                    continue
                #build a new stroke object
                stroke = Stroke(self, obj)
                #// build the stroke on the initial side
                stroke.one_side_stroke(0, attributeNames, deviatAngle, deviatSum);
                #// build the stroke on the final side
                stroke.one_side_stroke(-1, attributeNames, deviatAngle, deviatSum);
                #// add the stroke to the strokes set
                self.strokes.append(stroke);
    
    def __reconstruct_strokes(self, strokes):
        array = []
        for i, stroke in enumerate(strokes):
            listline=[]
            section=""
            for j, seg in enumerate(stroke.features):
                listline+=[seg["geom"]]
                section += str(seg["id"])
                if j < len(stroke.features):
                    section += ","
            print(listline)
            multi_line = MultiLineString(listline)
            merged_line = ops.linemerge(multi_line)
            array += [[i,merged_line,section]]
        return array

    def save_strokes_shp(self,path):
        #Save a shapefile of stroke eometry in the desired folder
        strokes_list = self.strokes
        array = self.__reconstruct_strokes(strokes_list)
        gdf = gpd.GeoDataFrame(array,  columns = ['id', 'geom',"section"],crs="epsg:2154",geometry="geom")   
        gdf.to_file(path, driver='ESRI Shapefile')

    def get_strokes(self):
        return self.__reconstruct_strokes(self.strokes)

                



        
