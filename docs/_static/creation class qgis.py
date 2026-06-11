import re

# Helpstring to Html
def transfoHtml(helpstring, doc_url):

    lines = [line.strip() for line in helpstring.splitlines() if line.strip()]

    description = []
    parameters = []
    nbParam = 0


    for line in lines:

        if re.findall("[–]", line):
            nbParam +=1
            if nbParam != 1:
                name, descrip = line.split("–", 1)
                name, type = name.split("(", 1)
                # print(f"nom : {name} et description {descrip}")
                parameters.append((name, descrip))
        else :
            description.append(line)

    for i in range(3):
        description.remove(description[-1])

    description[0] = description[0] + '\n'
    html = " ".join(description)

    if parameters:
        html += "\n         <h3> Parameters: </h3>\n        <ul>\n"

        for name, descrip in parameters:
            html += f'          <li> - <em>{name.capitalize()}</em> : {descrip} </li>\n'

        html += "        </ul>"

    if doc_url:
        html += (
            f"""
            For more see <a href="{doc_url}">help online</a>.
            """)
    
    return(html)

#Main Process
def creationQgisProcess(nom:str, helpstring:str, group:str, url:str, cheminFichier:str):
    """
    Fonction qui créer le texte pour la Classe qui fait tourner les algorithmes
    de CartAgen dans QGis
    """

    nom = nom.strip()

    #Nom avec un espace
    nbMaj = len(re.findall('[A-Z]',nom))
    listNom = re.split('([A-Z])', nom)
    midpoint = len(listNom)//nbMaj+1
    if nbMaj == 2:
        nomEsp = listNom[0:midpoint] + [" "] + listNom[midpoint:]
    if nbMaj == 3:
        nomEsp = listNom[0:midpoint] + [" "] + listNom[midpoint:midpoint+2] + [" "] + listNom[midpoint+2:]
    nomEsp = "".join(nomEsp)

    #Nom de sortie
    nomFin = listNom[0:midpoint] + ["ed "] + listNom[midpoint:]
    nomFin = "".join(nomFin)

    #Nom avec des tirets
    listNom = re.split('([A-Z])', nom)
    midpoint = len(listNom)//2+1
    nomTiret = listNom[0:midpoint] + ["_"] + listNom[midpoint:]
    nomTiret = "".join(nomTiret)
    nomTiret = nomTiret.lower()

    #Création du fichier
    chemin = f"{cheminFichier}{nomTiret}.py"
    fichier = open(chemin, "x")    

    #Nettoyage de l'helpstring
    tableauHelp = helpstring.splitlines()
    listHelp = []
    for ligne in tableauHelp:
        ligne = ligne.strip()
        if ligne != "":
            listHelp.append(ligne)

    #Récupération des paramètres et de leur type
    parameters = []
    nbParam = 0
    paramType = []
    listDefaultValue = []
    for i, ligne in enumerate(listHelp):
        if re.findall("[–]", ligne):
            nbParam = nbParam + 1
            if nbParam != 1 and i!=len(listHelp):
                name, description = ligne.split("–", 1)

                #Isolation de la valeur par défaut
                try:
                    resid, defaultValue = description.split("Default")
                    defaultValue = ''.join(defaultValue)
                    defaultValue = defaultValue.strip()
                    resid, defaultValue = defaultValue.split(" ", 1)
                    defaultValue = re.sub(r'[.]', '', defaultValue)
                    listDefaultValue.append(defaultValue)
                except:
                    listDefaultValue.append(1)

                #Isolation nom paramètre
                name, type = name.split("(", 1)
                parameters.append(name.strip())

                #Isolation du type
                if len(type.split(",",1))==2:
                    type, residus = type.split(",",1)
                    paramType.append((name.strip(), type.strip()))
                elif len(type.split(",",1))==1:
                    type, residus = type.split(")",1)
                    paramType.append((name.strip(), type.strip()))

            else :
                #Isolation du nom de la géométrie d'entrée
                name, description = ligne.split("–", 1)
                name, type = name.split("(", 1)
                type = type.strip()
                type = type[:-1]
                type = type.split()
                nomTypeGeom = (type, description)

    #Ecriture des types de paramètres du Sink
    typeGeom = nomTypeGeom[0]
    compteurLine = 0
    compteurPolyg = 0
    compteurPoint = 0
    typeSink = []
    for i, geom in enumerate(typeGeom):
        if re.findall('Point', geom):
            if compteurPoint == 0:
                typeSink.append('QgsProcessing.TypeVectorPoint')
            compteurPoint+=1
        if re.findall('Line', geom):
            if compteurLine == 0:
                typeSink.append('QgsProcessing.TypeVectorLine')
            compteurLine+=1
        if re.findall('Polygon', geom):
            if compteurPolyg == 0:
                typeSink.append('QgsProcessing.TypeVectorPolygon')
            compteurPolyg+=1
    typeSink = ','.join(typeSink)

    #Param pour les algos
    listParamAlgo = []
    for parameter in parameters :
        texte = f"{parameter}={parameter}"
        listParamAlgo.append(texte)
    paramAlgo=",".join(listParamAlgo)

    #Creation de la partie avec les paramètres en majuscule
    paramUp=""
    for parameter in parameters:
        up = parameter.upper()
        paramUp = paramUp+f"    {up}='{up}'\n"

    #Creation de l'helpstring telle quelle sera affichée
    helpStringHtml = transfoHtml(helpstring, url)

    initParam = f"""                
        # We add the input vector features source.
        input = QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('{nomTypeGeom[1]} :'),
                [{typeSink}]
            )
        self.addParameter(input)
        """
    for i, test in enumerate(paramType) :
        parameter, type = test
        if type == "float":
            paramTxt = f"""
        {parameter} = QgsProcessingParameterNumber(
            self.{parameter.upper()},
            self.tr('{parameter.capitalize()} :'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue={listDefaultValue[i]},
            optional=False
        )
        self.addParameter({parameter})
        """
            initParam += paramTxt

        elif type == "bool":
            paramTxt = f"""
        {parameter} = QgsProcessingParameterBoolean(
            self.{parameter.upper()},
            self.tr('{parameter.capitalize()} ?'),
            defaultValue=True,
            optional=False
        )
        self.addParameter({parameter})
            """
            initParam += paramTxt
        
        elif type == "int":
            paramTxt = f"""
        {parameter} = QgsProcessingParameterNumber(
            self.{parameter.upper()},
            self.tr('{parameter.capitalize()} :'),
            type=QgsProcessingParameterNumber.Integer,
            defaultValue={listDefaultValue[i]},
            optional=False
        )
        self.addParameter({parameter})
            """
            initParam += paramTxt

        elif type == "str":
            paramTxt = f"""
        {parameter}s = []
        {parameter} = QgsProcessingParameterEnum(
            self.{parameter.upper()},
            self.tr('{parameter.capitalize()} ?'),
            {parameter}s
            defaultValue=30.0,
        )
        self.addParameter({parameter})
            """
            initParam += paramTxt

        else :
            paramTxt = f"""
        {parameter} = QgsProcessingParameterFeatureSource(
            name=self.{parameter.upper()},
            description="{parameter.capitalize()} :",
            types=[QgsProcessing.TypeVector{type.split()[-1]}]
        )
        self.addParameter({parameter})
            """
            initParam += paramTxt

    #Creation des paramètres pour le Process
    processParam=""
    for parameter, type in paramType :
        if type == "float":
            paramTxt = f"""
        {parameter} = self.parameterAsDouble(parameters, self.{parameter.upper()}, context)
            """
            processParam += paramTxt

        elif type == "bool":
            paramTxt = f"""
        {parameter} = self.parameterAsBoolean(parameters, self.{parameter.upper()}, context)
            """
            processParam += paramTxt
        
        elif type == "int":
            paramTxt = f"""
        {parameter} = self.parameterAsInt(parameters, self.{parameter.upper()}, context)
            """
            processParam += paramTxt

        elif type == "str":
            paramTxt = f"""
        {parameter} = self.parameterAsEnum(parameters, self.{parameter.upper()}, context)
            """
            processParam += paramTxt

        else :
            paramTxt = f"""
        {parameter} = self.parameterAsSource(parameters, self.{parameter.upper()}, context)"""
            processParam += paramTxt

    modele = f'''
from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessing, QgsFeatureSink, QgsProcessingAlgorithm,
    QgsFeature, QgsGeometry, QgsProcessingParameterDefinition,
    QgsProcessingException, QgsWkbTypes,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterNumber,
    QgsProcessingParameterDistance,
)    

class {nom} (QgsProcessingAlgorithm):
    """
            {helpstring}
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    OUTPUT = 'OUTPUT'
    INPUT = 'INPUT'
{paramUp}
    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return '{nomEsp}'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr(self.name())

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return '{group}'

    def icon(self):
        """
        Should return a QIcon which is used for your provider inside
        the Processing toolbox.
        """
        from cartagen4qgis import get_plugin_icon
        return get_plugin_icon()

    def shortDescription(self):
        """
        Returns an optional translated short description of the algorithm. This 
        should be at most a single sentence, e.g. “Converts 2D features to 3D by 
        sampling a DEM raster.”
        """
        first_line = self.shortHelpString().strip().splitlines()[0]
        description = self.tr(first_line)
        
        return(description)

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        helpstring = """
        {helpStringHtml}
        """
        
        return self.tr(helpstring)
    
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return {nom}()

    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
{initParam}
        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).   
        output = QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('{nomFin}'))
        self.addParameter(output)

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        import geopandas as gpd
        import pandas
        from cartagen import {nomTiret}
        from cartagen4qgis.src.tools import list_to_qgis_feature_2

        # Retrieve the feature source and sink. The 'dest_id' variable is used
        # to uniquely identify the feature sink, and must be included in the
        # dictionary returned by the processAlgorithm function.
        source = self.parameterAsSource(parameters, self.INPUT, context)
        gdf = gpd.GeoDataFrame.from_features(source.getFeatures())
        
        # retrieve the other parameters values
{processParam}
        
        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / source.featureCount() if source.featureCount() else 0
        features = source.getFeatures()

        dp = gdf.copy()
        for i in range(len(gdf)):
            dp.loc[i,'geometry'] = {nomTiret} (list(gdf.geometry)[i], {paramAlgo})

            # Update the progress bar
            feedback.setProgress(int(i * total))

        res = dp.to_dict('records')
        res = list_to_qgis_feature_2(res,source.fields())

        # Create the output sink    
        (sink, dest_id) = self.parameterAsSink(parameters, self.OUTPUT,
                context, res[0].fields(), source.wkbType(), source.sourceCrs())
        
        # Add a feature in the sink
        sink.addFeatures(res, QgsFeatureSink.FastInsert)
        '''
    modele += r'''
        return {
            self.OUTPUT: dest_id
            }
        '''
    with open(chemin, 'w') as f:
        f.write(modele)    

    return modele

#LineString of the name of the algorithm. It needs to be written this way : GeneralizeAreaPatches, SpinalizePolygon
name = "SimplifyRaposo"

#A simple copy and paste of the documentation. The parenthesis (), hyphen –, and "Default is" are needed for the algorithm to fonctionne. 
helpstring = """
Simplify a line or a polygon using an hexagonal tessellation.

This algorithm proposed by Raposo simplifies lines based on a hexagonal tessellation. The algorithm also works for the simplification of the border of a polygon object. The idea of the algorithm is to put a hexagonal tessellation on top of the line to simplify, the size of the cells depending on the targeted granularity of the line. Similarly to the Li-Openshaw algorithm, only one vertex is kept inside each cell. This point can be the centroid of the removed vertices, or a projection on the initial line of this centroid. The shapes obtained with this algorithm are less sharp than the ones obtained with other algorithms such as Douglas-Peucker.

The algorithm is dedicated to the smooth simplification of natural features such as rivers, forests, coastlines, lakes.

Parameters:

        geometry (LineString, MultiLineString, Polygon, MultiPolygon, LinearRing) – The line to simplify.

        initial_scale (float) – Initial scale of the provided line (25000.0 for 1:25000 scale).

        final_scale (float) – Final scale of the simplified line.

        centroid (bool, optional) – If true, uses the center of the hexagonal cells as the new vertex, if false, the center is projected on the nearest point in the initial line.

        tobler (bool, optional) – If True, compute cell resolution based on Tobler’s formula, else uses Raposo’s formula

Returns:

    LineString, MultiLineString, Polygon, MultiPolygon, LinearRing
    """

#The group in which you want to put the algorithm. You can check in QGis what groups are avalailable
group = "Simplify lines and patches"

#The url of the fonction documentation
urlDoc = "https://cartagen.readthedocs.io/en/latest/reference/cartagen.simplify_raposo.html#cartagen.simplify_raposo"

#Where you want the end file to be. The name of the file is created during the process. 
cheminFichier = "C:/Users/Desktop/"

print(creationQgisProcess(name, helpstring, group, urlDoc, cheminFichier))
