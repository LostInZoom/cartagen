Processes
=========

AGENT
~~~~~

Agent-based generalization follows a step by step, knowledge based approach,
which arises from the seminal model from Ruas & Plazanet :footcite:p:`ruas:1996`.
The aim is not to find a global solution to generalize the whole map at once,
but to identify local situations that require generalization,
and try generalization algorithms to improve the legibility of this situation.

This process is a port of the version present in the
`CartAGen Java application <https://ignf.github.io/CartAGen/>`_.
For a better understanding of the AGENT process, you can refer to
Barrault *et al.* :footcite:p:`barrault:2001`, Ruas & Duchêne :footcite:p:`ruas:2007`
and Duchêne *et al.* :footcite:p:`duchene:2018`

Micro-agents
------------

Micro-agents are the simpler kind of agents as they aim at generalising single
cartographic features without consideration for their surroundings.

- Before you can initiate the AGENT process, you must first prepare the data:

.. plot::
    :nofigs:
    :context: reset
    :include-source: True
    :show-source-link: False

    import geopandas as gpd
    import cartagen as c4
    from shapely.wkt import loads

    # Given the following buildings
    buildings = [
        loads('Polygon ((395038.7 6272970.9, 395030.4 6272984, 395025.3 6272982, 395023.2 6272983.7, 395020 6272981.3, 395016.9 6272985.9, 395021.8 6272990.7, 395020.6 6272993.7, 395024.7 6272997.2, 395028.5 6272994.5, 395032.8 6272988.2, 395038.1 6272991.6, 395044.9 6272979.1, 395047.1 6272980.4, 395049.5 6272976.8, 395038.7 6272970.9))'),
        loads('Polygon ((394999.5 6272975, 395006.7 6272962.4, 395010.6 6272957.5, 394996.6 6272944.4, 394991 6272949, 394999.2 6272956.3, 394996.1 6272959.7, 394998.3 6272961.3, 394992 6272969.4, 394999.5 6272975))'),
        loads('Polygon ((395007.3 6272975.8, 395013.2 6272981, 395021.2 6272969.6, 395024.2 6272971.9, 395031 6272963.8, 395020.8 6272957.4, 395007.3 6272975.8))'),
        loads('Polygon ((395082.3 6272967.4, 395089.9 6272958, 395071.9 6272945.9, 395068.4 6272950.6, 395066 6272949, 395056.3 6272962, 395058.5 6272963.5, 395056.40000000002328306 6272966.8, 395059.4 6272969.9, 395056.9 6272972.6, 395054.5 6272968.3, 395049.6 6272973.4, 395058.4 6272981.6, 395073.6 6272962.5, 395082.3 6272967.4))')
    ]

    # Create a GeoDataFrame from the buildings
    gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(buildings))

- Then, for every building inside the GeoDataFrame, we can create an agent along with
  constraints that will guide the generalisation process.

.. list-table:: Building micro constraints
    :widths: 50 20 50
    :header-rows: 1

    * - name
      - property
      - actions
    * - BuildingSizeConstraint
      - area
      - enlarge, delete
    * - BuildingGranularityConstraint
      - granularity
      - simplify, simplify to rectangle
    * - BuildingSquarenessConstraint
      - squareness
      - squaring

.. plot::
    :nofigs:
    :context:
    :include-source: True
    :show-source-link: False

    agents = []
    for i, building in gdf.iterrows():
        # Create the agent
        agent = c4.BuildingAgent(building)

        # Adding a size constraint of 250 square meters to enlarge the building if needed
        size = c4.BuildingSizeConstraint(agent, importance=1, min_area=250)

        # Adding a granularity constraint that will simplify the building if
        # an edge is above the minimum allowed length
        granularity = c4.BuildingGranularityConstraint(agent, importance=1, min_length=6)

        # Adding a squareness constraint to square the building angles if needed
        squareness = c4.BuildingSquarenessConstraint(agent, importance=1)
        
        agent.constraints.append(size)
        agent.constraints.append(squareness)
        agent.constraints.append(granularity)
        agents.append(agent)
    
    c4.run_agents(agents)    

.. plot::
    :caption: Resulting generalisation using the micro-agents
    :context:

    import numpy
    from matplotlib import pyplot as plt
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    fig = plt.figure(1, (10, 6))
    sub1 = fig.add_subplot(111)
    sub1.axes.get_xaxis().set_visible(False)
    sub1.axes.get_yaxis().set_visible(False)

    for b in buildings:
        poly = Path.make_compound_path(Path(numpy.asarray(b.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in b.interiors])
        sub1.add_patch(PathPatch(poly, facecolor="lightgray", edgecolor='none'))

    for b in gdf.geometry:
        poly = Path.make_compound_path(Path(numpy.asarray(b.exterior.coords)[:, :2]),*[Path(numpy.asarray(ring.coords)[:, :2]) for ring in b.interiors])
        sub1.add_patch(PathPatch(poly, facecolor="none", edgecolor='red'))
    
    sub1.autoscale_view()
    plt.show()