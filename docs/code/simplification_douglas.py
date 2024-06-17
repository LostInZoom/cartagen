from matplotlib import pyplot as plt
import cartagen4py as c4
import geopandas as gp

# Loading the lines as a GeoDataFrame
lines = gp.read_file('data/test_lines.geojson')

# Creating an empty list to store simplified lines
slines = []

# Looping through the geometry of each lines
for line in lines.geometry:
    # Simplifying lines using the Douglas Peucker algorithm with an area of 200
    simplified = c4.douglas_peucker(line, 5)
    # Adding the simplified line to list
    slines.append(simplified)

# Creating a new GeoSeries from the list of simplified lines
generalized = gp.GeoSeries(slines)

# Plotting the original and the simplified lines
original = lines.plot(edgecolor='gray',linewidth=1)
generalized.plot(ax=original, edgecolor='red', linewidth=1)
original.axes.get_xaxis().set_visible(False)
original.axes.get_yaxis().set_visible(False)

plt.show()