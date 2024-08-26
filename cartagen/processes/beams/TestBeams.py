import geopandas

zipfile = "zip://C:/Users/gtouya/workspace/CartAGen4Py/data/beams_example_rivers.zip"
df = geopandas.read_file(zipfile)
print(df)
df.plot()