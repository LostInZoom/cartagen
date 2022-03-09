from evaluation.constraint_monitors import building_size_monitor
import geopandas

def test_building_monitors():
    filename = "D:\Donnees\Benchmarks\HuezTopo\IGN_BUILDING.zip"
    df = geopandas.read_file(filename, encoding='latin-1')
    for geom in df['geometry']:
        monitor = building_size_monitor.BuildingSizeMonitor(40.0,geom,geom,1)
        print(monitor.get_satisfaction())
    return

def main():
    test_building_monitors()

if __name__ == "__main__":
    main()
