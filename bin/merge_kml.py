import geopandas as gpd
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
import matplotlib.pyplot as plt

# frame_file = '/Users/qi/OneDrive - University of Leeds/projects/insar/velocity_analysis/kml/TS_asc.kml'
frame_file = '/Users/qi/Desktop/129_twoFrames.kml'
# load frame kml into geopandas dataframe
frames = gpd.read_file(frame_file, driver='KML')
frames['track'] = frames['Name'].str[:4]
tracks = frames.dissolve(by='track', aggfunc='sum')
# tracks.to_file('/Users/qi/OneDrive - University of Leeds/projects/insar/velocity_analysis/kml/TS_asc_track.kml', driver='KML')
track = tracks.loc[tracks['Name'].str[:4] == "129A"]
wkt = [geom.wkt for geom in track.geometry]

