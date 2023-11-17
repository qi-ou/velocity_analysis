from merge_tif import *
import geopandas as gpd
import shapely.speedups
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
shapely.speedups.enable()
from pykrige.uk import UniversalKriging
from pykrige.kriging_tools import write_asc_grid
from matplotlib import cm
import numpy as np
import itertools
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from 3_invert_gps_overlap_with_dummy import *

# def load_chris():
#     """ Load GPS data into a geopandas dataframe
#     - only need lon, lat, vel and sig for the purpose of interpolation"""
#     gps_gdf = gpd.GeoDataFrame(crs="EPSG:4326")
#     gps_gdf['geometry'] = None
#     index = 0
#     gps_2d_file = "../gps/tianshan_tol1.5_minocc2.5_dist2_2D_dec2022.dat"
#
#     fl = open(gps_2d_file, "r").readlines()
#     for line in fl:
#         if not line.startswith(('*', "Table")):
#             lon, lat, Ve_E, Vn_E, dVe, dVn, Cen, sta = line.split()
#             gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
#             gps_gdf.loc[index, 'station'] = sta[:4]
#             gps_gdf.loc[index, 've'] = float(Ve_E)  # eastern velocity in mm/yr in fixed eurasia reference frame
#             gps_gdf.loc[index, 'vn'] = float(Vn_E)  # northern velocity in mm/yr in fixed eurasia reference frame
#             gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
#             gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
#             gps_gdf.loc[index, 'cen'] = float(Cen)  # correlation between east and northern velocities
#             gps_gdf.loc[index, 'vu'] = float(0)  # up velocity in mm/yr in ITRF20082008 reference frame
#             gps_gdf.loc[index, 'su'] = float(0)  # sigma vu in ITRF20082008 reference frame
#             index += 1
#
#     gps_3d_file = "../gps/tianshan_tol1.5_minocc2.5_dist2_3D_dec2022.dat"
#     fl = open(gps_3d_file, "r").readlines()
#     for line in fl:
#         if not line.startswith(('*')):
#             lon, lat, Ve_E, Vn_E, Vu, dVe, dVn, dVu, _, _, _, sta = line.split()
#             if sta in gps_gdf['station'].values:  # only input vertical if Wang has the same station
#                 index = gps_gdf.index[gps_gdf['station'] == sta].tolist()[0]  # get index in dataframe based on station
#                 gps_gdf.loc[index, 'vu'] = float(Vu)  # up velocity in mm/yr relative to stable northern neighbour
#                 gps_gdf.loc[index, 'su'] = float(dVu)  # sigma vu relative to stable northern neighbour
#     return gps_gdf
#
#
# def load_compiled_vu_su():
#     """ Load GPS data into a geopandas dataframe
#     - only need lon, lat, vel and sig for the purpose of interpolation"""
#     gps_gdf = gpd.GeoDataFrame(crs="EPSG:4326")
#     gps_gdf['geometry'] = None
#     index = 0
#     gps_file = "/Users/qi/OneDrive - University of Leeds/projects/insar/velocity_analysis/decompose/compiled_vu.txt"
#
#     fl = open(gps_file, "r").readlines()
#     for line in fl:
#         if not line.startswith(('lon', 'Lon')):
#             lon, lat, Vu, Su = line.split()
#             gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
#             gps_gdf.loc[index, 'vu'] = float(Vu)  # eastern velocity in mm/yr in fixed eurasia reference frame
#             gps_gdf.loc[index, 'su'] = float(Su)  # sigma ve
#             index += 1
#     return gps_gdf


if __name__ == "__main__":

    # define map extent with longitude and latitude in degrees
    west = 63
    east = 98
    south = 36
    north = 52
    out_dir = "../los_weighted/decompose/"

    # fault
    fault_file = "/Users/qi/OneDrive - University of Leeds/datasets/vectors/faults/gem-global-active-faults-master/kml/gem_active_faults_harmonized.kml"
    faults = gpd.read_file(fault_file, driver='KML')
    fig, ax = plt.subplots(3, 2, sharex='all', sharey='all')
    for x in ax.flatten():
        faults.plot(ax=x, linewidth=0.5, color='grey')
        x.set_xlim((west, east))
        x.set_ylim((south, north))

    # load and compare InSAR and GPS Ve
    ve = OpenTif("../los_weighted/decompose/ve_tuned.tif")
    ve.data = ve.data
    culling_thresh = 0.7
    all_gps = load_chris()
    map_box = Polygon([(west, south), (east, south), (east, north), (west, north)])
    gps_in_box = all_gps[all_gps.within(map_box)]
    # data culling, get rid of contradictory clustering values by removing data with high uncertainties
    gps = gps_in_box[gps_in_box["se"] <= culling_thresh]
    gps.loc[:, 'InSAR_ve'] = [ve.extract_pixel_value(point.x, point.y)[0] for point in gps['geometry']]
    gps["InSAR-GNSS Ve"] = gps['InSAR_ve'] -gps["ve"]
    gps.dropna(inplace=True)

    # plot Ve
    gps.plot("ve", ax=ax[0,0], vmin=-10, vmax=10, cmap='bwr', legend=True)
    ax[0,0].set_title("GNSS Ve")
    gps.plot("InSAR_ve", ax=ax[1,0], vmin=-10, vmax=10, cmap='bwr', legend=True)
    ax[1,0].set_title("InSAR Ve")
    gps.plot("InSAR-GNSS Ve", ax=ax[2,0], vmin=-10, vmax=10, cmap='bwr', legend=True)
    ax[2,0].set_title("InSAR-GNSS Ve")
    # plt.show()

    # load and compare InSAR and GPS Vu
    vu = OpenTif("../los_weighted/decompose/vu_tuned.tif")
    vu.data = vu.data
    culling_thresh = 1.5
    # all_gps = load_chris()
    map_box = Polygon([(west, south), (east, south), (east, north), (west, north)])
    gps_in_box = all_gps[all_gps.within(map_box)]
    # data culling, get rid of contradictory clustering values by removing data with high uncertainties
    gps3d = gps_in_box[gps_in_box["su"] <= culling_thresh]
    gps3d.loc[:, 'InSAR_vu'] = [vu.extract_pixel_value(point.x, point.y)[0] for point in gps3d['geometry']]
    gps3d["InSAR-GNSS Vu"] = gps3d['InSAR_vu'] -gps3d["vu"]
    gps3d.dropna(inplace=True)

    # plot Vu
    gps3d.plot("vu", ax=ax[0,1], vmin=-10, vmax=10, cmap='bwr', legend=True)
    ax[0,1].set_title("GNSS Vu")
    gps3d.plot("InSAR_vu", ax=ax[1,1], vmin=-10, vmax=10, cmap='bwr', legend=True)
    ax[1,1].set_title("InSAR Vu")
    gps3d.plot("InSAR-GNSS Vu", ax=ax[2,1], vmin=-10, vmax=10, cmap='bwr', legend=True)
    ax[2,1].set_title("InSAR-GNSS Vu")

    plt.show()

    fig, ax = plt.subplots(1,2, figsize=(6.4, 1.8))
    # plot histogram for Vu
    ax[1].set_title("InSAR-GNSS Vu")
    count, bin_edge, _ = ax[1].hist(gps3d["InSAR-GNSS Vu"], bins=np.arange(-10, 10, 0.5))
    nanmode = (bin_edge[np.argmax(count)] + bin_edge[np.argmax(count)+1] ) /2
    ax[1].annotate("mean: %.1f \n std: %.1f " % (np.nanmean(gps3d["InSAR-GNSS Vu"]), np.nanstd(gps3d["InSAR-GNSS Vu"])),
                xy=(0.03, 0.93), xycoords='axes fraction',
                ha="left", va="top",
                fontsize=12)
    ax[1].set_xlabel("mm/yr")

    # plot histogram for Ve
    count, bin_edge, _ = ax[0].hist(gps["InSAR-GNSS Ve"], bins=np.arange(-10, 10, 0.5))
    nanmode = (bin_edge[np.argmax(count)] + bin_edge[np.argmax(count)+1] ) /2
    ax[0].set_title("InSAR-GNSS Ve")
    ax[0].annotate("mean: %.1f \n std: %.1f " % (np.nanmean(gps["InSAR-GNSS Ve"]), np.nanstd(gps["InSAR-GNSS Ve"])),
                xy=(0.03, 0.93), xycoords='axes fraction',
                ha="left", va="top",
                fontsize=12)
    ax[0].set_xlabel("mm/yr")
    plt.show()