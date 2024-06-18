from merge_tif import *
import geopandas as gpd
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
import matplotlib.pyplot as plt
from 4_decompose_into_ve_vu_using_interpolated_vn_and_NEU import *
from compare_gps_insar import *
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


def add_NEU_los(gps_df, N_coef_map, E_coef_map, U_coef_map):
    # extract pixel values from inc and heading based on gps lon lat and insert into geodataframe
    gps_df.loc[:, 'n_coef'] = [N_coef_map.extract_pixel_value(point.x, point.y)[0] for point in gps_df['geometry']]
    gps_df.loc[:, 'e_coef'] = [E_coef_map.extract_pixel_value(point.x, point.y)[0] for point in gps_df['geometry']]
    gps_df.loc[:, 'u_coef'] = [U_coef_map.extract_pixel_value(point.x, point.y)[0] for point in gps_df['geometry']]
    # calculate los and los uncertainty per gps point
    gps_df.loc[:, 'los'] = [ vn * n + ve * e + vu * u for vn, ve, vu, n, e, u
                            in gps_df[['vn', 've', 'vu', 'n_coef', 'e_coef', 'u_coef']].to_numpy()]
    gps_df.loc[:, 'los_sigma'] = [ abs(sn * n) + abs(se * e) + abs(su * u) for sn, se, su, n, e, u
                            in gps_df[['sn', 'se', 'su', 'n_coef', 'e_coef', 'u_coef']].to_numpy()]


if __name__ == "__main__":
    # load frame kml into geopandas dataframe
    frame_file = '../kml/TS_des.kml'
    # frame_file = '../kml/TS_asc.kml'
    frames = gpd.read_file(frame_file, driver='KML')
    frames['track'] = frames['Name'].str[:4]
    tracks = frames.dissolve(by='track', aggfunc='sum')
    outdir = "../los_tracks_no_plate_shifted_to_2dGPS/"

    # load gps
    all_gps = load_chris()
    # map_box = Polygon([(west, south), (east, south), (east, north), (west, north)])
    # gps_in_box = all_gps[all_gps.within(map_box)]

    for track in tracks['Name'].str[:4]:
        tif = OpenTif("../los_tracks_no_plate/" + track + "__mode.tif")
        sigfile = OpenTif("../corrected_vstd_tracks/" + track + "__none.tif")
        N = OpenTif("../NEU_tracks/" + track + "__N.tif")
        E = OpenTif("../NEU_tracks/" + track + "__E.tif")
        U = OpenTif("../NEU_tracks/" + track + "__U.tif")

        poly = tracks[tracks['Name'].str[:4] == track]
        poly.reset_index(drop=True, inplace=True)
        gps = all_gps[all_gps.within(poly.loc[0, 'geometry'])]
        add_NEU_los(gps, N, E, U)

        # breakpoint()
        gps.loc[:, 'InSAR'] = [tif.extract_pixel_value(point.x, point.y)[0] for point in gps['geometry']]
        gps['InSAR-GNSS(LOS)'] = gps['InSAR'] - gps['los']
        gps.dropna(inplace=True)
        gps_3D = gps[gps['vu'] != 0]

        fig, ax = plt.subplots(1, 4, sharey='all', figsize=(8,3))
        poly.plot(ax=ax[0], facecolor='None', edgecolor='gray')
        poly.plot(ax=ax[1], facecolor='None', edgecolor='gray')
        poly.plot(ax=ax[2], facecolor='None', edgecolor='gray')
        poly.plot(ax=ax[3], facecolor='None', edgecolor='gray')

        gps.plot("InSAR", ax=ax[0], vmin=-10, vmax=10, cmap='bwr')
        ax[0].set_title("InSAR")

        gps.plot("los", ax=ax[1], vmin=-10, vmax=10, cmap='bwr')
        gps_3D.plot(ax=ax[1], facecolor='None', edgecolor='gray')
        ax[1].set_title("GPS LOS")

        gps_3D.plot("InSAR-GNSS(LOS)", ax=ax[2], vmin=-10, vmax=10, cmap='bwr')
        ax[2].set_title("InSAR-3D_GNSS(LOS)")

        im = ax[3].scatter(gps.geometry.x, gps.geometry.y, c=gps["InSAR-GNSS(LOS)"], vmin=-10, vmax=10, cmap='bwr')
        ax[3].set_title("InSAR-GNSS(LOS)")
        cax = inset_axes(
            ax[3],
            width="5%",  # width: 5% of parent_bbox width
            height="100%",  # height: 50%
            loc="lower left",
            bbox_to_anchor=(1.05, 0., 1, 1),
            bbox_transform=ax[3].transAxes,
            borderpad=0,
        )
        fig.colorbar(im, cax=cax)

        plt.suptitle(track)
        fig.savefig("../gps/{}_los_gps.png".format(track), format='PNG', dpi=300, bbox_inches='tight')
        plt.show()

        if len(gps["InSAR-GNSS(LOS)"]) > 0:
            shift_mode = mode(np.around(gps["InSAR-GNSS(LOS)"], decimals=1))[0][0]
        else:
            shift_mode = 0
        export_tif(tif.data - shift_mode, tif, outdir + track + "__shifted.tif")
