import geopandas as gpd
import pandas as pd
from mpl_toolkits.axes_grid1 import ImageGrid
from shapely.geometry import Point
from shapely.geometry import box

from merge_tif import *


def add_NEU_los(gps_df, N_coef_map, E_coef_map, U_coef_map):
    """calculate GPS LOS at gps lon lat from NEU and GPS component velocities"""
    # extract NEU values at gps lat lon from the NEU maps
    gps_df.loc[:, 'n_coef'] = [N_coef_map.extract_pixel_value(point.x, point.y)[0] for point in gps_df['geometry']]
    gps_df.loc[:, 'e_coef'] = [E_coef_map.extract_pixel_value(point.x, point.y)[0] for point in gps_df['geometry']]
    gps_df.loc[:, 'u_coef'] = [U_coef_map.extract_pixel_value(point.x, point.y)[0] for point in gps_df['geometry']]
    # calculate los and los uncertainty per gps point
    gps_df.loc[:, 'los'] = [vn * n + ve * e + vu * u for vn, ve, vu, n, e, u
                            in gps_df[['vn', 've', 'vu', 'n_coef', 'e_coef', 'u_coef']].to_numpy()]
    gps_df.loc[:, 'los_sigma'] = [abs(sn * n) + abs(se * e) + abs(su * u) for sn, se, su, n, e, u
                                  in gps_df[['sn', 'se', 'su', 'n_coef', 'e_coef', 'u_coef']].to_numpy()]


def export_tif(data, df, filename):
    # Export data to tif format.
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(filename, df.xsize, df.ysize, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([df.left, df.xres, 0, df.top, 0, df.yres])  ##sets same geotransform as input
    # outdata.SetProjection(df.projection)  ##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(data)
    outdata.FlushCache()
    outdata.FlushCache()  # need to flush twice to export the last tif properly, otherwise it stops halfway.


def non_nan_merge(big_data, small_data, nodata_test, x_shift, y_shift, xsize, ysize):
    masked_data = np.choose(nodata_test,  # False = not nan = pick from 0th entry; True = nan = pick from 1st entry
                            (small_data, big_data[y_shift:y_shift + ysize, x_shift:x_shift + xsize]))
    big_data[y_shift:y_shift + ysize, x_shift:x_shift + xsize] = masked_data


def load_dummy_gps():
    # gps = gpd.read_file('../gps/dummy.shp')
    gps = gpd.read_file('../gps/dummy3.shp')
    gps['ve'] = 0
    gps['vn'] = 0
    gps['vu'] = 0
    gps['se'] = 1
    gps['sn'] = 1
    gps['su'] = 1
    del gps['id']
    gps.dropna(inplace=True)
    # # try:
    ve = OpenTif('../gps/ve_interpolated.tif')
    se = OpenTif('../gps/se_interpolated.tif')
    vn = OpenTif('../gps/vn_interpolated.tif')
    sn = OpenTif('../gps/sn_interpolated.tif')
    gps['ve'] = [ve.extract_pixel_value(point.x, point.y, 20)[0] for point in gps['geometry']]
    gps['se'] = [se.extract_pixel_value(point.x, point.y, 20)[0] for point in gps['geometry']]
    gps['vn'] = [vn.extract_pixel_value(point.x, point.y, 20)[0] for point in gps['geometry']]
    gps['sn'] = [sn.extract_pixel_value(point.x, point.y, 20)[0] for point in gps['geometry']]
    # # except:
    # #     pass
    return gps


def load_chris():
    """ Load GPS data into a geopandas dataframe
    - only need lon, lat, vel and sig for the purpose of interpolation"""
    gps_gdf = gpd.GeoDataFrame(crs="EPSG:4326")
    gps_gdf['geometry'] = None
    index = 0
    # gps_file = "/Users/qi/leeds/datasets/vectors/chris_tianshan_28_May_2023/tianshan_tol1_minocc2.5_dist6_edges_2D_fake3D_su1.dat"
    # gps_file = "/Users/qi/leeds/softwares/VELMAP/github/velmap/gps/tianshan_tol1_minocc2.5_dist6_edges_3Dandfake3D_su1.dat"
    gps2d_file = "../gps/chris_tianshan_7_july_2023/tianshan_tol2.1_minocc2.5_dist1_edges_2D_7July2023_weighted.dat"
    gps3d_file = "../gps/chris_tianshan_7_july_2023/tianshan_tol2.1_minocc2.5_dist1_edges_3D_7July2023_weighted.dat"

    fl = open(gps2d_file, "r").readlines()
    for line in fl:
        lon, lat, Ve, Vn, dVe, dVn, Cen, sta = line.split()
        gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
        gps_gdf.loc[index, 've'] = float(Ve)  # eastern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vn'] = float(Vn)  # northern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vu'] = 0  # northern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
        gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
        gps_gdf.loc[index, 'su'] = 1  # sigma vn
        index += 1

    fl = open(gps3d_file, "r").readlines()
    for line in fl:
        lon, lat, Ve, Vn, Vu, dVe, dVn, dVu, Cen, Ceu, Cnu, sta = line.split()
        gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
        gps_gdf.loc[index, 've'] = float(Ve)  # eastern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vn'] = float(Vn)  # northern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vu'] = float(Vu)  # northern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
        gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
        gps_gdf.loc[index, 'su'] = float(dVu)  # sigma vn
        index += 1

    # fl = open(gps_file, "r").readlines()
    # for line in fl:
    #     if not line.startswith(('*', "Table")):
    #         lon, lat, Ve_E, Vn_E, Vu_E, dVe, dVn, dVu, Cen, Ceu, Cnu, sta = line.split()
    #         gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
    #         gps_gdf.loc[index, 've'] = float(Ve_E)  # eastern velocity in mm/yr in fixed eurasia reference frame
    #         gps_gdf.loc[index, 'vn'] = float(Vn_E)  # northern velocity in mm/yr in fixed eurasia reference frame
    #         gps_gdf.loc[index, 'vu'] = float(Vu_E)  # northern velocity in mm/yr in fixed eurasia reference frame
    #         gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
    #         gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
    #         gps_gdf.loc[index, 'su'] = float(dVu)  # sigma vn
    #         index += 1

    return gps_gdf


def load_russia():
    """ Load GPS data into a geopandas dataframe
    - only need lon, lat, vel and sig for the purpose of interpolation"""
    gps_gdf = gpd.GeoDataFrame(crs="EPSG:4326")
    gps_gdf['geometry'] = None
    index = 0
    # gps_file = "/Users/qi/leeds/datasets/vectors/GNSS/from_Chris_Rollins/casia_tol1.5_minocc2.5_dist2_2D_dec2022_corner.dat"
    gps_file = "/Users/qi/leeds/datasets/vectors/chris_tianshan_3_july_2023/russia.txt"
    fl = open(gps_file, "r").readlines()
    for line in fl:
        if not line.startswith(('*', "Table")):
            lon, lat, Ve_E, Vn_E, _, dVe, dVn, _, Cen, _, _, sta = line.split()
            gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
            gps_gdf.loc[index, 've'] = float(Ve_E)  # eastern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vn'] = float(Vn_E)  # northern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vu'] = 0  # northern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
            gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
            gps_gdf.loc[index, 'su'] = 1  # sigma vn
            index += 1
    return gps_gdf


def load_tarim():
    """ Load GPS data into a geopandas dataframe
    - only need lon, lat, vel and sig for the purpose of interpolation"""
    gps_gdf = gpd.GeoDataFrame(crs="EPSG:4326")
    gps_gdf['geometry'] = None
    index = 0
    # gps_file = "/Users/qi/leeds/projects/insar/velocity_analysis/gps/tarim.txt"
    gps_file = "/Users/qi/leeds/projects/insar/velocity_analysis/gps/tarim_manual.txt"

    fl = open(gps_file, "r").readlines()
    for line in fl:
        if not line.startswith(('*', "Table")):
            lon, lat, Ve_E, Vn_E, dVe, dVn, _, _, _, _ = line.split()
            gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
            gps_gdf.loc[index, 've'] = float(Ve_E)  # eastern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vn'] = float(Vn_E)  # northern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vu'] = 0  # northern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'se'] = 1  # sigma ve
            gps_gdf.loc[index, 'sn'] = 1  # sigma vn
            gps_gdf.loc[index, 'su'] = 2  # sigma vn
            index += 1

    return gps_gdf


if __name__ == "__main__":
    ###################
    # input parameters:
    ###################
    # reference = '_Eurasia'  # Eurasia or ITRF2008
    input_dir = '../los_full/vel_no_plate/'
    input_suffix = '*no_plate.tif'
    sigma_dir = '../los_full/vstd_corrected/'
    sigma_suffix = '*scaled_vstd.tif'
    NEU_dir = '../NEU/track/'
    N_suffix = '_N.tif'
    E_suffix = '_E.tif'
    U_suffix = '_U.tif'
    output_dir = '../los_full/referenced_by_gps_overlap/'
    output_suffix = '_3'  # _jointly or #_weighted
    frame_3d_gps_threshold = 10

    # set up output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # load gps
    dummy_gps = load_dummy_gps()
    chris_gps = load_chris()
    # russia_gps = load_russia()
    tarim_gps = load_tarim()
    # all_gps = pd.concat([dummy_gps, chris_gps], ignore_index=True) #russia_gps,

    all_gps = pd.concat([dummy_gps, chris_gps, tarim_gps], ignore_index=True)
    # all_gps.plot('ve')
    plt.show()
    # breakpoint()

    # 1. list all frames with the correct suffix in the directory
    tifList = sorted(glob.glob(os.path.join(input_dir, '107*' + input_suffix)))
    print(tifList)

    # 2. define tracks from frame names
    trackList = sorted(set([os.path.basename(t)[:4] for t in tifList]))

    for track in trackList:

        # load insar los maps
        frameList = sorted(glob.glob(os.path.join(input_dir, track + input_suffix)))
        sigfileList = sorted(glob.glob(os.path.join(sigma_dir, track + sigma_suffix)))
        boxList = []
        ds = []
        sigma = []

        for index, frame in enumerate(frameList):
            tif = OpenTif(frame)
            sigtif = OpenTif(sigfileList[index])
            ds.append(tif)
            sigma.append(sigtif)
            boxList.append(box(tif.left, tif.bottom, tif.right, tif.top))

        # define the geographic boundary of the track and create an empty array of the size of the merged track
        poly = gpd.GeoSeries([gpd.GeoSeries(boxList).unary_union])
        [[track_left, track_bottom, track_right, track_top]] = poly.bounds.values
        track_xres = tif.xres
        track_yres = tif.yres
        track_xsize = int((track_right - track_left) / track_xres + 1.5)
        track_ysize = int((track_bottom - track_top) / track_yres + 1.5)

        # load NEU coefficients
        n_coef = OpenTif(os.path.join(NEU_dir, track + N_suffix))
        e_coef = OpenTif(os.path.join(NEU_dir, track + E_suffix))
        u_coef = OpenTif(os.path.join(NEU_dir, track + U_suffix))

        # Step 1: populate the insar gps rows of inversion matrix
        # LOS1-GPS1 = [y x 1 0 0 0 0 0 0][a1 b1 c1 a2 b2 c3 a3 b3 c3]T
        # LOS2-GPS2 = [0 0 0 y x 1 0 0 0][a1 b1 c1 a2 b2 c3 a3 b3 c3]T
        # x and y are col and row in the merged index frame

        # shortlist gps in track
        gps_in_track = all_gps[all_gps.within(poly.loc[0])]

        # load gps insar offsets for inversion
        g_offset = []
        g = []
        sig = []
        ys = []
        xs = []
        gps_los = []
        for index, tif in enumerate(ds):
            frame = tif.basename[:17]
            print(frame)
            # select gps data within the frame

            # to only use 3D gps if there are enough of them.
            gps = gps_in_track[gps_in_track.within(boxList[index])]

            # calculate gps los
            add_NEU_los(gps, n_coef, e_coef, u_coef)

            # extract insar los for gps locations
            gps['insar'] = [ds[index].extract_pixel_value(point.x, point.y, 60)[0] for point in gps['geometry']]
            gps['insar_sigma'] = [ds[index].extract_pixel_value(point.x, point.y, 60)[1] for point in gps['geometry']]
            gps['gps_sigma'] = np.sqrt(
                np.square(gps['insar_sigma'].to_numpy()) + np.square(gps['los_sigma'].to_numpy()))
            gps.dropna(inplace=True)
            gps = gps[gps['insar'] > -3]  # to avoid using gps data in rapidly subsiding fields for referencing
            gps = gps[gps['insar'] < 3]  # to avoid using gps data in rapidly uplifting regions due to earthquakes
            gps = gps[gps['vu'] < 3]  # to avoid using gps data in rapidly uplifting regions due to earthquakes
            gps = gps[gps['vu'] > -3]  # to avoid using gps data in rapidly uplifting regions due to earthquakes

            gps_2D = gps[gps['vu'] == 0]
            gps_3D = gps[gps['vu'] != 0]
            if len(gps_3D) < frame_3d_gps_threshold:
                # if there are too few 3D stations,
                # use 2D stations with Vu = 0, and Su = the largest absolute vu
                gps_2D.loc[:, 'su'] = gps_3D['vu'].abs().max()
                gps = pd.concat([gps_2D, gps_3D])
            else:
                gps = gps_3D
            sig.append(gps['gps_sigma'].to_numpy())
            x = [int((point.x - track_left) / ds[index].xres + 0.5) for point in gps['geometry']]
            y = [int((point.y - track_top) / ds[index].yres + 0.5) for point in gps['geometry']]

            # populate offset
            gps_los.append(gps['los'].to_numpy())
            g_offset.append(gps['insar'].to_numpy() - gps['los'].to_numpy())
            xs.append(x)
            ys.append(y)

            # populate g matrix
            G = np.zeros((gps.shape[0], 3 * len(frameList)))
            G[:, 3 * index + 0] = y
            G[:, 3 * index + 1] = x
            G[:, 3 * index + 2] = 1
            g.append(G)

        all_gps_los = np.concatenate(gps_los)
        all_x = np.concatenate(xs)
        all_y = np.concatenate(ys)
        all_gps_offset = np.concatenate(g_offset)
        all_gps_g = np.concatenate(g)
        all_gps_sig = np.concatenate(sig)
        all_gps_offset_with_sig = all_gps_offset / all_gps_sig
        all_gps_g_with_sig = all_gps_g / all_gps_sig[:, None]  # [:, None] used to allow broadcasting

        # count how many points there are in gps overlaps
        gps_offset_number = len(all_gps_offset)
        all_gps_len_list = list(map(len, g_offset))

        # Step 2: populate the insar overlap rows of inversion matrix
        # LOS1-LOS2 = [y x 1 -y -x -1 0 0 0][a1 b1 c1 a2 b2 c3 a3 b3 c3]T
        # LOS2-LOS3 = [0 0 0 y x 1 -y -x -1][a1 b1 c1 a2 b2 c3 a3 b3 c3]T
        # x and y are col and row in the merged index frame

        if len(frameList) > 1:

            # load insar overlap for inversion
            offset = []
            g = []
            sig = []
            yxs = []
            for i in range(len(frameList) - 1):
                over = Overlap(ds[i], ds[i + 1], vmin=-5, vmax=5)
                sigover = Overlap(sigma[i], sigma[i + 1])
                x_shift = int((over.left - track_left) / ds[i].xres + 0.5)
                y_shift = int((over.top - track_top) / ds[i].yres + 0.5)
                non_nan_yx = np.argwhere(~np.isnan(over.diff_array))
                non_nan_diff = over.diff_array[~np.isnan(over.diff_array)]
                if over.diff_array.shape != sigover.vector_sum_array.shape:
                    (row, col) = over.diff_array.shape
                    sigover.vector_sum_array = sigover.vector_sum_array[:row, :col]  # one case in 135D encountered dimension mismatch
                non_nan_sig = sigover.vector_sum_array[
                    ~np.isnan(over.diff_array)]  # extract the same pixels from sigma overlap
                non_nan_sig[np.isnan(non_nan_sig)] = np.nanmean(non_nan_sig)  # in case of void, fill with nanmean
                offset.append(non_nan_diff)
                sig.append(non_nan_sig)
                non_nan_yx += [y_shift, x_shift]
                yxs.append(non_nan_yx)
                G = np.zeros((non_nan_yx.shape[0], 3 * len(frameList)))
                G[:, 3 * i + 0: 3 * i + 2] = non_nan_yx
                G[:, 3 * i + 3: 3 * i + 5] = -non_nan_yx
                G[:, 3 * i + 2] = 1
                G[:, 3 * i + 5] = -1
                g.append(G)
                # overlap_map[non_nan_yx[:,0], non_nan_yx[:,1]] = non_nan_diff  # add for plotting

            all_overlap_offset = np.concatenate(offset)
            all_overlap_g = np.concatenate(g)
            all_overlap_sig = np.concatenate(sig)
            all_overlap_yxs = np.concatenate(yxs)
            # count how many point there are in insar overlaps
            insar_offset_number = len(all_overlap_offset)
            insar_offset_len_list = list(map(len, offset))

            # # Weighting
            weight = np.sqrt(insar_offset_number / gps_offset_number)
            all_overlap_offset_with_sig = all_overlap_offset / all_overlap_sig
            all_overlap_g_with_sig = all_overlap_g / all_overlap_sig[:, None]  # [:, None] used to allow broadcasting
            all_offset = np.concatenate([all_overlap_offset_with_sig, weight * all_gps_offset_with_sig])
            all_g = np.concatenate([all_overlap_g_with_sig, weight * all_gps_g_with_sig])
        else:
            weight = 1
            all_offset = all_gps_offset_with_sig
            all_g = all_gps_g_with_sig

        # inversion
        coefs, res, rank, singular = np.linalg.lstsq(all_g, all_offset, rcond=None)
        # calc residuals
        model = np.dot(all_g, coefs)
        all_res = all_offset - model
        all_gps_res = all_res[-gps_offset_number:] * all_gps_sig / weight

        # calculate residuals
        if len(frameList) > 1:
            all_insar_res = all_res[:insar_offset_number] * all_overlap_sig
        else:
            all_insar_res = []

        insar_offset_number = len(all_overlap_offset)
        insar_offset_len_list = list(map(len, offset))

        # populate the track with los map arrays for plotting
        track_map = np.ones((track_ysize, track_xsize)) * np.nan
        # populate the track with los map arrays for plotting
        track_map_proj = np.ones((track_ysize, track_xsize)) * np.nan
        ramps_map = np.ones((track_ysize, track_xsize)) * np.nan

        # reconstruct the full ramp image and export the projected insar los velocity map
        for i in np.arange(len(frameList)):
            # reconstruct ramp
            full_ramp_array = np.zeros(ds[i].data.shape)
            all_yx = np.argwhere(~np.isnan(full_ramp_array))
            x_shift = int((ds[i].left - track_left) / ds[i].xres + 0.5)
            y_shift = int((ds[i].top - track_top) / ds[i].yres + 0.5)
            all_yx += [y_shift, x_shift]
            G = np.ones((len(all_yx), 3))
            G[:, :2] = all_yx
            full_ramp_array = np.dot(G, coefs[3 * i:3 * i + 3]).reshape(ds[i].ysize, ds[i].xsize)

            # calc projected insar
            ds[i].data_projected = ds[i].data - full_ramp_array

            # populate the new track map for plotting
            nodata_test = np.isnan(ds[i].data_projected)
            non_nan_merge(track_map, ds[i].data, nodata_test, x_shift, y_shift, ds[i].xsize, ds[i].ysize)
            non_nan_merge(track_map_proj, ds[i].data_projected, nodata_test, x_shift, y_shift, ds[i].xsize, ds[i].ysize)
            non_nan_merge(ramps_map, full_ramp_array, nodata_test, x_shift, y_shift, ds[i].xsize, ds[i].ysize)

            # Export projected insar to tif format.
            export_tif(ds[i].data_projected, ds[i], output_dir + ds[i].basename[:17] + output_suffix + '.tif')

        # start a 2x4 figure for displays
        fig = plt.figure(figsize=(7, 5.6))
        axs = ImageGrid(fig, 111,
                        nrows_ncols=(2, 3),
                        axes_pad=0.2,
                        cbar_mode='single',
                        cbar_location='left',
                        cbar_pad=-0.1,
                        cbar_size="4%",
                        share_all=True
                        )

        # remove the x and y ticks
        for ax in axs:
            ax.set_xticks([])
            ax.set_yticks([])

        vmin = np.nanpercentile(track_map, 0.5)
        vmax = np.nanpercentile(all_gps_offset, 99.5)
        vmin = max(vmin, -10)
        # vmin = -2
        # vmax = 2
        im = axs[0].imshow(track_map, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        axs[0].set_title("InSAR LOS", pad=2)
        cbar = ax.cax.colorbar(im)
        cbar.set_label_text("LOS velocity / mm per yr")

        # populate the gps map with empty array followed by gps overlaps
        gps_map = np.ones((track_ysize, track_xsize)) * np.nan
        axs[1].imshow(gps_map)  # empty map to define the dimension
        axs[1].scatter(all_x, all_y, c=all_gps_offset, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        axs[1].set_title("GNSS Offset", pad=2)

        axs[3].imshow(ramps_map, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        axs[3].set_title("Inverted Ramps", pad=2)

        axs[4].imshow(gps_map)  # empty map to define the dimension
        axs[4].scatter(all_x, all_y, c=all_gps_res, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        axs[4].set_title("GNSS Residual", pad=2)

        if len(frameList) > 1:
            # plot raw los maps as a track
            overlap_map = np.ones((track_ysize, track_xsize)) * np.nan
            overlap_map[all_overlap_yxs[:, 0], all_overlap_yxs[:, 1]] = all_overlap_offset
            axs[2].imshow(overlap_map, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
            axs[2].set_title("Overlap Offset", pad=2)

            overlap_res_map = np.ones((track_ysize, track_xsize)) * np.nan
            overlap_res_map[all_overlap_yxs[:, 0], all_overlap_yxs[:, 1]] = all_insar_res
            axs[5].imshow(overlap_res_map, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
            axs[5].set_title("Overlap Residual", pad=2)

        plt.suptitle(track)
        plt.show()
        fig.savefig(os.path.join(output_dir, track + "_stitch.png"), format='PNG', dpi=300, bbox_inches='tight')

        # plot histograms
        number_of_hists = 2 * len(
            frameList)  # for each track, there are gps before and after, and overlaps before and after
        fig, axs = plt.subplots(number_of_hists, 1, sharex=True, figsize=(2, 5.6))
        for ax in axs:
            ax.set_yticks([])
            ax.set_xlim([-10, 10])

        if len(frameList) == 5:
            data_list = [all_gps_offset,
                         offset[0],
                         offset[1],
                         offset[2],
                         offset[3],
                         all_gps_res,
                         all_insar_res[:insar_offset_len_list[0]],
                         all_insar_res[insar_offset_len_list[0]:insar_offset_len_list[0] + insar_offset_len_list[1]],
                         all_insar_res[insar_offset_len_list[0] + insar_offset_len_list[1]:insar_offset_len_list[0] +
                                                                                           insar_offset_len_list[1] +
                                                                                           insar_offset_len_list[2]],
                         all_insar_res[insar_offset_len_list[0] + insar_offset_len_list[1] + insar_offset_len_list[2]:
                                       insar_offset_len_list[0] + insar_offset_len_list[1] + insar_offset_len_list[
                                           2] + + insar_offset_len_list[3]]
                         ]

        if len(frameList) == 4:
            data_list = [all_gps_offset,
                         offset[0],
                         offset[1],
                         offset[2],
                         all_gps_res,
                         all_insar_res[:insar_offset_len_list[0]],
                         all_insar_res[insar_offset_len_list[0]:insar_offset_len_list[0] + insar_offset_len_list[1]],
                         all_insar_res[insar_offset_len_list[0] + insar_offset_len_list[1]:insar_offset_len_list[0] +
                                                                                           insar_offset_len_list[1] +
                                                                                           insar_offset_len_list[2]]
                         ]

        if len(frameList) == 3:
            data_list = [all_gps_offset,
                         offset[0],
                         offset[1],
                         all_gps_res,
                         all_insar_res[:insar_offset_len_list[0]],
                         all_insar_res[insar_offset_len_list[0]:insar_offset_len_list[0] + insar_offset_len_list[1]]
                         ]

        if len(frameList) == 2:
            data_list = [all_gps_offset, offset[0],
                         all_gps_res, all_insar_res[:insar_offset_len_list[0]]]

        if len(frameList) == 1:
            data_list = [all_gps_offset, all_gps_res]

        for i, ax in enumerate(axs):
            data = data_list[i]
            sns.histplot(data, kde=True, ax=ax)
            try:
                nanmode = mode(np.around(data[~np.isnan(data)], decimals=1))[0][0]
            except:
                nanmode = 0
            ax.annotate("m: %.1f \n s: %.1f " % (nanmode, np.nanstd(data)),
                        xy=(0.03, 0.93), xycoords='axes fraction',
                        ha="left", va="top",
                        fontsize=12)
        ax.set_xlabel("mm/yr")
        axs[0].set_title("Histograms", pad=2)
        plt.show()
        fig.savefig(os.path.join(output_dir, track + "_stitch_hist.png"), format='PNG', dpi=300, bbox_inches='tight')

        vmin = np.nanpercentile(track_map_proj, 0.5)
        vmin = max(vmin, -10)
        vmax = np.nanpercentile(track_map_proj, 99.5)
        plt.imshow(track_map_proj, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        plt.show()
