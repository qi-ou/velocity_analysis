from merge_tif import *
import os


if __name__ == "__main__":
    # loop through tifs along track and combine two at a time
    #############################

    # Input parameters:
    reference = ''  # can be empty string, _Eurasia or _ITRF2008
    frame_style = '__mode'  # can be empty string, _*range, _*mode, or _*azimuth
    track_style = "mode"   # can add none, range, mode, or azimuth, cannot be an empty string
    input_dir = '../los_tracks_no_plate/'
    output_dir = '../los_tracks_no_plate_merge/'
    input_suffix = frame_style + ".tif"

    ascList = [os.path.join(input_dir, '144A' + input_suffix),
               os.path.join(input_dir, '071A' + input_suffix),
               os.path.join(input_dir, '173A' + input_suffix),
               os.path.join(input_dir, '100A' + input_suffix),
               os.path.join(input_dir, '027A' + input_suffix),
               os.path.join(input_dir, '129A' + input_suffix),
               os.path.join(input_dir, '056A' + input_suffix),
               os.path.join(input_dir, '158A' + input_suffix),
               os.path.join(input_dir, '085A' + input_suffix),
               os.path.join(input_dir, '012A' + input_suffix),
               os.path.join(input_dir, '114A' + input_suffix),
               os.path.join(input_dir, '041A' + input_suffix)]
    dscList = [os.path.join(input_dir, '049D' + input_suffix),
               os.path.join(input_dir, '151D' + input_suffix),
               os.path.join(input_dir, '078D' + input_suffix),
               os.path.join(input_dir, '005D' + input_suffix),
               os.path.join(input_dir, '107D' + input_suffix),
               os.path.join(input_dir, '034D' + input_suffix),
               os.path.join(input_dir, '136D' + input_suffix),
               os.path.join(input_dir, '063D' + input_suffix),
               os.path.join(input_dir, '165D' + input_suffix),
               os.path.join(input_dir, '092D' + input_suffix),
               os.path.join(input_dir, '019D' + input_suffix),
               os.path.join(input_dir, '121D' + input_suffix)]

    for trackList in [dscList ]: #, dscList ascList,
        print(trackList)
        count = 0
        tif1 = trackList[count]
        for count in range(len(trackList)-1):
            print(count)
            tif2 = trackList[count+1]
            # if tif2 == os.path.join(input_dir, '092D' + input_suffix):
            #     merged_basename = combine2tifs(tif1, tif2, track_style, reference, output_dir, rows=750)
            # elif tif2 == os.path.join(input_dir, '019D' + input_suffix):
            #     merged_basename = combine2tifs(tif1, tif2, track_style, reference, output_dir, rows=750)
            # else:
            merged_basename = combine2tifs(tif1, tif2, track_style, reference, output_dir, rows=250)
            tif1 = os.path.join(output_dir, merged_basename)
            count += 1


# 106D-033D-135D-062D-164D_jointly.tif
# 106D-033D-135D-062D-164D_weighted.tif
# 026A-128A-055A_jointly.tif
# 026A-128A-055A_weighted.tif
#
# gdal_calc.py -A 106D-033D-135D-062D-164D_jointly.tif -B ../gps_plots/106D-033D-135D-062D-164D_Eurasia_projectednone.tif --calc=A-B --outfile=jointly-projectednone_dsc.tif
# plot_tif.py -f jointly-projectednone_dsc.tif
# gdal_calc.py -A 026A-128A-055A_jointly.tif -B ../gps_plots/026A-128A-055A_Eurasia_projectednone.tif --calc=A-B --outfile=jointly-projectednone_asc.tif
# plot_tif.py -f jointly-projectednone_asc.tif

# gdal_calc.py -A 106D-033D-135D-062D-164D_Eurasia_range_projected.tif -B ../gps_plots/106D-033D-135D-062D-164D_Eurasia_projectednone.tif --calc=A-B --outfile=optimise-projected_dsc.tif
# plot_tif.py -f optimise-projected_dsc.tif
# gdal_calc.py -A 026A-128A-055A_Eurasia_range_projected.tif -B ../gps_plots/026A-128A-055A_Eurasia_projectednone.tif --calc=A-B --outfile=optimise-projected_asc.tif
# plot_tif.py -f optimise-projected_asc.tif

# gdal_calc.py -A 106D-033D-135D-062D-164D_weighted.tif -B 106D-033D-135D-062D-164D_corner_two_d.tif --calc=A-B --outfile=weighted-corner2d_dsc.tif
# plot_tif.py -f weighted-corner2d_dsc.tif
# gdal_calc.py -A 026A-128A-055A_weighted.tif -B 026A-128A-055A_corner_two_d.tif --calc=A-B --outfile=weighted-corner2d_asc.tif
# plot_tif.py -f weighted-corner2d_asc.tif

# gdal_calc.py -A ../decomposition/final_result/ve.tif -B ve.tif --calc=A-B --outfile=diff_ve.tif
# plot_tif.py -f diff_ve.tif
# gdal_calc.py -A 026A-128A-055A_weighted.tif -B 026A-128A-055A_corner_two_d.tif --calc=A-B --outfile=weighted-corner2d_asc.tif
# plot_tif.py -f weighted-corner2d_asc.tif
