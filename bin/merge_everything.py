from merge_tif import *
from 4_decompose_into_ve_vu_using_interpolated_vn_and_NEU import *

if __name__ == "__main__":
    #############################
    # Input parameters:
    input_dir = '../../glacier/avg/'  # where input velocity frames are stored
    output_dir = '../../glacier/avg/'
    input_suffix = '_avg.tif'
    output_suffix = 'rate_merge'
    # output_suffix = 'err_merge'


    # 1. list all frames with the correct suffix in the directory
    tifList = glob.glob(os.path.join(input_dir, '*'+input_suffix))
    print(tifList)

    # load data
    df_list = []
    for tif in tifList:
        df_list.append(OpenTif(tif))

    # define the geographic boundary of the canvas
    left = []
    right = []
    top = []
    bottom = []
    xres = []
    yres = []

    for df in df_list:
        left.append(df.left)
        right.append(df.right)
        top.append(df.top)
        bottom.append(df.bottom)
        xres.append(df.xres)
        yres.append(df.yres)

    logger.info("define the geographic boundary of the map")
    canvas = Raster(north=max(top), south=min(bottom), west=min(left), east=max(right), x_step=xres[0], y_step=yres[0])
    merge_map = np.ones((canvas.ysize, canvas.xsize)) * np.nan

    # merge
    for df in df_list:
        x_shift = int((df.left - canvas.left) / canvas.xres + 0.5)
        y_shift = int((df.top - canvas.top) / canvas.yres + 0.5)
        non_nan_merge(merge_map, df.data, np.isnan(df.data), x_shift, y_shift, df.xsize, df.ysize)

    export_tif(merge_map, canvas, output_dir + '{}.tif'.format(output_suffix))
