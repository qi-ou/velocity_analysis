from merge_tif import *

if __name__ == "__main__":

    reference = ''  # like a prefix in the input filenames
    input_dir = '../los_weighted/vel_no_plate/'  # where input velocity frames are stored
    input_suffix = '.vel.no_plate.tif'

    # 1. list all frames with the correct suffix in the directory
    tifList = sorted(glob.glob(os.path.join(input_dir, '151D_04454_121212*'+input_suffix)))
    print(tifList)

    # 2. define tracks from frame names
    trackList = set([os.path.basename(t)[:5] for t in tifList])
    for track in trackList:
        frameList = sorted(glob.glob(os.path.join(input_dir, track+'*'+input_suffix)))
        north_frame = frameList[0]
        tif = OpenTif(north_frame)
        std = np.nanstd(tif.data)
        mode = plot_hist(tif.data, input_dir+tif.basename[:-4]+"hist.png", tif.basename[:-4])
        export_tif(tif.data - mode, tif, north_frame)
        # with open(input_dir+'{}.txt'.format(track), 'w') as f:
        #     print(std, file=f)