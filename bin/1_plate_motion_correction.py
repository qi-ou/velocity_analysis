from merge_tif import *

if __name__ == "__main__":
    #############################
    # Input parameters:
    style = '_mode'  # can be empty string, _*range, _*mode, or _*azimuth _range

    reference = ''  # like a prefix in the input filenames
    input_dir = '../los_weighted/vel/'  # where input velocity frames are stored
    input_suffix = '129A_04939_131313.weighted.vel.tif'
    out_dir = '../los_weighted/vel_no_plate/'
    NEU_dir = '../NEU/'
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    ###############################
    # 0. load ITRF plate motion
    plate_north = OpenTif("../plate_motion/plate_motion_north.tif")
    plate_east = OpenTif("../plate_motion/plate_motion_east.tif")

    # 1. list all frames with the correct suffix in the directory
    tifList = glob.glob(os.path.join(input_dir, '*'+input_suffix))

    for file in tifList:
        frame = os.path.splitext(os.path.basename(file))[0][:-13]
        tif = OpenTif(file,
                      N="../NEU/"+frame+".geo.N.ml10.tif",
                      E="../NEU/"+frame+".geo.E.ml10.tif")
        outfile = os.path.join(out_dir, frame+'.vel.no_plate.tif')
        outpng = os.path.join(out_dir, frame+'.vel.no_plate.png')
        plate_motion_file = os.path.join(out_dir, frame+'.vel.plate_motion.tif')

        north = Overlap(plate_north, tif).d1array
        east = Overlap(plate_east, tif).d1array

        # %% Set preliminaly reference
        reffile = os.path.join('../refs/'+frame+'_ref.txt')
        with open(reffile, "r") as f:
            refarea = f.read().split()[0]  # str, x1/x2/y1/y2
        refx1, refx2, refy1, refy2 = [int(s) for s in re.split('[:/]', refarea)]

        try:
            plate_motion = north * tif.N + east * tif.E
            ref = np.nanmean(plate_motion[refy1: refy2, refx1: refx2])
            plate_motion-=ref
            fig, ax = plt.subplots(1,3,sharey='all', sharex='all', figsize=(6.4, 2.8))
            vmin=np.nanpercentile(tif.data, 1)
            vmax = np.nanpercentile(tif.data, 99)
            im0 = ax[0].imshow(tif.data, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
            # vmin=np.nanpercentile(plate_motion, 1)
            # vmax = np.nanpercentile(plate_motion, 99)
            im1 = ax[1].imshow(plate_motion, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
            corrected = tif.data - plate_motion
            # vmin=np.nanpercentile(corrected, 1)
            # vmax = np.nanpercentile(corrected, 99)
            im2 = ax[2].imshow(corrected, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
            plt.colorbar(im0, ax=ax, aspect=40, orientation="horizontal")
            # plt.colorbar(im1, ax=ax[1], orientation="horizontal")
            # plt.colorbar(im2, ax=ax[2], orientation="horizontal")
            ax[0].set_title("los, std={:.1f}".format(np.nanstd(tif.data.flatten())))
            ax[1].set_title("plate_motion, std={:.1f}".format(np.nanstd(plate_motion.flatten())))
            ax[2].set_title("corrected, std={:.1f}".format(np.nanstd(corrected.flatten())))
            ax[0].plot([refx1, refx2, refx2, refx1, refx1], [refy1, refy1, refy2, refy2, refy1], color='red', linewidth=0.5)#, fillstyle="none", markerfacecolor='none')
            ax[1].plot([refx1, refx2, refx2, refx1, refx1], [refy1, refy1, refy2, refy2, refy1], color='red', linewidth=0.5)#, fillstyle="none", markerfacecolor='none')
            ax[2].plot([refx1, refx2, refx2, refx1, refx1], [refy1, refy1, refy2, refy2, refy1], color='red', linewidth=0.5)#, fillstyle="none", markerfacecolor='none')
            plt.suptitle(frame)
            plt.savefig(outpng, format='PNG', dpi=300, bbox_inches='tight', transparent=True)
            plt.show()

            export_tif(corrected, tif, outfile)
            export_tif(plate_motion, tif, plate_motion_file)
        except:
            print("an exception has occurred at frame{}".format(frame))