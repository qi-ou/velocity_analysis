from merge_tif import *

if __name__ == "__main__":
    dir1 = '../los_weighted/vstd/'
    dir2 = '../los_weighted/'
    out_dir = '../los_weighted/'
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    tifList = glob.glob(os.path.join(dir1, '041A_04587_131313.weighted.vstd.tif'))
    for file in tifList:
        filename = os.path.basename(file)
        tif1 = OpenTif(file)
        tif2 = OpenTif(os.path.join(dir2, filename))
        overlap = Overlap(tif1, tif2)
        fig, ax = plt.subplots(1, 3, sharey='all')
        vmin=np.nanpercentile(overlap.d1array, 1)
        vmax=np.nanpercentile(overlap.d1array, 99)
        ax[0].imshow(overlap.d1array, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        im = ax[1].imshow(overlap.d2array, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax[:2], orientation='horizontal')
        vmin=np.nanpercentile(overlap.diff_array, 1)
        vmax=np.nanpercentile(overlap.diff_array, 99)
        im = ax[2].imshow(overlap.diff_array, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax[2], orientation='horizontal')
        ax[0].set_title(dir1)
        ax[1].set_title(dir2)
        ax[2].set_title('Diff Panel 1-2')
        plt.suptitle(filename)
        fig.savefig(out_dir+filename[:-4]+'.png', format='PNG', dpi=300, bbox_inches='tight', transparent=True)
        plt.show()


# pw=`pwd`
# for i in * ; do
#     cd $pw/$i
#     LiCSBAS132_3D_correction.py -o GEOCml10GACOS4 --suffix 3 -s 0.15
# done
