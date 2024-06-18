from merge_tif import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == "__main__":
    dir1 = '../los_full/decompose/'
    dir2 = '../los_weighted/decompose/'
    out_dir = '../los_full/decompose/'
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    tifList = glob.glob(os.path.join(dir1, 'vu.tif'))
    for file in tifList:
        filename = os.path.basename(file)
        # tif1 = OpenTif(file)
        # tif2 = OpenTif(os.path.join(dir2, filename))
        tif1 = OpenTif('../los_full/decompose/ve_3.tif')
        tif2 = OpenTif('/Users/qi/leeds/datasets/rasters/tianshan_permafrost.tif')
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
        # ax[0].set_title(dir1)
        # ax[1].set_title(dir2)
        # ax[2].set_title('Diff Panel 1-2')
        # plt.suptitle(filename)
        # fig.savefig(out_dir+filename[:-4]+'.png', format='PNG', dpi=300, bbox_inches='tight', transparent=True)
        # plt.tight_layout()
        plt.show()


# pw=`pwd`
# for i in * ; do
#     cd $pw/$i
#     LiCSBAS132_3D_correction.py -o GEOCml10GACOS4 --suffix 3 -s 0.15
# done
fig, ax = plt.subplots(figsize=(10,5))
im=ax.imshow(overlap.diff_array, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
plt.yticks([])
plt.xticks([])
ax.set_title("full_3-full")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)
c = plt.colorbar(im, ax=ax, cax=cax)
plt.show()
fig.savefig(out_dir+"full_3-full.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)


fig, ax = plt.subplots(figsize=(3.5,2.5))
plt.yticks([])
data = overlap.diff_array.flatten()
plt.hist(data, bins=np.arange(vmin,vmax,0.01))
nanmode = mode(np.around(data, decimals=2))[0][0]
plt.annotate("mode = %.2f \n std = %.2f " % (nanmode, np.nanstd(data)),
                 xy=(0.985, 0.965), xycoords='axes fraction',
                 ha="right", va="top",
                 fontsize=14)
plt.show()
fig.savefig(out_dir+"full_3-full_hist.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)


export_tif(overlap.diff_array, tif1, "../los_full/decompose/full-dt13.tif")