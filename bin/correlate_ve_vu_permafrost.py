from merge_tif import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy
out_dir = '../los_full/decompose/'

tif1 = OpenTif('../los_full/decompose/ve_3_remove_median_filter_100km.tif')
tif2 = OpenTif('/Users/qi/leeds/datasets/rasters/tianshan_permafrost_0.005.tif')
overlap = Overlap(tif1, tif2)

fig, ax = plt.subplots(1, 3, sharey='all')
vmin = np.nanpercentile(overlap.d1array, 1)
vmax = np.nanpercentile(overlap.d1array, 99)
overlap.d2array[overlap.d2array<0]=np.nan
ax[0].imshow(overlap.d1array, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
ax[1].imshow(overlap.d2array, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
ve_under_permafrost = copy.copy(overlap.d1array)
ve_under_permafrost[np.isnan(overlap.d2array)] = np.nan
im = ax[2].imshow(ve_under_permafrost, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
plt.colorbar(im, ax=ax, orientation='horizontal')
#
# vmin = np.nanpercentile(overlap.diff_array, 1)
# vmax = np.nanpercentile(overlap.diff_array, 99)
# im = ax[2].imshow(overlap.diff_array, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
# plt.colorbar(im, ax=ax[2], orientation='horizontal')
# ax[0].set_title(dir1)
# ax[1].set_title(dir2)
# ax[2].set_title('Diff Panel 1-2')
# plt.suptitle(filename)
# fig.savefig(out_dir+filename[:-4]+'.png', format='PNG', dpi=300, bbox_inches='tight', transparent=True)
# plt.tight_layout()
plt.show()


fig, ax = plt.subplots(figsize=(10,5))
im=ax.imshow(ve_under_permafrost, cmap=cm.roma.reversed(), vmin=vmin, vmax=vmax)
plt.yticks([])
plt.xticks([])
ax.set_title("ve_under_permafrost")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)
c = plt.colorbar(im, ax=ax, cax=cax)
plt.show()
fig.savefig(out_dir+"ve_under_permafrost.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)


fig, ax = plt.subplots(figsize=(3.5,2.5))
plt.yticks([])
data = ve_under_permafrost.flatten()
plt.hist(data, bins=np.arange(-10,10,0.01))
nanmode = mode(np.around(data, decimals=1))[0][0]
plt.annotate("mode = %.1f \n std = %.1f " % (nanmode, np.nanstd(data)),
                 xy=(0.985, 0.965), xycoords='axes fraction',
                 ha="right", va="top",
                 fontsize=14)
plt.show()
fig.savefig(out_dir+"ve_under_permafrost_hist.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)


# Export merged data to tif format.
driver = gdal.GetDriverByName("GTiff")
outdata = driver.Create(out_dir+'ve_under_permafrost.tif', ve_under_permafrost.shape[1], ve_under_permafrost.shape[0], 1, gdal.GDT_Float32)
outdata.SetGeoTransform([overlap.left, 0.005, 0, overlap.top, 0, -0.005])  ##sets same geotransform as input
outdata.SetProjection(overlap.d1.projection)  ##sets same projection as input
outdata.GetRasterBand(1).WriteArray(ve_under_permafrost)

