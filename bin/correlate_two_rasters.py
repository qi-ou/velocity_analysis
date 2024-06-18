from merge_tif import *

tif1 = OpenTif('../../glacier/rate_merge.tif')
tif2 =  OpenTif('../decompose/vu.tif')
overlap = Overlap(tif1, tif2)

d1 = overlap.d1array.flatten()*1000
d2 = overlap.d2array.flatten()

# plt.scatter(d1, d2)
vmin1 = np.nanpercentile(d1,5)
vmax1 = np.nanpercentile(d1,95)
plt.hist(d1[np.logical_and(d1>vmin1, d1<vmax1)], bins=100)
plt.title("glacier mm/y")
plt.show()
vmin2 = np.nanpercentile(d2,5)
vmax2 = np.nanpercentile(d2,95)
plt.hist(d2[np.logical_and(d2>vmin2, d2<vmax2)], bins=100)
plt.title("insar vu mm/y")
plt.show()

mask_d1 = np.logical_and(d1>vmin1, d1<vmax1)
mask_d2 = np.logical_and(d2>vmin2, d2<vmax2)

mask = np.logical_and(mask_d1, mask_d2)
plt.scatter(d1[mask], d2[mask])
plt.hist2d(d1[mask], d2[mask], bins=100)
plt.xlabel("glacier vu, mm/yr")
plt.ylabel("insar vu, mm/y")
plt.show()

