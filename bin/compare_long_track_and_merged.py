import matplotlib.pyplot as plt

from merge_tif import *
input_dir = "/Users/qi/OneDrive - University of Leeds/projects/insar/long_frame_test/long_with_gacos/"

long_vstd = OpenTif(input_dir+"truncate.vstd.tif")
merged_vstd = OpenTif(input_dir+"085A__none.tif")
fig, ax = plt.subplots(1,2, figsize=(6,3))
im = ax[0].imshow(long_vstd.data, vmin=0, vmax=2)
ax[0].set_title("Long track vstd")
plt.colorbar(im, ax=ax[0])
im = ax[1].imshow(merged_vstd.data, vmin=0, vmax=2)
ax[1].set_title("Frame vstd merged")
plt.colorbar(im, ax=ax[1])
plt.show()

long_track = OpenTif(input_dir+"truncate.vel.geo.tif")
merged_track = OpenTif(input_dir+"085A__mode.tif")
overlap = Overlap(long_track, merged_track)
nanmode = plot_hist(overlap.diff_array, "diff_hist.png", "Difference", input_dir)

fig, ax = plt.subplots(1,4, figsize=(13,3))
im = ax[0].imshow(overlap.d1array, vmin=-17.5+nanmode, vmax=2.5+nanmode)
ax[0].set_title("Long track processed")
plt.colorbar(im, ax=ax[0])
im = ax[1].imshow(overlap.d2array, vmin=-17.5, vmax=2.5)
ax[1].set_title("Frames merged")
plt.colorbar(im, ax=ax[1])
im = ax[2].imshow(overlap.diff_array, vmin=-2.5, vmax=2.5)
plt.colorbar(im, ax=ax[2])
ax[2].set_title("Difference")
im = ax[3].imshow(overlap.diff_array - nanmode, vmin=-2.5, vmax=2.5)
plt.colorbar(im, ax=ax[3])
ax[3].set_title("Difference de-mode")
plt.show()

demode = overlap.diff_array - nanmode
nanmode = plot_hist(demode[:1250, :] , "diff_hist.png", "Difference", input_dir)


