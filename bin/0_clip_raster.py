from merge_tif import *
import copy
# input_list = "/Users/qi/leeds/projects/insar/velocity_analysis/corrected_vstd/056A_04414_191919_lower.fringe.vstd_scaled_vstd.tif"
# output_dir = "/Users/qi/leeds/projects/insar/velocity_analysis/corrected_vstd/"
# #
# input_list = "../los_weighted/vel_no_plate/092D_04475_131313.vel.no_plate.tif"
# output_dir = "../los_weighted/vel_no_plate/"
#
# tifList = glob.glob(input_list)
# print(tifList)
# for file in tifList:
#     tif=OpenTif(file)
#     pix_lin, pix_col = np.indices((tif.ds.RasterYSize, tif.ds.RasterXSize))
#
#     data = copy.copy(tif.data)
#     data[420:] = np.nan
#     # data[data < -2]=np.nan
#     # data[data > 1] = np.nan
#
#     fig, ax = plt.subplots(1,2)
#     vmin = np.nanpercentile(data, 1)
#     vmax = np.nanpercentile(data, 99)
#     im = ax[0].imshow(tif.data, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
#     ax[1].imshow(data, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
#     plt.colorbar(im, ax=ax, label='LOS, mm/yr', orientation='bottom')
#     plt.suptitle(tif.basename)
#     plt.show()
#
#     export_tif(data, tif, file[:-4]+"_masked.tif")
#
# #
# # mask tif1 with overlapping pixels in tif2
# tif1 = OpenTif("../los_weighted/vel_no_plate/092D_04475_131313.vel.no_plate.tif")
# tif2 = OpenTif("../los_weighted/vel_no_plate/092D_04679_141413.vel.no_plate.tif")
# overlap = Overlap(tif1, tif2)
# d2nonnan_mask = ~np.isnan(overlap.d2array)
# overlap.d1.data[overlap.d1yt: overlap.d1yb, overlap.d1xl: overlap.d1xr][d2nonnan_mask] = np.nan
# export_tif(overlap.d1.data, tif1, "../los_weighted/vel_no_plate/092D_04475_131313.vel.no_plate_masked.tif")
#
#
#
# input_list = "../los_full/vel_no_plate/056A_04661_131313.vel.no_plate.tif"
# output_dir = "../los_full/vel_no_plate/vel_no_plate/"
#
# tifList = glob.glob(input_list)
# print(tifList)
# for file in tifList:
#     tif=OpenTif(file)
#     data = copy.copy(tif.data)
#     data[490:] = np.nan
#     # data[data < -2]=np.nan
#     # data[data > 1] = np.nan
#
#     fig, ax = plt.subplots(1,2)
#     vmin = np.nanpercentile(data, 1)
#     vmax = np.nanpercentile(data, 99)
#     im = ax[0].imshow(tif.data, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
#     ax[1].imshow(data, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
#     plt.colorbar(im, ax=ax, label='LOS, mm/yr', orientation='bottom')
#     plt.suptitle(tif.basename)
#     plt.show()
#
#     export_tif(data, tif, file[:-4]+"_masked.tif")
#
# #
# #
# mask_file = "../los_weighted/vstd_corrected/107D_04294_131313.weighted.vstd_scaled_vstd.tif"
# mask_dir = "../mask_ts/"
# tif_file = '../los_weighted/referenced_by_gps_overlap/107D_04294_131313_improved_few_tarim.tif'
# tif = OpenTif(tif_file)
# # mask = OpenTif(mask_file)
# # fig, ax = plt.subplots(1, 2)
# # plt.imshow(mask.data)
# # plt.colorbar()
# # plt.show()
#
# fig, ax = plt.subplots(1, 3)
# ax[0].imshow(tif.data)
# # ax[1].imshow(mask.data)
# tif.data[:350] = np.nan
# ax[2].imshow(tif.data)
# plt.show()
#
# export_tif(tif.data, tif, tif_file)
#
# # tif.data[mask.data>]=np.nan
# export_tif(tif.data, tif, '../los_weighted/referenced_by_gps_overlap/107D_04294_131313_masked.tif')
#
# #
#
# vel = OpenTif('los_weighted/vel/056A_04414_191919_upper.weighted.vel.tif')
# for i in "NEU":
#     infile = OpenTif('NEU/056A_04414_191919.geo.{}.ml10.tif'.format(i))
#     outfile = 'NEU/056A_04414_191919_upper.geo.{}.ml10.tif'.format(i)
#     overlap = Overlap(vel, infile)
#     export_tif(overlap.d2array, vel, outfile)
#
#
#
# tif_file = '../los_improved_weighted/vstd/129A_04474_212121.weighted.vstd.tif'
# mask_file = '../los_improved_weighted/vel/129A_04474_212121_weighted.vel.tif'
#
# tif = OpenTif(tif_file)
# mask = OpenTif(mask_file)
# tif.data[np.isnan(mask.data)]=np.nan
# export_tif(tif.data, tif, tif_file)


#
# input_list = "../los_weighted/referenced_by_gps_overlap/*_improved_few_tarim.tif"
# mask_dir = "../mask_ts/"
#
# tifList = glob.glob(input_list)
# print(tifList)
# for file in tifList:
#     try:
#         tif = OpenTif(file)
#         frame = tif.basename[:17]
#         mask = OpenTif(mask_dir+frame+'_mask.tif')
#         print(frame)
#         fig, ax = plt.subplots(1,3)
#         ax[0].imshow(tif.data)
#         ax[1].imshow(mask.data)
#         tif.data[np.isnan(mask.data)] = np.nan
#         ax[2].imshow(tif.data)
#         plt.suptitle(frame)
#         plt.show()
#         export_tif(tif.data, tif, '../los_weighted/referenced_by_gps_overlap/'+frame+'_masked.tif')
#     except:
#         print(frame+" cannot be masked properly")
#         fig, ax = plt.subplots(1,2)
#         ax[0].imshow(tif.data)
#         ax[1].imshow(mask.data)
#         plt.title(frame)
#         plt.show()
#         export_tif(tif.data, tif, '../los_weighted/referenced_by_gps_overlap/' + frame + '_masked.tif')
#
#
# infile = "../los_weighted/decompose/ve_masked.tif"
# tif = OpenTif(infile)
# fig, ax = plt.subplots(2,1)
# ax[0].imshow(tif.data, vmin=-2.5, vmax=3.5, cmap=cm.roma.reversed())
# tif.data[tif.data>4] = np.nan
# tif.data[tif.data<-2.5] = np.nan
# ax[1].imshow(tif.data, vmin=-2.5, vmax=3.5, cmap=cm.roma.reversed())
# plt.show()
# export_tif(tif.data, tif, '../los_weighted/decompose/ve_masked_clipped.tif')
#
#
# infile = "../los_weighted/decompose/ve_unmasked.tif"
# tif = OpenTif(infile)
# fig, ax = plt.subplots(2,1)
# ax[0].imshow(tif.data, vmin=-2.5, vmax=3.5, cmap=cm.roma.reversed())
# tif.data[tif.data>4] = np.nan
# tif.data[tif.data<-2.5] = np.nan
# ax[1].imshow(tif.data, vmin=-2.5, vmax=3.5, cmap=cm.roma.reversed())
# plt.show()
# export_tif(tif.data, tif, '../los_weighted/decompose/ve_unmasked_clipped.tif')
#


#
# input_list = "../los_full/vel/063D_04245_121111.full.weighted.vel.tif"
# mask_dir = "../los_full/vel/"
# output_dir = "../los_full/vel_final/"
#
# tifList = glob.glob(input_list)
# print(tifList)
# for file in tifList:
#     # try:
#         tif = OpenTif(file)
#         frame = tif.basename[:17]
#         mask = OpenTif(mask_dir+frame+'.full.vel.mskd.tif')
#         print(frame)
#         fig, ax = plt.subplots(1,3)
#         ax[0].imshow(tif.data)
#         ax[1].imshow(mask.data)
#         tif.data[np.isnan(mask.data)] = np.nan
#         ax[2].imshow(tif.data)
#         plt.suptitle(frame)
#         plt.savefig(output_dir+frame+'_mask.png')
#         plt.show()
#         export_tif(tif.data, tif, output_dir+frame+'_masked.tif')
    # except:
    #     print(frame+" cannot be masked properly")
    #     fig, ax = plt.subplots(1,2)
    #     ax[0].imshow(tif.data)
    #     ax[1].imshow(mask.data)
    #     plt.title(frame)
    #     plt.show()
    #     # export_tif(tif.data, tif, '../los_weighted/referenced_by_gps_overlap/' + frame + '_masked.tif')
# #
# # 107D_04294_131313 cannot be masked properly
# # 056A_04414_191919 cannot be masked properly
#

#
input_list = "../los_full/vel/063D_04245_121111*vstd.tif"
mask_dir = "../los_full/vel_final/"
output_dir = "../los_full/vstd_final/"

tifList = glob.glob(input_list)
print(tifList)
for file in tifList:
    try:
        tif = OpenTif(file)
        frame = tif.basename[:17]
        mask = OpenTif(mask_dir+frame+'.full.weighted.vel.tif')
        print(frame)
        fig, ax = plt.subplots(1,3)
        ax[0].imshow(tif.data)
        ax[1].imshow(mask.data)
        tif.data[np.isnan(mask.data)] = np.nan
        ax[2].imshow(tif.data)
        plt.suptitle(frame)
        plt.savefig(output_dir+frame+'_mask.png')
        plt.show()
        export_tif(tif.data, tif, output_dir+frame+'.vstd.masked.tif')
    except:
        print(frame+" cannot be masked properly")
        fig, ax = plt.subplots(1,2)
        ax[0].imshow(tif.data)
        ax[1].imshow(mask.data)
        plt.title(frame)
        plt.show()
        # export_tif(tif.data, tif, '../los_weighted/referenced_by_gps_overlap/' + frame + '_masked.tif')


tif="/Users/qi/leeds/datasets/rasters/tianshan_permafrost_0.005.tif"
nan="/Users/qi/leeds/datasets/rasters/tianshan_permafrost_0.005_nan.tif"
df=OpenTif(tif)
df.data[df.data<0]=np.nan
export_tif(df.data, df, nan)
