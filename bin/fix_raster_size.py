from merge_tif import *

file="../los_tracks_no_plate/092D__mode.tif"
final_xsize = 1220
final_ysize = 2248

tif = OpenTif(file)
final_data = tif.data[:final_ysize, :final_xsize]
# Export data to tif format.
driver = gdal.GetDriverByName("GTiff")
outdata = driver.Create(file, final_xsize, final_ysize, 1, gdal.GDT_Float32)
outdata.SetGeoTransform([tif.left, tif.xres, 0, tif.top, 0, tif.yres])  ##sets same geotransform as input
# outdata.SetProjection(df.projection)  ##sets same projection as input
outdata.GetRasterBand(1).WriteArray(final_data)
outdata.FlushCache()
outdata.FlushCache()  # need to flush twice to export the last tif properly, otherwise it stops halfway.

