import glob
import os
import re
import gdal
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import mode
gdal.UseExceptions()
# import matplotlib.cm as cm
from cmcrameri import cm
import copy
import shutil
import gdal
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)

class OpenTif:
    """ a Class that stores the band array and metadata of a Gtiff file."""
    def __init__(self, filename, sigfile=None, incidence=None, heading=None, N=None, E=None, U=None):
        self.ds = gdal.Open(filename)
        self.basename = os.path.splitext(os.path.basename(filename))[0]
        self.band = self.ds.GetRasterBand(1)
        self.data = self.band.ReadAsArray()
        self.xsize = self.ds.RasterXSize
        self.ysize = self.ds.RasterYSize
        self.left = self.ds.GetGeoTransform()[0]
        self.top = self.ds.GetGeoTransform()[3]
        self.xres = self.ds.GetGeoTransform()[1]
        self.yres = self.ds.GetGeoTransform()[5]
        self.right = self.left + self.xsize * self.xres
        self.bottom = self.top + self.ysize * self.yres
        self.projection = self.ds.GetProjection()
        pix_lin, pix_col = np.indices((self.ds.RasterYSize, self.ds.RasterXSize))
        self.lat, self.lon = self.top + self.yres*pix_lin, self.left+self.xres*pix_col

        # convert 0 and 255 to NaN
        self.data[self.data==0.] = np.nan
        self.data[self.data==255] = np.nan
        self.data[self.data == -9999] = np.nan

        if sigfile is not None:
            self.dst = gdal.Open(sigfile)
            self.bandt = self.dst.GetRasterBand(1)
            self.sigma = self.bandt.ReadAsArray()
            self.sigma[self.sigma == 0] = np.nan
            if self.dst.RasterXSize != self.xsize or self.dst.RasterYSize != self.ysize:
                try:
                    self.sigma = self.sigma[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Sigma and Velocity file not the same size!')
                    print('sig has size = ' + str(self.dst.RasterXSize) + ', ' + str(self.dst.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))
                # self.sigma = np.ones((self.ysize, self.xsize))
        # else:
        #     self.sigma = np.ones((self.ysize, self.xsize))

        if incidence is not None:
            self.ds_inc = gdal.Open(incidence)
            self.band_inc = self.ds_inc.GetRasterBand(1)
            self.inc = np.deg2rad(self.band_inc.ReadAsArray())
            self.inc[self.inc == 0] = np.nan
            if self.ds_inc.RasterXSize != self.xsize or self.ds_inc.RasterYSize != self.ysize:
                try:
                    self.inc = self.inc[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Inc and Velocity file not the same size!')
                    print('inc has size = ' + str(self.ds_inc.RasterXSize) + ', ' + str(self.ds_inc.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if heading is not None:
            self.ds_head = gdal.Open(heading)
            self.band_head = self.ds_head.GetRasterBand(1)
            self.head = np.deg2rad(self.band_head.ReadAsArray())
            self.head[self.head == 0] = np.nan
            if self.ds_head.RasterXSize != self.xsize or self.ds_head.RasterYSize != self.ysize:
                try:
                    self.head = self.head[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Heading and Velocity file not the same size!')
                    print('head has size = ' + str(self.ds_head.RasterXSize) + ', ' + str(self.ds_head.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if N is not None:
            self.ds_N = gdal.Open(N)
            self.band_N = self.ds_N.GetRasterBand(1)
            self.N = self.band_N.ReadAsArray()
            self.N[self.N == 0] = np.nan
            if self.ds_N.RasterXSize > self.xsize or self.ds_N.RasterYSize > self.ysize:
                self.N = self.N[:self.ysize, :self.xsize]
            if self.ds_N.RasterXSize < self.xsize or self.ds_N.RasterYSize < self.ysize:
                self.N = self.N[:self.ysize, :self.xsize]
                tmp = np.ones((self.ysize, self.xsize))*np.nan
                tmp[:self.ds_N.RasterYSize, :self.ds_N.RasterXSize] = self.N
                self.N = tmp

        if E is not None:
            self.ds_E = gdal.Open(E)
            self.band_E = self.ds_E.GetRasterBand(1)
            self.E = self.band_E.ReadAsArray()
            self.E[self.E == 0] = np.nan
            if self.ds_E.RasterXSize != self.xsize or self.ds_E.RasterYSize != self.ysize:
                try:
                    self.E = self.E[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Heading and Velocity file not the same size!')
                    print('head has size = ' + str(self.ds_E.RasterXSize) + ', ' + str(self.ds_E.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if U is not None:
            self.ds_U = gdal.Open(U)
            self.band_U = self.ds_U.GetRasterBand(1)
            self.U = self.band_U.ReadAsArray()
            self.U[self.U == 0] = np.nan
            if self.ds_U.RasterXSize != self.xsize or self.ds_U.RasterYSize != self.ysize:
                try:
                    self.U = self.U[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Heading and Velocity file not the same size!')
                    print('head has size = ' + str(self.ds_U.RasterXSize) + ', ' + str(self.ds_U.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if np.logical_and(incidence is not None, heading is not None):
            # VLos = VE (-cos(head)sin(inc) + VN (sin(head)sin(inc)) + VU cos(inc)
            #      = VE (-cos(head)sin(inc) + VUN (sqrt(1 - sin^2(inc)cos^2(head)))
            self.E = - np.cos(self.head) * np.sin(self.inc)
            self.N = np.sin(self.head) * np.sin(self.inc)
            self.U = np.cos(self.inc)

    # def shape_array_to_self(self, arr):
    #     if arr.shape[1] > self.xsize or arr.shape[0] > self.ysize:
    #         arr = arr[:self.ysize, :self.xsize]
    #     if arr.shape[1] < self.xsize or arr.shape[0] < self.ysize:
    #         tmp = np.ones((self.ysize, self.xsize))*np.nan
    #         tmp[:arr.shape[0], :arr.shape[1]] = arr
    #         return tmp

    def clean_by_sigma(self, threshold):
        clean_data = copy.copy(self.data)
        clean_data[self.sigma > threshold] = np.nan
        clean_sigma = copy.copy(self.sigma)
        clean_sigma[np.isnan(clean_data)] = np.nan

        vmin = np.nanpercentile(self.data, 1)
        vmax = np.nanpercentile(self.data, 99)
        smin = np.nanpercentile(self.sigma, 1)
        smax = np.nanpercentile(self.sigma, 99)

        fig, ax = plt.subplots(2, 2, sharey='all', sharex='all')
        im = ax[0, 0].imshow(self.data, vmin=vmin, vmax=vmax)
        ax[0, 0].set_title("Raw Data")
        ax[0, 1].imshow(clean_data, vmin=vmin, vmax=vmax)
        ax[0, 1].set_title("Clean Data")
        fig.colorbar(im, ax=ax[0, :], shrink=0.8, label="LOS / mm/yr")

        im = ax[1, 0].imshow(self.sigma, vmin=smin, vmax=smax)
        ax[1, 0].set_title("Raw Sigma")
        ax[1, 1].imshow(clean_sigma, vmin=smin, vmax=smax)
        ax[1, 1].set_title("Clean Sigma")
        fig.colorbar(im, ax=ax[1, :], shrink=0.8, label="SIGMA / mm/yr")

        plt.suptitle(self.basename)
        plt.show()

        return clean_data, clean_sigma

    def extract_pixel_value(self, lon, lat, max_width=200):
        x = int((lon - self.left)/self.xres + 0.5)
        y = int((lat - self.top) / self.yres + 0.5)
        # increase window size in steps of 2 until there are non-nan values in the window
        # starting from 2 with 5x5 window because if 1x1 window, stdev will be zero
        # if we use the std of values instead of the corresponding sigma files as stdev
        for n in np.arange(2, max_width+1, 2):
            pixel_values = self.data[y - n: y + n + 1, x - n: x + n + 1]
            index = np.nonzero(~np.isnan(pixel_values))
            if len(index[0]) > 10:
                # print(n, pixel_values)
                break
        pixel_value = np.nanmean(pixel_values)
        # weighted mean with sigma
        # pixel_sigma = self.sigma[y - n: y + n + 1, x - n: x + n + 1]
        # print(pixel_sigma)
        # pixel_value = np.nansum(pixel_values[index] * pixel_sigma[index]) / np.nansum(pixel_sigma[index])   # takes care of the case where sigma is nan but data is not nan
        stdev = np.nanstd(pixel_values)  # by using nanstd(pixel_values), we are not taking into account the quality of the pixels here.
        return pixel_value, stdev

    def extract_inc(self, lon, lat):
        x = int((lon-self.left)/self.xres+0.5)
        y = int((lat - self.top) / self.yres + 0.5)
        inc = self.inc[y, x]
        return inc

    def extract_head(self, lon, lat):
        x = int((lon-self.left)/self.xres+0.5)
        y = int((lat - self.top) / self.yres + 0.5)
        head = self.head[y, x]
        return head


class Overlap:
    """ a class that calculates and stores the overlapping boundaries"""
    def __init__(self, d1, d2, vmin=None, vmax=None):
        self.d1 = d1
        self.d2 = d2

        ''' a method to calculate the overlapping area between two images'''
        self.left = np.maximum(self.d1.left, self.d2.left)
        self.right = np.minimum(self.d1.right, self.d2.right)
        self.top = np.minimum(self.d1.top, self.d2.top)
        self.bottom = np.maximum(self.d1.bottom, self.d2.bottom)
        if self.left > self.right or self.top < self.bottom :
            raise ValueError('Two images do not overlap.')

        ''' a method to crop arrays from two images into overlapping area'''
        self.d1xl = int((self.left-self.d1.left)/self.d1.xres+0.5)
        self.d1xr = int((self.right-self.d1.left)/self.d1.xres+0.5)
        self.d1yt = int((self.top-self.d1.top)/self.d1.yres+0.5)
        self.d1yb = int((self.bottom-self.d1.top)/self.d1.yres+0.5)
        self.d2xl = int((self.left-self.d2.left)/self.d2.xres+0.5)
        self.d2xr = int((self.right-self.d2.left)/self.d2.xres+0.5)
        self.d2yt = int((self.top-self.d2.top)/self.d2.yres+0.5)
        self.d2yb = int((self.bottom-self.d2.top)/self.d2.yres+0.5)
        miss = (self.d1xr - self.d1xl) - (self.d2xr - self.d2xl)
        # if miss != 0:
        #     print("x-dimensions are different...")
        #     print(self.d1xr, self.d1xl, self.d2xr, self.d2xl)
        if miss < 0:
            self.d2xr += miss
        elif miss > 0:
            self.d1xr -= miss
        miss = (self.d1yb - self.d1yt) - (self.d2yb - self.d2yt)
        # if miss != 0:
        #     print("y-dimensions are different...")
        #     print(self.d1yb, self.d1yt, self.d2yb, self.d2yt)
        if miss < 0:
            self.d2yb += miss
        elif miss > 0:
            self.d1yb -= miss
        self.d1array = self.d1.data[self.d1yt: self.d1yb, self.d1xl: self.d1xr]
        self.d2array = self.d2.data[self.d2yt: self.d2yb, self.d2xl: self.d2xr]

        if hasattr(d1, 'sigma') and hasattr(d1, 'sigma'):
            self.d1sigma = self.d1.sigma[self.d1yt: self.d1yb, self.d1xl: self.d1xr]
            self.d2sigma = self.d2.sigma[self.d2yt: self.d2yb, self.d2xl: self.d2xr]
            try:
                self.sigma = np.sqrt(self.d1sigma.flatten() ** 2 + self.d2sigma.flatten() ** 2)
            except:
                print(self.d1xr, self.d1xl, self.d2xr, self.d2xl)
                print(self.d1yb, self.d1yt, self.d2yb, self.d2yt)
                print(self.d1sigma.shape)
                print(self.d2sigma.shape)
                breakpoint()
            # try:
            #     self.sigma = np.sqrt(self.d1sigma.flatten()**2 + self.d2sigma.flatten()**2)
            # except:
            #     self.d1xr = self.d1xl + self.d1sigma.shape[1]
            #     self.d1yb = self.d1yt + self.d1sigma.shape[0]
            #     self.d2xr = self.d2xl + self.d1sigma.shape[1]
            #     self.d2yb = self.d2yt + self.d1sigma.shape[0]
            #     # breakpoint()
            #     self.d1sigma = self.d1.sigma[self.d1yt: self.d1yb, self.d1xl: self.d1xr]
            #     self.d2sigma = self.d2.sigma[self.d2yt: self.d2yb, self.d2xl: self.d2xr]
            #     self.sigma = np.sqrt(self.d1sigma.flatten() ** 2 + self.d2sigma.flatten() ** 2)
            #     self.d1array = self.d1.data[self.d1yt: self.d1yb, self.d1xl: self.d1xr]
            #     self.d2array = self.d2.data[self.d2yt: self.d2yb, self.d2xl: self.d2xr]

        # print(self.d1.data.shape, self.d1array.shape, self.d2.data.shape, self.d2array.shape)

        self.diff_array = np.subtract(self.d1array, self.d2array)
        if vmin:
            self.diff_array[self.d1array < vmin] = np.nan
            self.diff_array[self.d2array < vmin] = np.nan
        if vmax:
            self.diff_array[self.d1array > vmax] = np.nan
            self.diff_array[self.d2array > vmax] = np.nan

        self.sum_array = np.add(self.d1array, self.d2array)
        self.vector_sum_array = np.sqrt(np.add(np.square(self.d1array), np.square(self.d2array)))
        self.mode = mode(np.around(self.diff_array[~np.isnan(self.diff_array)], decimals=2))[0][0]
        # print('mode=%.2f' % self.mode)
        pix_lin, pix_col = np.indices(self.d1array.shape)
        self.lat, self.lon = self.top + self.d1.yres*pix_lin, self.left+self.d1.xres*pix_col


class Merge:
    """ a class that combines two imagery through a constant offset """
    def __init__(self, d1, d2):
        self.d1 = d1
        self.d2 = d2

        # define geographic coordinates of merged matrix and an empty matrix of the merged shape
        self.left = np.minimum(self.d1.left, self.d2.left)
        self.right = np.maximum(self.d1.right, self.d2.right)
        self.top = np.maximum(self.d1.top, self.d2.top)
        self.bottom = np.minimum(self.d1.bottom, self.d2.bottom)
        self.xsize = int((self.right - self.left) / min(self.d1.xres, self.d2.xres) + 1.5)  # 1.5 to ensure the merged matrix has enough columns to host two matrix
        self.ysize = int((self.top - self.bottom) / -min(self.d1.yres, self.d2.yres) + 1.5)  # 1.5 to ensure the merged matrix has enough rows to host two matrix
        self.array = np.ones((self.ysize, self.xsize), np.int32) * np.nan

        # place two matrices into the right positions within the merged matrix
        self.d1xl = int((self.d1.left-self.left)/self.d1.xres+0.5)
        # self.d1xr = int((self.d1.right-self.left)/self.d1.xres+0.5)
        self.d1xr = self.d1xl + self.d1.data.shape[1]
        self.d1yt = int((self.d1.top-self.top)/self.d1.yres+0.5)
        # self.d1yb = int((self.d1.bottom-self.top)/self.d1.yres+0.5)
        self.d1yb = self.d1yt + self.d1.data.shape[0]
        self.d2xl = int((self.d2.left-self.left)/self.d2.xres+0.5)
        # self.d2xr = int((self.d2.right-self.left)/self.d2.xres+0.5)
        self.d2xr = self.d2xl + self.d2.data.shape[1]
        self.d2yt = int((self.d2.top-self.top)/self.d2.yres+0.5)
        # self.d2yb = int((self.d2.bottom-self.top)/self.d2.yres+0.5)
        self.d2yb = self.d2yt + self.d2.data.shape[0]

    def non_nan_merge(self):
        # Place the top frame in to the merged matrix first
        self.array[self.d1yt:self.d1yb, self.d1xl:self.d1xr] = self.d1.data

        # Choose the non-nan positions from d2 and fill in with d2 data,
        # keep existing values in merged matrix (d1 data) in positions where d2 is nan.
        nodata_test = np.isnan(self.d2.data)  # = True if nan, = False if not nan; True = 1, False = 0
        masked_data2 = np.choose(nodata_test, (self.d2.data, self.array[self.d2yt:self.d2yb, self.d2xl:self.d2xr]))
        # masked_data2 = np.choose(nodata_test, (self.d2.data, self.array[self.d2yt:self.d2yt++self.d2.data.shape[0], self.d2xl:self.d2xl+self.d2.data.shape[1]]))
        # False = not nan = pick from 0th entry; True = nan = pick from 1st entry;
        # 1st entry already contains value from the top matrix
        self.array[self.d2yt:self.d2yb, self.d2xl:self.d2xr] = masked_data2

    def count_band(self):
        # Place the top frame in to the merged matrix first
        nodata_test = np.isnan(self.d1.data)
        band_count1 = np.choose(nodata_test, (1, 0))  # band_count = 1 if not nan, = 0 if nan.
        nodata_test = np.isnan(self.d2.data)
        band_count2 = np.choose(nodata_test, (1, 0))  # band_count = 1 if not nan, = 0 if nan.

        self.array[self.d1yt:self.d1yb, self.d1xl:self.d1xr] = 0  # first set to 0 to avoid nan+1=nan problem
        self.array[self.d2yt:self.d2yb, self.d2xl:self.d2xr] = 0  # first set to 0 to avoid nan+1=nan problem
        self.array[self.d1yt:self.d1yb, self.d1xl:self.d1xr] += band_count1
        self.array[self.d2yt:self.d2yb, self.d2xl:self.d2xr] += band_count2

    def export_tif(self, export_title):
        # Export merged data to tif format.
        driver = gdal.GetDriverByName("GTiff")
        outdata = driver.Create(export_title+'.tif', self.xsize, self.ysize, 1, gdal.GDT_Float32)
        outdata.SetGeoTransform([self.left, self.d1.xres, 0, self.top, 0, self.d1.yres])  ##sets same geotransform as input
        outdata.SetProjection(self.d1.projection)  ##sets same projection as input
        outdata.GetRasterBand(1).WriteArray(self.array)


def plot_hist(array, hist_fname, plot_title):
    """ plot histogram and the array in two subplots, top and bottom if the array is fat, left and right if thin """
    if array.shape[0] < array.shape[1]:
        fig, axes = plt.subplots(
            nrows=2, ncols=1, sharex=False, sharey=False,
            gridspec_kw={'width_ratios': [1]}, figsize=(3.2, 2.4)
        )

        axes[0].set_title(plot_title, fontsize=12)
    else:
        fig, axes = plt.subplots(
            nrows=1, ncols=2, sharex=False, sharey=False,
            gridspec_kw={'width_ratios': [4, 3]},
            figsize=(3.8, 3.)
        )
        fig.suptitle(plot_title, fontsize=14)  #[-16:]


    vmin = np.nanpercentile(array, 0.5)
    vmax = np.nanpercentile(array, 99.5)

    axes[0].set_xlim(-10, 10)
    axes[0].tick_params(axis='both', which='major', labelsize=12)
    axes[0].yaxis.set_visible(False)
    axes[1].yaxis.set_visible(False)
    axes[1].xaxis.set_visible(False)
    data = array.flatten()
    sns.histplot(data[(data > vmin) & (data < vmax)], kde=True,
                 bins=100, color='darkblue',
                 ax=axes[0])

    # annotate with stats
    nanmode = mode(np.around(array[~np.isnan(array)], decimals=2))[0][0]

    axes[0].annotate("m=%.1f \n s=%.1f " % (
        nanmode, np.nanstd(array)),
                     xy=(0.97, 0.93), xycoords='axes fraction',
                     ha="right", va="top",
                     fontsize=14)
    print('mode=%.2f' % nanmode)
    axes[0].axvline(x=nanmode, c='white', lw=2, linestyle='-')
    axes[0].axvline(x=nanmode - np.nanstd(array), c='white', lw=1, linestyle='--')
    axes[0].axvline(x=nanmode + np.nanstd(array), c='white', lw=1, linestyle='--')

    # plot overlapping area in the lower panel
    im2 = axes[1].imshow(array, vmax=vmax, vmin=vmin, cmap=cm.roma.reversed())

    axes[1].set_xlim(0, array.shape[1] * 1.05)
    # colorbar for the overlapping area
    if array.shape[0] < array.shape[1]:
        axins1 = inset_axes(axes[1], width="3%", height="100%", loc='right', )
        c = fig.colorbar(im2, cax=axins1)
        c.ax.tick_params(labelsize=12)
    else:
        c = fig.colorbar(im2)
        c.ax.tick_params(labelsize=12)
        c.set_label('mm/yr')
    plt.show()
    fig.savefig(hist_fname, format='PNG', dpi=300, bbox_inches='tight', transparent=True)
    return nanmode


def plot_merge(array, plot_title):
    vmin = np.nanpercentile(array, 0.5)
    vmax = np.nanpercentile(array, 99.5)
    vmin = max(vmin, -10)
    fig1, ax = plt.subplots(1,1)
    im = ax.imshow(array, vmax=vmax, vmin=vmin, cmap=cm.roma.reversed())
    ax.set_title(os.path.basename(plot_title))
    c = plt.colorbar(im)
    c.set_label('LOS mm/yr')
    plt.show()
    fig1.savefig(plot_title+'.png', format='PNG', dpi=300, bbox_inches='tight')


def deramp_along_range(overlap_array, theta_deg, title):
    """Invert ramp along theta direction and calculate residual"""
    fig1, axs = plt.subplots(3, 1, sharex='col', constrained_layout=True)

    # Find the left most index with non-nan value in the overlapping region ###
    vmin = np.nanpercentile(overlap_array, 0.5)
    vmax = np.nanpercentile(overlap_array, 99.5)
    im = axs[0].imshow(overlap_array, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
    axs[0].set_title(title[-21:] + "_offset")

    non_nan_indices = np.argwhere(~np.isnan(overlap_array))  # non-nan rows (ys), columns (xs)
    data_phi = overlap_array[~np.isnan(overlap_array)]
    (left_most_y, left_most_x) = non_nan_indices[np.argmin(non_nan_indices[:, 1])] # a long matrix of (row, column) = ([y, x])
    swap_row_matrix = [[0, 1], [1, 0]]
    data_xy = np.dot(swap_row_matrix, non_nan_indices.T)
    data_xy0 = data_xy - [[left_most_x], [left_most_y]]  # shift origin to be the left most non-nan pixel

    theta = theta_deg * np.pi / 180

    R = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]  # rotational matrix
    data_xy_primes = np.dot(R, data_xy0)  # rotate axis
    data_size = data_xy_primes.shape[1]

    # Invert for solution
    G = np.ones((data_size, 2))
    x_prime = data_xy_primes[0]
    G[:, 0] = x_prime  # only with non-nan data
    (m, c), res, rank, singular = np.linalg.lstsq(G, data_phi, rcond=None)
    print('theta=%.2f \nm=%.5f \nc=%.2f' % (theta_deg, m, c))

    # Reconstruct model solution where overlap contains value
    model = np.dot(G, (m, c))
    (height, width) = overlap_array.shape
    model_array = np.empty((height * width))
    model_array[:] = np.nan
    model_array[~np.isnan(overlap_array.flatten())] = model
    model_array_2d = model_array.reshape(height, width)
    axs[1].imshow(model_array_2d, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())

    # Calculate residual
    residual = overlap_array - model_array_2d
    axs[2].imshow(residual, vmin=vmin, vmax=vmax, cmap=cm.roma.reversed())
    axs[2].set_title('residual')

    fig1.colorbar(im, ax=axs)
    plt.show()
    # fig1.savefig('../png/' + title[-21:] + '_deramp.png', format='PNG', dpi=300, bbox_inches='tight')

    return residual, m, c, left_most_y, left_most_x, R


def calc_ramp(array, m, c, left_most_y, left_most_x, rotational_matrix):
    """Reconstruct model solution to cover the whole array with a ramp"""
    (height, width) = array.shape
    x = np.arange(0, width)
    y = np.arange(0, height)
    x0 = x - left_most_x
    y0 = y - left_most_y
    xs = np.tile(x0, height)
    ys = np.repeat(y0, width)
    xy = np.asarray([xs, ys])
    xy_rotated = np.dot(rotational_matrix, xy)
    xs_rotated = xy_rotated[0]
    G_full = np.ones((height * width, 2))
    G_full[:, 0] = xs_rotated
    model_array_2d = np.dot(G_full, (m, c)).reshape(height, width)
    return model_array_2d


def plot_data_plus_ramp(data, model, title):
    """To visualise the effect of adding a ramp to the original array"""
    fig, axs = plt.subplots(1, 3, sharey='row', constrained_layout=True, figsize=(6.4, 3))
    fig.suptitle(title)
    vmin = np.nanpercentile(data, 0.5)
    vmax = np.nanpercentile(data, 99.5)
    im = axs[0].imshow(data, alpha=0.8, vmin=vmin, vmax=vmax)
    axs[0].set_title("data")
    axs[1].imshow(model, alpha=0.8, vmin=vmin, vmax=vmax)
    axs[1].set_title("ramp")
    axs[2].imshow(data+model, alpha=0.8, vmin=vmin, vmax=vmax)
    axs[2].set_title("plus_ramp")
    plt.colorbar(im, ax=axs, orientation='horizontal')
    plt.show()
    # fig.savefig('../png/' + title + '_plus_ramp.png', format='PNG', dpi=300, bbox_inches='tight')


def get_heading(frame):
    """get the heading angle either based on letter A or D in the frame name or from the corresponding parfile"""
    # set default
    if "A" in frame:
        heading = float(13.5)
    else:
        heading = float(166.5)

    # get more accurate heading from parameter file if available
    par = "../para/" + frame + ".par"
    if os.path.exists(par):
        for line in open(par):
            if "heading" in line:
                heading = -float(line.split()[1])  # as the y axis increases downwards, need to flip the sign of heading
    print("heading =", heading)
    return heading


def combine2tifs(t1, t2, stitch_style, rf, outdir, rows=0):   #rf is defined with an "_" in front
    """Master function that calls the other functions"""
    ds1 = OpenTif(t1)
    ds2 = OpenTif(t2)
    overlap = Overlap(ds1, ds2)
    base_title = os.path.basename(t2)[:17]
    # base_title = make_title(os.path.basename(t1), os.path.basename(t2), stitch_style, rf)  # generate output name
    hist_file = os.path.join(outdir, base_title + rf + stitch_style + "_offset.png")  # to be saved as filename
    if rows == 0:
        diff_array = overlap.diff_array
        offset = overlap.mode
        hist_title = base_title  # + "_offset"  # to be displayed in plot
    else:
        nonnan_rows = overlap.diff_array[~np.isnan(overlap.diff_array).all(axis=1)]
        diff_array = nonnan_rows[:rows]
        offset = mode(np.around(diff_array[~np.isnan(diff_array)], decimals=2))[0][0]
        hist_title = base_title+str(rows)
    plot_hist(diff_array, hist_file, hist_title)
    if np.isnan(offset):
        offset = 0

    if stitch_style[-4:] == 'mode':  # stitch by removing a mode from offset array
        ds2.data += offset

    if (stitch_style[-5:] == 'range') or (stitch_style[-7:] == 'azimuth'):  # stitch by deramp along range or azimuth
        frame = os.path.basename(t2)[:17]
        if stitch_style[-5:] == 'range':
            heading = get_heading(frame)
        else:
            heading = get_heading(frame)+90
        residual, m, c, left_most_y, left_most_x, rotation = deramp_along_range(overlap.diff_array, heading, base_title)
        res_file = os.path.join(outdir, base_title + rf + stitch_style + "_residual.png")  # to be saved as filename
        res_title = base_title[-21:] #+ "_residual"  # to be displayed in plot
        plot_hist(residual, res_file, res_title)
        # originShiftX = overlap.d2xl + left_most_x
        # originShiftY = overlap.d2yt + left_most_y
        # model = calc_ramp2d(ds2.data, a, b, c, originShiftY, originShiftX, rotation)
        # model = calc_ramp(ds2.data, m, c, originShiftY, originShiftX, rotation)
        model = calc_ramp(ds2.data, m, c, left_most_y, left_most_x, rotation)
        plot_data_plus_ramp(ds2.data, model, base_title)  # to visualise the change before and after the ramp
        ds2.data += model

    # if (stitch_style[-6:] == 'planar'):  # stitch by fitting a planar ramp, ax+by+c, to the offset
    #     # over = Overlap(ds[i], ds[i + 1])
    #     # sigover = Overlap(sigma[i], sigma[i + 1])
    #     # x_shift = int((overlap.left - track_left) / ds[i].xres + 0.5)
    #     # y_shift = int((overlap.top - track_top) / ds[i].yres + 0.5)
    #     non_nan_yx = np.argwhere(~np.isnan(overlap.diff_array))
    #     non_nan_diff = overlap.diff_array[~np.isnan(overlap.diff_array)]
    #     # if over.diff_array.shape != sigover.sum_array.shape:
    #     #     (row, col) = over.diff_array.shape
    #     #     sigover.sum_array = sigover.sum_array[:row, :col]  # one case in 135D encountered dimension mismatch
    #     # non_nan_sig = sigover.sum_array[~np.isnan(over.diff_array)]  # extract the same pixels from sigma overlap
    #     # non_nan_sig[np.isnan(non_nan_sig)] = np.nanmean(non_nan_sig)  # in case of void, fill with nanmean
    #     # non_nan_yx += [y_shift, x_shift]
    #     G = np.zeros((non_nan_yx.shape[0], 6))
    #     G[:, 0: 2] = non_nan_yx
    #     # G[:, 3: 5] = -non_nan_yx
    #     G[:, 2] = 1
    #     # G[:, 5] = -1
    #
    #     # inversion
    #     coefs, res, rank, singular = np.linalg.lstsq(G, non_nan_diff, rcond=None)
    #
    #     # reconstruct the full ramp image and export the projected insar los velocity map
    #     # reconstruct ramp
    #     full_ramp_array = np.zeros(ds[i].data.shape)
    #     all_yx = np.argwhere(~np.isnan(full_ramp_array))
    #     x_shift = int((ds[i].left - track_left) / ds[i].xres + 0.5)
    #     y_shift = int((ds[i].top - track_top) / ds[i].yres + 0.5)
    #     all_yx += [y_shift, x_shift]
    #     G = np.ones((len(all_yx), 3))
    #     G[:, :2] = all_yx
    #     full_ramp_array = np.dot(G, coefs[3*i:3*i+3]).reshape(ds[i].ysize, ds[i].xsize)
    #
    #     # calc projected insar
    #     projected_insar_array = ds[i].data - full_ramp_array
    # breakpoint()

    merge = Merge(ds1, ds2)
    merge.non_nan_merge()  # first plot tif1, then plot tif2, filling nan pixels in overlapping areas with tif2 values.

    output_basename = base_title + rf + stitch_style  # generate output filename
    merge.export_tif(os.path.join(outdir, output_basename))  # export merged tif
    plot_merge(merge.array, os.path.join(outdir, output_basename))  # plot merged array
    return output_basename + '.tif'


def make_title(t1, t2, merge_style, ref):
    """Combine two input file names to generate the output filename"""
    t1 = re.sub(ref, '', t1)
    t2 = re.sub(ref, '', t2)
    t1 = re.sub(merge_style, '', t1)
    t2 = re.sub(merge_style, '', t2)
    t1 = re.sub('_131313', '', t1)
    t2 = re.sub('_131313', '', t2)
    regex1 = re.compile('[a-z.]')  # remove any lower case letters
    t1 = regex1.sub('', t1)
    t2 = regex1.sub('', t2)
    regex2 = re.compile('_*$')  # remove any trailing underscores
    t1 = regex2.sub('', t1)
    t2 = regex2.sub('', t2)
    basetitle = (t1 + '-' + t2)
    print('combine ' + t1 + ' with ' + t2 + ' for ' + basetitle)
    return basetitle


def export_tif(data, df, filename):
    # Export data to tif format.
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(filename, df.xsize, df.ysize, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([df.left, df.xres, 0, df.top, 0, df.yres])  ##sets same geotransform as input
    # outdata.SetProjection(df.projection)  ##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(data)
    outdata.FlushCache()
    outdata.FlushCache()  # need to flush twice to export the last tif properly, otherwise it stops halfway.


if __name__ == "__main__":
    #############################
    # Input parameters:
    style = '_mode'  # can be empty string, _*range, _*mode, or _*azimuth _range

    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../../dehua/insardata20230413_1.5_m_Dehua/frames/'  # where input velocity frames are stored
    # input_suffix = '_U.tif'
    # merge_dir = input_dir + 'merge/'
    # output_dir = input_dir + 'tracks/'

    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../los_frames/'  # where input velocity frames are stored
    # input_suffix = '.fringe.vel.geo.tif'
    # merge_dir = '../los_frames/merge/'
    # output_dir = '../los_frames/track/'

    # reference = ''  # like a prefix in the input filenames
    # input_dir = '/Users/qi/Downloads/gacos_vel/tif/'  # where input velocity frames are stored
    # input_suffix = '_gacos_vel.tif'
    # merge_dir = '/Users/qi/Downloads/gacos_vel/tif/merge/'
    # output_dir = '/Users/qi/Downloads/gacos_vel/tif/output/'

    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../vel_gacos/'  # where input velocity frames are stored
    # input_suffix = '_gacos_vel.tif'
    # merge_dir = '../vel_gacos/merge/'
    # output_dir = '../vel_gacos/tracks/'

    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../vel_nogacos/'  # where input velocity frames are stored
    # input_suffix = '_vel_nogacos.tif'
    # merge_dir = '../vel_nogacos/merge/'
    # output_dir = '../vel_nogacos/tracks/'



    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../vel_mskd/'  # where input velocity frames are stored
    # input_suffix = '.fringe.vel.mskd.tif'
    # merge_dir = '../vel_mskd_merge/'
    # output_dir = '../vel_mskd_tracks/'
    #
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../without_plate_motion/'  # where input velocity frames are stored
    # input_suffix = '.vel.no_plate.tif'
    # merge_dir = '../los_frame_merge_no_plate/'
    # output_dir = '../los_tracks_no_plate/'
    #
    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../corrected_vstd/'  # where input velocity frames are stored
    # input_suffix = '_scaled_vstd.tif'
    # merge_dir = '../corrected_vstd_tracks_merge/'
    # output_dir = '../corrected_vstd_tracks/'

    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../corrected_vstd/'  # where input velocity frames are stored
    # input_suffix = 'vstd_scaled_vstd.tif'
    # merge_dir = '../corrected_vstd/merge/'
    # output_dir = '../corrected_vstd/frame/'

    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../referenced_by_gps_overlap/'  # where input velocity frames are stored
    # input_suffix = '_tarim.tif'
    # merge_dir = '../referenced_by_gps_overlap/merge/'
    # output_dir = '../referenced_by_gps_overlap/track/'

    style = '_N'  # can be empty string, _*range, _*mode, or _*azimuth _range
    reference = ''  # like a prefix in the input filenames
    input_dir = '../NEU/'  # where input velocity frames are stored
    input_suffix = '.N.ml10.tif'
    merge_dir = '../NEU/merge/'
    output_dir = '../NEU/track/'

    # style = '_hgt'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../hgt/'  # where input velocity frames are stored
    # input_suffix = '.hgt.ml10.tif'
    # merge_dir = '../hgt_tracks_merge/'
    # output_dir = '../hgt_tracks/'

    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../without_plate_motion/'  # where input velocity frames are stored
    # input_suffix = '.vel.plate_motion.tif'
    # merge_dir = '../plate_motion_merge/'
    # output_dir = '../plate_motion_tracks/'


    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../los_frames/ml5/merge_pair/final_3/'  # where input velocity frames are stored
    # input_suffix = '.fringe.vel.geo.tif'
    # merge_dir = '../los_frame_merge/'
    # output_dir = '../los_tracks/'


    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../los_weighted/referenced_by_gps_overlap'  # where input velocity frames are stored
    # input_suffix = '*_improved_few_tarim.tif'
    # merge_dir = '../los_weighted/referenced_by_gps_overlap/merge'
    # output_dir = '../los_weighted/referenced_by_gps_overlap/track_unmasked'

    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../los_full/referenced_by_gps_overlap'  # where input velocity frames are stored
    # input_suffix = '*_3.tif'
    # merge_dir = '../los_full/referenced_by_gps_overlap/merge'
    # output_dir = '../los_full/referenced_by_gps_overlap/track'
    # #

    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../los_weighted/referenced_by_gps_overlap'  # where input velocity frames are stored
    # input_suffix = '*masked.tif'
    # merge_dir = '../los_weighted/referenced_by_gps_overlap/merge'
    # output_dir = '../los_weighted/referenced_by_gps_overlap/track_masked'

    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../los_weighted/vstd_corrected/'  # where input velocity frames are stored
    # input_suffix = '*.tif'
    # merge_dir = '../los_weighted/vstd_corrected/merge'
    # output_dir = '../los_weighted/vstd_corrected/track'

    # style = '_none'  # can be empty string, _*range, _*mode, or _*azimuth _range
    # reference = ''  # like a prefix in the input filenames
    # input_dir = '../los_full/vstd_corrected/'  # where input velocity frames are stored
    # input_suffix = '*.tif'
    # merge_dir = '../los_full/vstd_corrected/merge'
    # output_dir = '../los_full/vstd_corrected/track'

    ###############################

    # 1. list all frames with the correct suffix in the directory
    tifList = sorted(glob.glob(os.path.join(input_dir, '107*'+input_suffix)))
    print(tifList)

    Path(merge_dir).mkdir(parents=True, exist_ok=True)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # 2. define tracks from frame names
    trackList = set([os.path.basename(t)[:4] for t in tifList])

    # 3. for each track, loop through all frames along track,
    # do combine2tifs on consecutively merged frames with the next frame
    for track in trackList:
        frameList = sorted(glob.glob(os.path.join(input_dir, track+'*'+input_suffix)))
        # frameList = ['../los_frames/092D_04050_141516.fringe.vel.geo.tif','../los_frames/092D_04271_141413.fringe.vel.geo.tif','../los_frames/092D_04475_131313.fringe.vel.geo.tif']
        # frameList = ['../los_frames/092D_04679_141413.fringe.vel.geo.tif','../los_frames/092D_04877_111214.fringe.vel.geo.tif']
        # frameList = ['../los_tracks/092D__mode_1.tif','../los_tracks/092D__mode_2.tif']
        # frameList = glob.glob(os.path.join(input_dir, track+'*04294*'+input_suffix))

        print(frameList)
        count = 0
        tif1 = frameList[count]
        for count in range(len(frameList)-1):
            print(count)
            tif2 = frameList[count+1]
            merged_basename = combine2tifs(tif1, tif2, style, reference, merge_dir)
            tif1 = os.path.join(merge_dir, merged_basename)
            count += 1
            if count == len(frameList)-1:
                shutil.copy(tif1, os.path.join(output_dir, track+style+".tif"))
                shutil.copy(tif1[:-4]+".png", os.path.join(output_dir, track+style+".png"))
                for f in glob.glob(os.path.join(merge_dir, "*.tif")):
                    os.remove(f)
                for f in glob.glob(os.path.join(merge_dir, "*"+style+".png")):
                    os.remove(f)

#
# dir = "../los_frames/ml5/merge_pair/"
# tif1 = dir + '056A_04661_131313.fringe.vel.geo.tif'
# tif2 = dir + '056A_04947_282019.fringe.vel.geo.tif'
# tif1 = dir + '056A_04414_191919_merged_mskd.tif'
# tif2 = dir + "056A_04947_282019_mode.tif"
# combine2tifs(tif1, tif2, '_mode', "", ".")


#
# dir = "../without_plate_motion/"
# tif1 = dir + '092D_04271_141413.vel.no_plate.tif'
# tif2 = dir + '092D_04475_131313.vel.no_plate.tif'
# outdir = "../manual_stitching/"
# combine2tifs(tif1, tif2, '_mode', "", outdir)
#
# dir = "../without_plate_motion/"
# tif1 = dir + '092D_04679_141413.vel.no_plate.tif'
# tif2 = dir + '092D_04877_111214.vel.no_plate.tif'
# outdir = "../manual_stitching/"
# combine2tifs(tif1, tif2, '_mode', "", outdir)
#
# dir = "../manual_stitching/"
# tif1 = dir + '092D_04475_131313_mode.tif'
# tif2 = dir + '092D_04877_111214_mode.tif'
# outdir = "../manual_stitching/"
# combine2tifs(tif1, tif2, '_none', "", outdir)

#
# dir = "../without_plate_motion/"
# tif2 = dir + '092D_04050_141516.vel.no_plate.tif'
# dir = "../manual_stitching/"
# tif1 = dir + '092D_04877_111214_none.tif'
# outdir = "../manual_stitching/"
# combine2tifs(tif1, tif2, '_mode', "", outdir)