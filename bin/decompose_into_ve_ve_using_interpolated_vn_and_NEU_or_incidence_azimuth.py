#!/usr/bin/env python3
# -*- coding: utf-8 -*-

####################
# Decomposing ascending and descending LOS velocities into Ve and Vu by first removing contribution of interpolated
# GNSS Vn from the LOS.
# - lines 338-372: Change input parameters to read your own data files.
#
# Written by Qi Ou, University of Leeds, 27 Oct 2023
# email: q.ou@leeds.ac.uk
#
# If you use this script for your study, please cite:
# Ou, Q., Daout, S., Weiss, J. R., Shen, L., Lazecký, M., Wright, T. J., & Parsons, B. E. (2022). Large-scale interseismic strain mapping of the NE Tibetan Plateau from Sentinel-1 interferometry. Journal of Geophysical Research: Solid Earth, 127, e2022JB024176. https://doi.org/10.1029/2022JB024176
####################

import logging
# from merge_tif import *
import numpy as np
import os
import gdal
import copy
import matplotlib.pyplot as plt
from pathlib import Path
from cmcrameri import cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('decomposition.log')


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


def load_los_sigma_neu(group_dic, track_list):
    ##########################################################################
    # Loading input LOS, Sigma and coefs N, E, U into a dictionary of tracks.#
    # Manually define file directories and suffix here.                      #
    ##########################################################################

    for track in track_list:
        if NEU:
            group_dic[track] = OpenTif(los_dir + los_prefix + track + los_suffix,
                                       sigfile=sigma_dir + sigma_preffix + track + sigma_suffix,
                                       N=angle_dir + N_prefix + track + N_suffix,
                                       E=angle_dir + E_prefix + track + E_suffix,
                                       U=angle_dir + U_prefix + track + U_suffix)
        else:
            group_dic[track] = OpenTif(los_dir + los_prefix + track + los_suffix,
                                       sigfile=sigma_dir + sigma_preffix + track + sigma_suffix,
                                       incidence=angle_dir + inc_prefix + track + inc_suffix,
                                       heading=angle_dir + head_prefix + track + head_suffix)

        if flip_sign:
            group_dic[track].data = group_dic[track].data * -1

        if constant_los_sig:
            if track == 'A1' or track == 'A2':
                group_dic[track].sigma = ~np.isnan(group_dic[track].sigma) * asc_sig
            if track == 'D1' or track == 'D2':
                group_dic[track].sigma = ~np.isnan(group_dic[track].sigma) * dsc_sig


def export_tif(data, df, filename):
    # Export data to tif format.
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(filename, df.xsize, df.ysize, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([df.left, df.xres, 0, df.top, 0, df.yres])  ##sets same geotransform as input
    # outdata.SetProjection(df.projection)  ##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(data)
    outdata.FlushCache()
    outdata.FlushCache()  # need to flush twice to export the last tif properly, otherwise it stops halfway.


class Raster(object):
    """ a Class that stores the dimensional specification of a grid mesh."""
    def __init__(self, north=None, south=None, west=None, east=None, x_step=None, y_step=None, width=None, length=None):
        if west is not None:
            self.left = west
        if north is not None:
            self.top = north
        if east is not None:
            self.right = east
        if south is not None:
            self.bottom = south
        if x_step is not None:
            self.xres = x_step
        if y_step is not None:
            self.yres = y_step
        if width is not None:
            self.xsize = width
        else:
            self.xsize = int((self.right - self.left) / self.xres + 1.5)
        if length is not None:
            self.ysize = length
        else:
            self.ysize = int((self.bottom - self.top) / self.yres + 1.5)


def non_nan_merge(big_data, small_data, nodata_test, x_shift, y_shift, xsize, ysize):
    masked_data = np.choose(nodata_test,  # False = not nan = pick from 0th entry; True = nan = pick from 1st entry
        (small_data, big_data[y_shift:y_shift+ysize, x_shift:x_shift+xsize]))
    big_data[y_shift:y_shift+ysize, x_shift:x_shift+xsize] = masked_data


def remove_vn_contribution_from_los_and_propagate_sigma_n_to_sigma_los(gp, vn):
    for i in gp:
        vn_overlap = Overlap(gp[i], vn)
        gp[i].data = gp[i].data - vn_overlap.d2array * gp[i].N
        gp[i].sigma = np.sqrt(np.power(gp[i].sigma, 2) + np.power(vn_overlap.d2sigma * gp[i].N, 2))


def define_map_size(vn, asc_0, asc_1, dsc_0, dsc_1):
    """ Define size of the big raster grid, big enough to host all Vn and LOS tracks."""
    left = []
    right = []
    top = []
    bottom = []
    xres = []
    yres = []

    left.append(vn.left)
    right.append(vn.right)
    top.append(vn.top)
    bottom.append(vn.bottom)
    xres.append(vn.xres)
    yres.append(vn.yres)

    for group in [asc_0, asc_1, dsc_0, dsc_1]:
        for track in group:
            # print(track)
            left.append(group[track].left)
            right.append(group[track].right)
            top.append(group[track].top)
            bottom.append(group[track].bottom)
            xres.append(group[track].xres)
            yres.append(group[track].yres)

    if len(set(xres)) > 1 or len(set(yres)) > 1:
        raise Warning("tracks have different resolution")

    # define the geographic boundary of the map
    logger.info("define the geographic boundary of the map")
    big_map = Raster(north=max(top), south=min(bottom), west=min(left), east=max(right), x_step=xres[0], y_step=yres[0])
    return big_map


def focus_to_content(array, raster):
    """
    Zoom to non-nan region and give a map extent of the zoomed-in region
    :param array: 2D array
    :param raster: corresponding instance of the Raster class
    :return: sub-array and focused raster
    """
    nonnan_indices = np.argwhere(~np.isnan(array))
    x_min = np.min(nonnan_indices[0])
    x_max = np.max(nonnan_indices[0])
    y_min = np.min(nonnan_indices[1])
    y_max = np.max(nonnan_indices[1])
    west = raster.west + x_min * raster.x_step
    east = raster.east + x_min * raster.x_step


def input_track_los_sigma_E_U_to_map(gp, los_map, sig_map, E_map, U_map, maps):
    """ Add los, sig, E and U data from each group into four raster layers of size defined by maps. """
    print("input_track_los_sigma_E_U_to_map")
    for i in gp:
        print(i)
        x_shift = int((gp[i].left - maps.left) / maps.xres + 0.5)
        y_shift = int((gp[i].top - maps.top) / maps.yres + 0.5)
        #nodata_test = np.isnan(gp[i].data)  # = True if nan, = False if not nan; True = 1, False = 0
        non_nan_merge(los_map, gp[i].data, np.isnan(gp[i].data), x_shift, y_shift, gp[i].xsize, gp[i].ysize)
        non_nan_merge(sig_map, gp[i].sigma, np.isnan(gp[i].sigma), x_shift, y_shift, gp[i].xsize, gp[i].ysize)
        non_nan_merge(E_map, gp[i].E, np.isnan(gp[i].E), x_shift, y_shift, gp[i].xsize, gp[i].ysize)
        non_nan_merge(U_map, gp[i].U, np.isnan(gp[i].U), x_shift, y_shift, gp[i].xsize, gp[i].ysize)


if __name__ == "__main__":

    # #################################################
    # # Define output directories and file suffix here.#
    # #################################################
    output_dir = "/Users/qi/leeds/projects/insar/india_asia_results_mar24/"
    output_suffix = ""
    export_ve_vu = True
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    #################################################
    # Define vn, los, los_sigma directories and file formats here.#
    #################################################
    vn_dir = "/Users/qi/leeds/projects/insar/india_asia_results_mar24/velmap_results/"
    vn_file = vn_dir+"Vn.tif"
    sn_file = vn_dir+"Vn_std.tif"

    asc_0_list = ["A1"]
    asc_1_list = ["A2"]
    dsc_0_list = ["D1"]
    dsc_1_list = ["D2"]

    los_dir = "/Users/qi/leeds/projects/insar/india_asia_results_mar24/insar_los_mosaics/"
    los_prefix, los_suffix = "vel_eurasia_", ".tif"
    sigma_dir = "/Users/qi/leeds/projects/insar/india_asia_results_mar24/insar_los_mosaics/"
    sigma_preffix, sigma_suffix = "vstd_", ".tif"
    flip_sign = True  # True = input is from Velmap [+away from satellite]; False = input is from LiCSBAS [+towards satellite]

    #################################################
    # Define directory containing NEU coefs or incidence and heading angles in degrees
    #################################################
    angle_dir = "/Users/qi/leeds/projects/insar/india_asia_results_mar24/insar_los_mosaics/"
    NEU = False
    if NEU:
        N_prefix, N_suffix = "", "_N.tif"
        E_prefix, E_suffix = "", "_E.tif"
        U_prefix, U_suffix = "", "_U.tif"
    else:
        inc_prefix, inc_suffix = "incidence_", ".tif"
        head_prefix, head_suffix = "azimuth_", ".tif"

    ####################################
    # predict los based on one or more geometries of los turned off
    predict_los = False
    noA1 = False
    noA2 = False
    noD1 = False
    noD2 = False
    ####################################

    ####################################
    # Use constants as los sig or read LiCSBAS sigma outputs
    constant_los_sig = False
    if constant_los_sig:
        asc_sig = 1.7
        dsc_sig = 1.7
        output_suffix = output_suffix+"_asig{}_dsig{}".format(asc_sig, dsc_sig)
    ####################################


    ####################################
    # Start running
    ####################################

    # group tracks into staggered layers to avoid conflict
    asc_0 = {}
    asc_1 = {}
    dsc_0 = {}
    dsc_1 = {}
    print("loading asc_0...")
    load_los_sigma_neu(asc_0, asc_0_list)
    print("loading asc_1...")
    load_los_sigma_neu(asc_1, asc_1_list)
    print("loading dsc_0...")
    load_los_sigma_neu(dsc_0, dsc_0_list)
    print("loading dsc_1...")
    load_los_sigma_neu(dsc_1, dsc_1_list)

    # remove vn contribution from asc and dsc los and sigma(los)
    vn = OpenTif(vn_file, sigfile=sn_file)
    for group in [asc_0, asc_1, dsc_0, dsc_1]:
        print("remove vn from {}".format(group))
        remove_vn_contribution_from_los_and_propagate_sigma_n_to_sigma_los(group, vn)

    # define outer boundary of the data sets
    maps = define_map_size(vn, asc_0, asc_1, dsc_0, dsc_1)
    # empty_raster = np.ones((maps.ysize, maps.xsize)) * np.nan
    # # export_tif(empty_raster, maps, "empty.tif")
    # # empty = OpenTif("empty.tif")
    # #
    # x_shift = int((vn.left - maps.left) / maps.xres + 0.5)
    # y_shift = int((vn.top - maps.top) / maps.yres + 0.5)
    # # nodata_test = np.isnan(gp[i].data)  # = True if nan, = False if not nan; True = 1, False = 0
    # non_nan_merge(empty_raster, vn.data, np.isnan(vn.data), x_shift, y_shift, vn.xsize, vn.ysize)
    # export_tif(empty_raster, maps, "vn_qi.tif")
    #
    # sn = OpenTif(sn_file)
    # empty_raster = np.ones((maps.ysize, maps.xsize))
    # # non_nan_merge(empty_raster, sn.data, np.isnan(sn.data), x_shift, y_shift, sn.xsize, sn.ysize)
    # export_tif(empty_raster, maps, "vstd_qi.tif")
    # plt.imshow(empty_raster)
    # plt.show()
    # breakpoint()

    # create 4x4 layers of empty map arrays to host input
    logger.info("create 4x4 layers of empty map arrays to host input")
    asc_los_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    asc_los_1 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_los_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_los_1 = np.ones((maps.ysize, maps.xsize)) * np.nan
    asc_sig_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    asc_sig_1 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_sig_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_sig_1 = np.ones((maps.ysize, maps.xsize)) * np.nan
    asc_E_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    asc_E_1 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_E_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_E_1 = np.ones((maps.ysize, maps.xsize)) * np.nan
    asc_U_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    asc_U_1 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_U_0 = np.ones((maps.ysize, maps.xsize)) * np.nan
    dsc_U_1 = np.ones((maps.ysize, maps.xsize)) * np.nan

    # place groups of inputs in the right positions in the 4x4 layers of maps.
    logger.info("place groups of inputs in the right positions in the 4x4 layers of maps.")
    input_track_los_sigma_E_U_to_map(asc_0, asc_los_0, asc_sig_0, asc_E_0, asc_U_0, maps)
    input_track_los_sigma_E_U_to_map(asc_1, asc_los_1, asc_sig_1, asc_E_1, asc_U_1, maps)
    input_track_los_sigma_E_U_to_map(dsc_0, dsc_los_0, dsc_sig_0, dsc_E_0, dsc_U_0, maps)
    input_track_los_sigma_E_U_to_map(dsc_1, dsc_los_1, dsc_sig_1, dsc_E_1, dsc_U_1, maps)

    if predict_los:
        # turn one layer to nan
        if noA1:
            output_suffix = output_suffix + "_noA1"
            a0los = copy.copy(asc_los_0)
            asc_los_0 = asc_los_0 * np.nan
        if noA2:
            output_suffix = output_suffix + "_noA2"
            a1los = copy.copy(asc_los_1)
            asc_los_1 = asc_los_1 * np.nan
        if noD1:
            output_suffix = output_suffix + "_noD1"
            d0los = copy.copy(dsc_los_0)
            dsc_los_0 = dsc_los_0 * np.nan
        if noD2:
            output_suffix = output_suffix + "_noD2"
            d1los = copy.copy(dsc_los_1)
            dsc_los_1 = dsc_los_1 * np.nan

    print("plotting...")
    for map in [asc_los_0, asc_los_1, dsc_los_0, dsc_los_1, asc_sig_0, asc_sig_1, dsc_sig_0, dsc_sig_1,
                asc_E_0, asc_E_1, dsc_E_0, dsc_E_1, asc_U_0, asc_U_1, dsc_U_0, dsc_U_1]:
        plt.imshow(map, vmin=np.nanpercentile(map, 10), vmax=np.nanpercentile(map, 90))
        plt.colorbar()
        plt.show()
    breakpoint()

    # flatten the 2D maps into 1D arrays for inversion
    asc_los_0_flat = asc_los_0.flatten()
    asc_los_1_flat = asc_los_1.flatten()
    dsc_los_0_flat = dsc_los_0.flatten()
    dsc_los_1_flat = dsc_los_1.flatten()
    # mask only pixels with at least one asc and at least one dsc in the 4 layers.
    mask_test = np.logical_and(np.logical_or(~np.isnan(asc_los_0_flat), ~np.isnan(asc_los_1_flat)),
                               np.logical_or(~np.isnan(dsc_los_0_flat), ~np.isnan(dsc_los_1_flat)))
    # select elements from the 4x4 layers which have enough data to enter inversion.
    asc_los_0_mask = asc_los_0_flat[mask_test]
    asc_los_1_mask = asc_los_1_flat[mask_test]
    dsc_los_0_mask = dsc_los_0_flat[mask_test]
    dsc_los_1_mask = dsc_los_1_flat[mask_test]
    asc_sig_0_mask = asc_sig_0.flatten()[mask_test]
    asc_sig_1_mask = asc_sig_1.flatten()[mask_test]
    dsc_sig_0_mask = dsc_sig_0.flatten()[mask_test]
    dsc_sig_1_mask = dsc_sig_1.flatten()[mask_test]
    asc_E_0_mask = asc_E_0.flatten()[mask_test]
    asc_E_1_mask = asc_E_1.flatten()[mask_test]
    dsc_E_0_mask = dsc_E_0.flatten()[mask_test]
    dsc_E_1_mask = dsc_E_1.flatten()[mask_test]
    asc_U_0_mask = asc_U_0.flatten()[mask_test]
    asc_U_1_mask = asc_U_1.flatten()[mask_test]
    dsc_U_0_mask = dsc_U_0.flatten()[mask_test]
    dsc_U_1_mask = dsc_U_1.flatten()[mask_test]

    # empty lists to register inversion results
    ve_mask = []
    vu_mask = []
    ve_sig_mask = []
    vu_sig_mask = []

    # LOS = Ve * E + Vu * U

    # # inversion starts
    length = len(asc_los_0_mask)
    logger.info("inversion starts...total {} pixels".format(length))
    c = 1
    for la0, la1, ld0, ld1, lasig0, lasig1, ldsig0, ldsig1, \
        Gae0, Gae1, Gde0, Gde1, Gau0, Gau1, Gdu0, Gdu1 in \
            zip(asc_los_0_mask, asc_los_1_mask, dsc_los_0_mask, dsc_los_1_mask,
                asc_sig_0_mask, asc_sig_1_mask, dsc_sig_0_mask, dsc_sig_1_mask,
                asc_E_0_mask,   asc_E_1_mask,   dsc_E_0_mask,   dsc_E_1_mask,
                asc_U_0_mask,   asc_U_1_mask,   dsc_U_0_mask,   dsc_U_1_mask):

        d = np.array([[la0], [la1], [ld0], [ld1]])
        G = np.array([[Gae0, Gau0], [Gae1, Gau1], [Gde0, Gdu0], [Gde1, Gdu1]])
        sig = np.array([[lasig0], [lasig1], [ldsig0], [ldsig1]])
        sig[np.isnan(sig)] = 1 # in case there are empty sigma values where velocity have values
        sig = sig[~np.isnan(d).any(axis=1)]
        G = G[~np.isnan(d).any(axis=1)]
        d = d[~np.isnan(d).any(axis=1)]
        if len(d) != len(G) or len(d) != len(sig):
            raise Warning("pixel inversion inputs have unequal dimensions...")
        [[e], [u]] = np.linalg.lstsq(G/sig, d/sig, rcond=None)[0]
        cov_d = np.diag(np.square(sig).transpose()[0])
        cov_m = np.linalg.inv(np.dot(np.dot(G.transpose(), np.linalg.inv(cov_d)), G))
        if c % (length//100) == 0:
            logger.info("{}, {:.2f}%, {}, {}".format(c, c/length*100, e, u))
        ve_mask.append(e)
        vu_mask.append(u)
        ve_sig_mask.append(np.sqrt(cov_m[0, 0]))
        vu_sig_mask.append(np.sqrt(cov_m[1, 1]))
        c = c+1

    print('done with pixel wise inversion')

    # create empty 1D arrays to register inversion results
    length = len(asc_los_0_flat)
    ve_flat = np.ones(length) * np.nan
    vu_flat = np.ones(length) * np.nan
    ve_sig_flat = np.ones(length) * np.nan
    vu_sig_flat = np.ones(length) * np.nan
    # populate the masked pixels of the full flat arrays with inversion results
    ve_flat[mask_test] = ve_mask
    vu_flat[mask_test] = vu_mask
    ve_sig_flat[mask_test] = ve_sig_mask
    vu_sig_flat[mask_test] = vu_sig_mask
    # reshape the flat arrays into maps
    ve = ve_flat.reshape((maps.ysize, maps.xsize))
    vu = vu_flat.reshape((maps.ysize, maps.xsize))
    ve_sig = ve_sig_flat.reshape((maps.ysize, maps.xsize))
    vu_sig = vu_sig_flat.reshape((maps.ysize, maps.xsize))

    if predict_los:
        # plotting los differences and uncertainty estimates
        fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, sharex='all', sharey='all')
        plt.suptitle("LOS_diff"+output_suffix)
        # back calculate removed los and compare differences
        if noA1:
            asc_los_0 = ve * asc_E_0 + vu * asc_U_0
            asc_los_0_diff = a0los - asc_los_0
            asc_los_0_std = np.nanstd(asc_los_0_diff)
            # plot
            vmin = np.nanpercentile(asc_los_0_diff, 1)
            vmax = np.nanpercentile(asc_los_0_diff, 99)
            im = ax1.imshow(asc_los_0_diff, vmin=vmin, vmax=vmax)
            c = plt.colorbar(im, ax=ax1)
            c.set_label('mm/yr')
            ax1.set_title('A1 diff (std = {:.1f} mm/yr)'.format(asc_los_0_std))
            # export
            export_tif(asc_los_0_diff, maps, output_dir + 'a1_diff{}.tif'.format(output_suffix))

        if noA2:
            asc_los_1 = ve * asc_E_1 + vu * asc_U_1
            asc_los_1_diff = a1los - asc_los_1
            asc_los_1_std = np.nanstd(asc_los_1_diff)
            # plot
            vmin = np.nanpercentile(asc_los_1_diff, 1)
            vmax = np.nanpercentile(asc_los_1_diff, 99)
            im = ax2.imshow(asc_los_1_diff, vmin=vmin, vmax=vmax)
            c = plt.colorbar(im, ax=ax2)
            c.set_label('mm/yr')
            ax2.set_title('A2 diff (std = {:.1f} mm/yr)'.format(asc_los_1_std))
            # export
            export_tif(asc_los_1_diff, maps, output_dir + 'a2_diff{}.tif'.format(output_suffix))

        if noD1:
            dsc_los_0 = ve * dsc_E_0 + vu * dsc_U_0
            dsc_los_0_diff = d0los - dsc_los_0
            dsc_los_0_std = np.nanstd(dsc_los_0_diff)
            # plot
            vmin = np.nanpercentile(dsc_los_0_diff, 1)
            vmax = np.nanpercentile(dsc_los_0_diff, 99)
            im = ax3.imshow(dsc_los_0_diff, vmin=vmin, vmax=vmax)
            c = plt.colorbar(im, ax=ax3)
            c.set_label('mm/yr')
            ax3.set_title('D1 diff (std = {:.1f} mm/yr)'.format(dsc_los_0_std))
            # export
            export_tif(dsc_los_0_diff, maps, output_dir + 'd1_diff{}.tif'.format(output_suffix))

        if noD2:
            dsc_los_1 = ve * dsc_E_1 + vu * dsc_U_1
            dsc_los_1_diff = d1los - dsc_los_1
            dsc_los_1_std = np.nanstd(dsc_los_1_diff)
            # plot
            vmin = np.nanpercentile(dsc_los_1_diff, 1)
            vmax = np.nanpercentile(dsc_los_1_diff, 99)
            im = ax4.imshow(dsc_los_1_diff, vmin=vmin, vmax=vmax)
            c = plt.colorbar(im, ax=ax4)
            c.set_label('mm/yr')
            ax4.set_title('D2 diff (std = {:.1f} mm/yr)'.format(dsc_los_1_std))
            # export
            export_tif(dsc_los_1_diff, maps, output_dir + 'd2_diff{}.tif'.format(output_suffix))

        plt.show()
        fig.savefig(output_dir+"LOS_diff"+output_suffix+".png", format='PNG', dpi=300, bbox_inches='tight')

    if export_ve_vu:
        export_tif(ve, maps, output_dir+'ve{}.tif'.format(output_suffix))
        export_tif(vu, maps, output_dir+'vu{}.tif'.format(output_suffix))
        export_tif(ve_sig, maps, output_dir+'ve_sig{}.tif'.format(output_suffix))
        export_tif(vu_sig, maps, output_dir+'vu_sig{}.tif'.format(output_suffix))

        # plotting
        ve_vmin = np.nanpercentile(ve, 1)
        ve_vmax = np.nanpercentile(ve, 99)
        vu_vmin = np.nanpercentile(vu, 1)
        vu_vmax = np.nanpercentile(vu, 99)
        ve_sig_vmin = np.nanpercentile(ve_sig, 1)
        ve_sig_vmax = np.nanpercentile(ve_sig, 99)
        vu_sig_vmin = np.nanpercentile(vu_sig, 1)
        vu_sig_vmax = np.nanpercentile(vu_sig, 99)

        # ve_vmin=0
        # ve_vmax=20
        vu_vmin=-5
        vu_vmax=5
        ve_sig_vmin=0
        ve_sig_vmax=1
        vu_sig_vmin=0
        vu_sig_vmax=1

        fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, sharex='all', sharey='all', figsize=(6.8, 4))
        im = ax1.imshow(ve, vmin=ve_vmin, vmax=ve_vmax, cmap=cm.roma.reversed(), extent=(maps.left, maps.right, maps.bottom, maps.top))
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="3%", pad="2%")
        c = plt.colorbar(im, cax=cax)
        c.set_label('mm/yr')
        ax1.set_title('Ve')
        im = ax2.imshow(vu, vmin=vu_vmin, vmax=vu_vmax, cmap=cm.roma.reversed(), extent=(maps.left, maps.right, maps.bottom, maps.top))
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="3%", pad="2%")
        c = plt.colorbar(im, cax=cax)
        c.set_label('mm/yr')
        ax2.set_title('Vu')
        im = ax3.imshow(ve_sig, vmin=ve_sig_vmin, vmax=ve_sig_vmax, extent=(maps.left, maps.right, maps.bottom, maps.top))
        divider = make_axes_locatable(ax3)
        cax = divider.append_axes("right", size="3%", pad="2%")
        c = plt.colorbar(im, cax=cax)
        c.set_label('mm/yr')
        ax3.set_title(r'$\sigma$(Ve)')
        im = ax4.imshow(vu_sig, vmin=vu_sig_vmin, vmax=vu_sig_vmax, extent=(maps.left, maps.right, maps.bottom, maps.top))
        divider = make_axes_locatable(ax4)
        cax = divider.append_axes("right", size="3%", pad="2%")
        c = plt.colorbar(im, cax=cax)
        c.set_label('mm/yr')
        ax4.set_title(r'$\sigma$(Vu)')
        plt.show()
        fig.savefig(output_dir+'ve_vu{}.png'.format(output_suffix), format='PNG', dpi=300, bbox_inches='tight')

















