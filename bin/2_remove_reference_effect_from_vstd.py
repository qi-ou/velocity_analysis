#!/usr/bin/env python3
# -*- coding: utf-8 -*-

####################
# Removing the reference effect from the line-of-sight uncertainty maps by fitting a spherical model to the
# scatter between uncertainty and distance away from the reference center.
# Written by Qi Ou, University of Oxford, 1 May 2022
# email: qi.ou@earth.ox.ac.uk
####################
from merge_tif import *
from pylab import *
# from invert_gps_overlap import *
# gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
# shapely.speedups.enable()
import pyproj
from scipy import stats
import random
from lmfit.model import *
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


# define global variable
plot_variogram = True


def spherical(d, p, n, r):
    """
    Compute spherical variogram model
    @param d: 1D distance array
    @param p: partial sill
    @param n: nugget
    @param r: range
    @return: spherical variogram model
    """
    if r>d.max():
        r=d.max()-1
    return np.where(d > r, p + n, p * (3/2 * d/r - 1/2 * d**3 / r**3) + n)


def fit_model(model, dat, x, sigma):
    # fitting data with variogram model
    model.set_param_hint('p', value=dat[-1])  # guess last dat
    model.set_param_hint('n', value=dat[0])  # guess first dat
    model.set_param_hint('r', value=x[len(x)//2])  # guess mid point distance
    out = model.fit(dat, d=x, weights=1/sigma)
    return out


def exponential(d, p, n, r):
    """
    Compute exponential variogram model
    @param d: 1D distance array
    @param p: partial sill
    @param n: nugget
    @param r: range
    @return: spherical variogram model
    """
    return n + p * (1 - np.exp(-d * 3/r))


def define_UTM(opentif):
    # Define UTM zone based on the centre of the OpenTif object
    if 96 <= (opentif.left + opentif.right) / 2 < 102:
        utm = pyproj.Proj(proj='utm', zone=47)
    elif 102 <= (opentif.left + opentif.right) / 2 < 108:
        utm = pyproj.Proj(proj='utm', zone=48)
    elif 90 <= (opentif.left + opentif.right) / 2 < 96:
        utm = pyproj.Proj(proj='utm', zone=46)
    elif 84 <= (opentif.left + opentif.right) / 2 < 90:
        utm = pyproj.Proj(proj='utm', zone=45)
    elif 78 <= (opentif.left + opentif.right) / 2 < 84:
        utm = pyproj.Proj(proj='utm', zone=44)
    elif 72 <= (opentif.left + opentif.right) / 2 < 78:
        utm = pyproj.Proj(proj='utm', zone=43)
    elif 66 <= (opentif.left + opentif.right) / 2 < 72:
        utm = pyproj.Proj(proj='utm', zone=42)
    elif 60 <= (opentif.left + opentif.right) / 2 < 66:
        utm = pyproj.Proj(proj='utm', zone=41)
    return utm


def find_median_index_of_lowest_values(opentif):
    """
    Define the reference point of an uncertainty map based on the median indices of pixels with lowest 2% uncertainties
    @param opentif: OpenTif Class object
    @return: add a reference location as attribute to the OpenTif Class object
    """
    dat = opentif.data
    # get the indices of all pixels with lowest sigma
    ref_locs = np.nonzero(dat < np.nanpercentile(dat, 1))
    # get the median row, col index for the "centre of the reference window"
    ref_loc = np.median(ref_locs, axis=1)
    return ref_loc


def get_profile_data(tf):
    """
    Find reference point and turn sigma map to uncertainty/distance profiles relative to reference location
    @param tf: OpenTif class object
    @return: nonnan_mask, a sigma 1D array and a distance (km) 1D array for nonnan pixels
    """
    tf.ref_loc = find_median_index_of_lowest_values(tf)
    ref_lat, ref_lon = tf.top + tf.yres * tf.ref_loc[0], tf.left + tf.xres * tf.ref_loc[1]

    # get UTM coordinates of the frame and the reference pixel
    UTM = define_UTM(tf)
    ref_east, ref_north = UTM(ref_lon, ref_lat)
    east, north = UTM(tif.lon, tif.lat)
    # print(tf.basename, tf.ref_loc, ref_east, ref_north)
    print(tf.basename)

    # calc 1D arrays of uncertainties, easting, northing, distance for non-nan pixels
    sig = tf.data.flatten()
    e = east.flatten()
    n = north.flatten()
    mask = ~np.isnan(sig)
    dat_nonnan = sig[mask]
    e_nonnan = e[mask]
    n_nonnan = n[mask]
    dist_nonnan = np.sqrt((e_nonnan - ref_east) ** 2 + (n_nonnan - ref_north) ** 2) / 1000
    # print(dist_nonnan)

    return mask, dat_nonnan, dist_nonnan


def calc_median_std(x, y, bin_num=50):
    """
    Calculate binned statistics (median and std) of y scatters for 100 bins along x axis.
    @param bin_num: number of bins for median/std calculations, default to 100
    @param x: 1D array
    @param y: 1D array
    @return: 1D arrays of median, std and bincenters
    """
    # calc statistics of 1D scatter plot
    median_array, binedges = stats.binned_statistic(x, y, 'median', bins=bin_num)[:-1]
    std_array = stats.binned_statistic(x, y, 'std', bins=bin_num)[0]
    bincenter_array = (binedges[0:-1] + binedges[1:]) / 2
    return median_array, std_array, bincenter_array


def scale_value_by_variogram_ratio(y, x, model_result, model):
    """
    Scale values of y by theoretical ratio of y at x and y at sill
    @param y: 1D array
    @param x: 1D array
    @param model_result: ModelResult class object from fitting variogram model
    @param model: string describing model type, used to toggle between forward model calculations
    @return: scaled array of y
    """
    best = model_result.best_values
    sill = best['p'] + best['n']  # total sill
    if model == 'spherical':
        theoretical_y = spherical(x, best['p'], best['n'], best['r'])
    if model == 'exponential':
        theoretical_y = exponential(x, best['p'], best['n'], best['r'])
    scaling_factor = sill / theoretical_y
    y_scaled = y * scaling_factor
    return y_scaled


if __name__ == "__main__":
    ###################
    # input parameters:
    ###################
    # sigma_suffix = '*.fringe.vstd.tif'
    # sigma_dir = '../vstd/'
    sigma_suffix = '129A_04939_131313.weighted.vstd.tif'  # regular expression for vstd.tif
    sigma_dir = '../los_weighted/vstd'
    output_dir = '../los_weighted/vstd_corrected/'
    output_suffix = '_scaled_vstd.tif'
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    sigfileList = sorted(glob.glob(os.path.join(sigma_dir, sigma_suffix)))
    for f in sigfileList:
        tif = OpenTif(f)
        try:
            # prepare profile data
            nonnan_mask, sig_nonnan, dist = get_profile_data(tif)
            # calc stats along profile
            median, std, bincenters = calc_median_std(dist, sig_nonnan)
            # fit variogram model to stats along profile
            spherical_model = Model(spherical)
            exponential_model = Model(exponential)

            try:
                # take smaller range
                mask = bincenters < 150
                bincenters = bincenters[mask]
                std = std[mask]
                median = median[mask]
                # weight points nearer to ref more heavily
                result_spherical = fit_model(spherical_model, median, bincenters, std+np.power(bincenters/max(bincenters), 3))
                result_exponential = fit_model(exponential_model, median, bincenters, std+np.power(bincenters/max(bincenters), 3))
            except:
                try:
                    # take smaller range
                    mask = bincenters < 120
                    bincenters = bincenters[mask]
                    std = std[mask]
                    median = median[mask]

                    # weight points nearer to ref more heavily
                    result_spherical = fit_model(spherical_model, median, bincenters, std+np.power(bincenters/max(bincenters),3))
                    result_exponential = fit_model(exponential_model, median, bincenters, std+np.power(bincenters/max(bincenters),3))
                except:
                    # take smaller range
                    mask = bincenters < 100
                    bincenters = bincenters[mask]
                    std = std[mask]
                    median = median[mask]

                    # weight points nearer to ref more heavily
                    result_spherical = fit_model(spherical_model, median, bincenters, std+np.power(bincenters/max(bincenters),3))
                    result_exponential = fit_model(exponential_model, median, bincenters, std+np.power(bincenters/max(bincenters),3))

            if result_spherical.chisqr < result_exponential.chisqr:
                result = result_spherical
                model = 'spherical'
                print("Choosing spherical model")
            else:
                result = result_exponential
                model = 'exponential'
                print("Choosing exponential model")

            # scale sig_nonnan based on variogram model
            sig_nonnan_scaled = scale_value_by_variogram_ratio(sig_nonnan, dist, result, model)
            # populate scalled uncertainty map
            new_map = np.ones(nonnan_mask.shape)*np.nan
            new_map[nonnan_mask] = sig_nonnan_scaled
            new_map = new_map.reshape(tif.data.shape)
            #
            # median_scaled = scale_value_by_variogram_ratio(median, bincenters, result)
            # median_plus_scaled = scale_value_by_variogram_ratio(median+std, bincenters, result)
            # median_minus_scaled = scale_value_by_variogram_ratio(median-std, bincenters, result)
            fit_scaled = scale_value_by_variogram_ratio(result.best_fit, bincenters, result, model)

            # plotting
            fig, [[ax1, ax2], [ax3, ax4]] = subplots(2, 2, figsize=(6.4, 4.8))

            # plot scatter with stats of uncertainty values
            subset = random.sample(range(1, len(dist)), min(int(len(dist)/2), 50000))
            ax1.scatter(dist[subset], sig_nonnan[subset], c=sig_nonnan[subset], s=0.1, vmin=0, vmax=2)
            ax1.plot(bincenters, median, linewidth=2, c="gold")
            ax1.plot(bincenters, median+std, linewidth=1, c="gold")
            ax1.plot(bincenters, median-std, linewidth=1, c="gold")
            ax1.set_xlim((0, bincenters[-1]))
            ax1.set_ylim((0, min(2, max(median+3*std))))
            ax1.tick_params(labelbottom=False)
            ax1.set_ylabel("$\sigma$(LOS), mm/yr")
            ax1.set_title("Uncertainty Profile")
            ax1.plot(bincenters, result.best_fit, linewidth=2, c="red", label='n={:.1f}, s={:.1f}, r={:.0f}'
                     .format(result.best_values['n'], result.best_values['p']+result.best_values['n'], result.best_values['r']))
            ax1.legend(loc=4)


            # plot scatter with stats of uncertainty values
            ax3.scatter(dist[subset], sig_nonnan_scaled[subset], c=sig_nonnan_scaled[subset], s=0.1, vmin=0, vmax=2)
            # ax3.plot(bincenters, median_scaled, linewidth=2, c="gold")
            # ax3.plot(bincenters, median_plus_scaled, linewidth=1, c="gold")
            # ax3.plot(bincenters, median_minus_scaled, linewidth=1, c="gold")
            ax3.plot(bincenters, fit_scaled, linewidth=2, c="red")
            ax3.set_ylim((0, min(2, max(median+3*std))))
            ax3.set_xlabel("Distance to reference, km")
            ax3.set_ylabel("$\sigma$(LOS), mm/yr")
            ax3.set_title("Scaled Uncertainty Profile")
            ax3.set_xlim((0, bincenters[-1]))

            # plot uncertainty map with reference location
            im = ax2.imshow(tif.data, vmin=0, vmax=2)
            ax2.plot(tif.ref_loc[1], tif.ref_loc[0], marker="o", markersize=5, c='gold')
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, label="mm/yr")
            ax2.set_title("Uncertainty Map")

            # plot scaled uncertainty map
            im = ax4.imshow(new_map, vmin=0, vmax=2)
            divider = make_axes_locatable(ax4)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, label="mm/yr")
            ax4.set_title("Scaled Uncertainty Map")
            ax4.set_xlabel(tif.basename[:17], labelpad=15)

            for ax in ax2, ax4:
                ax.set_xticks([])
                ax.set_yticks([])
            plt.tight_layout()
            plt.show()
            fig.savefig(output_dir+tif.basename+".png", format='PNG', dpi=300, bbox_inches='tight')

            # Export scaled uncertainty map to tif format.
            driver = gdal.GetDriverByName("GTiff")
            outdata = driver.Create(output_dir + tif.basename + output_suffix, tif.xsize, tif.ysize, 1, gdal.GDT_Float32)
            outdata.SetGeoTransform([tif.left, tif.xres, 0, tif.top, 0, tif.yres])  ##sets same geotransform as input
            outdata.SetProjection(tif.projection)  ##sets same projection as input
            outdata.GetRasterBand(1).WriteArray(new_map)
            outdata.FlushCache()
            outdata.FlushCache()  # need to flush twice to export the last tif properly, otherwise it stops halfway.

        except:
            print("entering exception, plot before scaling only")
            # prepare profile data
            nonnan_mask, sig_nonnan, dist = get_profile_data(tif)
            # calc stats along profile
            median, std, bincenters = calc_median_std(dist, sig_nonnan)

            # plotting
            fig, [[ax1, ax2], [ax3, ax4]] = subplots(2, 2, figsize=(6.4, 4.8))

            # plot scatter with stats of uncertainty values
            subset = random.sample(range(1, len(dist)), min(int(len(dist)/2), 50000))
            ax1.scatter(dist[subset], sig_nonnan[subset], c=sig_nonnan[subset], s=0.1, vmin=0, vmax=2)
            ax1.plot(bincenters, median, linewidth=2, c="gold")
            ax1.plot(bincenters, median+std, linewidth=1, c="gold")
            ax1.plot(bincenters, median-std, linewidth=1, c="gold")
            ax1.set_xlim((0, bincenters[-1]))
            ax1.set_ylim((0, min(2, max(median+3*std))))
            ax1.tick_params(labelbottom=False)
            ax1.set_ylabel("$\sigma$(LOS), mm/yr")
            ax1.set_title("Uncertainty Profile")
            ax1.plot(bincenters, result.best_fit, linewidth=2, c="red", label='n={:.1f}, s={:.1f}, r={:.0f}'
                     .format(result.best_values['n'], result.best_values['p']+result.best_values['n'], result.best_values['r']))
            ax1.legend(loc=4)

            # plot uncertainty map with reference location
            im = ax2.imshow(tif.data, vmin=0, vmax=2)
            ax2.plot(tif.ref_loc[1], tif.ref_loc[0], marker="o", markersize=5, c='gold')
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, label="mm/yr")
            ax2.set_title("Uncertainty Map")

            for ax in [ax2]:
                ax.set_xticks([])
                ax.set_yticks([])
            plt.tight_layout()
            plt.show()


