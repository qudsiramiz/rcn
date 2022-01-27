# This the python version of IDL code named 'RX_model_batch.pro'
import datetime
import time
import warnings

import geopack.geopack as gp
import h5py as hf
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pyspedas as spd
import pytplot as ptt
from dateutil import parser
import scipy as sp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import meijering, sato, frangi, hessian
import multiprocessing as mp
from tabulate import tabulate

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

start = time.time()

today_date = datetime.datetime.today().strftime('%Y-%m-%d')


def get_shear(b_vec_1, b_vec_2, angle_unit="radians"):
    r"""
    Get the shear angle between two magnetic field lines.

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        Input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Input magnetic field vector.
    angle_unit : str, optional
        Preferred unit of angle returned by the code. Default is "radians".

    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
    angle: float
        Angle between the two vectors, in radians by default
    """
    unit_vec_1 = b_vec_1/np.linalg.norm(b_vec_1)
    unit_vec_2 = b_vec_2/np.linalg.norm(b_vec_2)
    angle = np.arccos(np.dot(unit_vec_1, unit_vec_2))

    #dp = b_vec_1[0] * b_vec_2[0] + b_vec_1[1] * b_vec_2[1] + b_vec_1[2] * b_vec_2[2]

    #mag1 = np.linalg.norm( b_vec_1)
    #mag2 = np.linalg.norm( b_vec_2)
    #angle = np.arccos(dp/(mag1*mag2))

    if (angle_unit == "radians"):
        return angle
    elif (angle_unit == "degrees"):
        return angle * 180/np.pi
    else:
        raise KeyError("angle_unit must be radians or degrees")


def get_rxben(b_vec_1, b_vec_2):
    r"""
    Get the reconnection energy between two magnetic field lines.

    It has the following mathematical expression:

    .. math:: rexben = 0.5 (|\vec{B_1}| + \vec{B_2}) (1 - \hat{B_1} \cdot \hat{B_2})

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        Input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Input magnetic field vector.

    Returns
    -------
    rxben : float
        Reconnection field energy density in nPa
    """

    alpha = - 14.87 * np.pi / 180  # radians (From Hesse2013)
    b_vec_1 = np.array(b_vec_1)
    b_vec_2 = np.array(b_vec_2)
    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    # The bisector vector of thwo input vectors
    unit_vec_bisec = (unit_vec_1 + unit_vec_2) / np.linalg.norm(unit_vec_1 + unit_vec_2)

    # Cross product of the two input vectors with the bisector vector to get the reconnection
    # component of the magnetic field.
    rx_b_1 = np.cross(b_vec_1, unit_vec_bisec)
    rx_b_2 = np.cross(b_vec_2, unit_vec_bisec)

    rx_b_mag_1 = np.linalg.norm(rx_b_1)
    rx_b_mag_2 = np.linalg.norm(rx_b_2)

    # The guide field of the reconnection
    gd_b_1 = b_vec_1 - rx_b_1
    gd_b_2 = b_vec_2 - rx_b_2

    gd_b = gd_b_1 + gd_b_2
    gd_b_mag = np.linalg.norm(gd_b)

    # The reconnection energy
    b_u_prime = rx_b_mag_1 * np.cos(alpha) + gd_b_mag * np.sin(alpha)
    b_d_prime = rx_b_mag_2 * np.cos(alpha) + gd_b_mag * np.sin(alpha)

    # Reconnection energy
    #rx_en = b_u_prime ** 2 * b_d_prime ** 2
    rx_en = rx_b_mag_1 ** 2 * rx_b_mag_2 ** 2
    #unit_vec_1 = b_vec_1/mag_vec_1
    #unit_vec_2 = b_vec_2/mag_vec_2
#
    #u_bisect = (unit_vec_1 + unit_vec_2) / 2
    #rx_bmag1 = np.dot(u_bisect, b_vec_1)
    #rx_bmag2 = np.dot(u_bisect, b_vec_2)
    #b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)
#
    ##rx_en = 0.5 * (mag_vec_1 * mag_vec_2) * (1 + b1_b2_dotp)
    #rx_en = (rx_bmag1**2 * rx_bmag2**2)   # nPa
    ##rx_en = (rx_bmag1**2 + rx_bmag2**2) * 1.03  # MJ/RE^3

    # Angle between the two vectors
    #angle_b1_b2 = get_shear(b_vec_1, b_vec_2, angle_unit="radians")
#
    #angle_bisect = angle_b1_b2/2.
#
    ## Reconnecting component of the magnetic field
    #rx_bmag1 = mag_vec_1 * np.sin(angle_bisect)
    #rx_bmag2 = mag_vec_2 * np.sin(angle_bisect)
#
    ## Non-reconnecting/guide field component of the magnetic field
    #gd_bmag1 = mag_vec_1 * np.cos(angle_bisect)
    #gd_bmag2 = mag_vec_2 * np.cos(angle_bisect)
#
    ##TODO: Check if this is correct
    ##gd_mag = np.sqrt(gd_bmag1**2 + gd_bmag2**2)  # Magnitude of guide field?
    #gd_mag = np.sqrt(gd_bmag1**2 + gd_bmag2**2)
    #alpha = 14.87 * np.pi/180.
#
    #b_u_prime = rx_bmag1 * np.cos(alpha)# + gd_mag * np.sin(alpha)
    #b_d_prime = rx_bmag2 * np.cos(alpha)# + gd_mag * np.sin(alpha)
#
    #rx_en = b_u_prime**2 * b_d_prime**2

    #angle_bisect = np.arccos(np.dot(b_vec_1, b_vec_2)/(np.linalg.norm(b_vec_1)*np.linalg.norm(b_vec_2)))

    return rx_en


def get_vcs(b_vec_1, b_vec_2, n_1, n_2):
    r"""
    Get vcs code.

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        Input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Input magnetic field vector.
    n_1 : float
        Density of first component
    n_2 : float
        Density of second component

    Returns
    -------
    vcs : float
        The exhaust velocity in km/s
    """
    va_p1 = 21.812  # conv. nT, m_P/cm ^ 3 product to km/s cassak-shay

    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)

    # bisector = mag_vec_1*b_vec_1 + mag_vec_2*b_vec_2 u_bisect = bisector /
    # np.linalg.norm(bisector)

    rx_mag_1 = mag_vec_1 * (1 + b1_b2_dotp)/2
    rx_mag_2 = mag_vec_2 * (1 + b1_b2_dotp)/2

    vcs = va_p1 * np.sqrt(rx_mag_1 * rx_mag_2 * (rx_mag_1 + rx_mag_2)/(rx_mag_1 * n_2 +
                                                                       rx_mag_2 * n_1))

    # vcs = va_p1 * np.sqrt(mag_vec_1 * mag_vec_2 * (mag_vec_1 + mag_vec_2)/(mag_vec_1 * n_2 +
    #                                                                         mag_vec_2 * n_1))

    return vcs


def get_bis(b_vec_1, b_vec_2):
    r"""
    Get the shear angle between two magnetic field lines.

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        First input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Second input magnetic field vector.

    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
    bis_field_1 : float
        The magnitude of bisected field line corresponding to the first input magnetic field vector.

    bis_field_2 : float
        The magnitude of bisected field line corresponding to the second input magnetic field vector.
    """
    b_vec_1 = np.array(b_vec_1)
    b_vec_2 = np.array(b_vec_2)

    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    unit_vec_bisect = (unit_vec_1 + unit_vec_2) / np.linalg.norm(unit_vec_1 + unit_vec_2)

    # Cross product of the two vectors with the bisector
    b_vec_1_cross = np.cross(b_vec_1, unit_vec_bisect)
    b_vec_2_cross = np.cross(b_vec_2, unit_vec_bisect)

    # Magnitude of the cross product vectors (should be the same for both for a symmetric
    # reconnection). These magnitude are the reconnecting components of the magnetic field.
    bis_field_1 = np.linalg.norm(b_vec_1_cross)
    bis_field_2 = np.linalg.norm(b_vec_2_cross)
    return bis_field_1, bis_field_2


def get_ca(b_vec, angle_unit="radians"):
    r"""
    Get ca.

    Parameters
    ----------
    b_vec : array of shape 1x3
        Input magnetic field vector.
    angle_unit : str, optional
        Preferred unit of angle returned by the code. Default is "radians".

    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
    angle : float Returns arctan of y- and z-component of the magnetic field vector.

    """
    angle = np.arctan(b_vec[1]/b_vec[2])

    if (angle_unit == "radians"):
        return angle
    elif (angle_unit == "degrees"):
        return angle * 180/np.pi
    else:
        raise KeyError("angle_unit must be radians or degrees")


def ridge_finder(
    image=None,
    t_range=['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
    xrange=[-15.1, 15],
    yrange=[-15.1, 15],
    dr=0.5,
    dipole_tilt_angle=None,
    sigma=2.2,
    mode="nearest",
    alpha=1.,
    vmin=None,
    vmax=None,
    cmap="viridis",
    draw_patch=False,
    draw_ridge=False,
    save_fig=True,
    fig_name="new",
    fig_format="pdf",
    c_label="none",
    c_unit="none",
    ):
    r"""
    Finds ridges in an image and plot the points with maximum ridge value on the given image.

    Parameters
    ----------
    image : ndarray
            The image to find ridges in. Default is None.
    xrange : list of floats, optional
            The range of x-values for image. Default is [-15.1, 15].
    yrange : list of floats, optional
            The range of y-values for image. Default is [-15.1, 15].
    dr : float, optional
            The step size for the grid. Default is 0.5.
    sigma : float
            The size of the filter. Default is 2.2.
    mode : str
            The mode of the filter. Can be 'nearest', 'reflect', 'constant', 'mirror', 'wrap' or
            'linear'. Default is 'nearest'.
    alpha : float
            The alpha value for the filter. Default is 0.5.
    vmin : float, optional
            The minimum value of the colorbar. Default is None.
    vmax : float, optional
            The maximum value of the colorbar. Default is None.
    cmap : str, optional
            The colormap to use. Default is 'viridis'.
    draw_patch : bool, optional
            Whether to draw the circular patch. Default is False.
    draw_ridge : bool, optional
            Whether to draw the ridge line. Default is False.
    save_fig : bool, optional
            Whether to save the figure. Default is True.
    fig_name : str, optional
            The name of the figure. Default is "new".
    fig_format : str, optional
            The format of the figure. Default is "pdf".
    c_label : str, optional
            The label for the colorbar. Default is "none".
    c_unit : str, optional
            The units for the colorbar label. Default is "none".

    Raises
    ------
    ValueError: If the image is not a numpy array.

    Returns
    -------
    ridge_points : ndarray
    """
    if image is None:
        raise ValueError("No image given")

    # NOTE: This is a hack to ensure that the output of shear angle, reconnection energy etc. agrees
    # with what has been reported in literature. Plot for shear angle seems to agree reasonably well
    # (for "trange = ['2016-12-07 05:11:00', '2016-12-07 05:21:00']") with the one reported by
    # FuselierJGR2019 (doi:10.1029/2019JA027143, see fig. 4).
    image_rotated = np.transpose(image)

    if cmap is None:
        cmap = "viridis"
    if(vmin is not None and vmax is not None):
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = plt.Normalize()

    #cmap = plt.cm.jet

    kwargs = {'sigmas': [3], 'black_ridges': False, 'mode': mode, 'alpha': 1}

    # Smoothen the image
    image_smooth = sp.ndimage.filters.gaussian_filter(image_rotated, sigma=[5, 5], mode=mode)
    result = meijering(image_smooth, **kwargs)

    x_len = image_rotated.shape[0]
    y_len = image_rotated.shape[1]

    y_val = np.full(y_len, np.nan)
    im_max_val = np.full(y_len, np.nan)
    for i in range(y_len):
        y_val[i] = np.argmax(result[:, i]) * dr + yrange[0]
        im_max_val[i] = np.argmax(image_rotated[:, i]) * dr + yrange[0]

    # plt.close('all')
    fig, axs1 = plt.subplots(1, 1, figsize=(8, 6))

    im1 = axs1.imshow(image_smooth, extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                      origin='lower', cmap=cmap, norm=norm)
    divider1 = make_axes_locatable(axs1)

    # Take rolling average of the y_val array
    y_val_avg = np.full(len(y_val), np.nan)
    for i in range(len(y_val)):
        y_val_avg[i] = np.nanmean(y_val[max(0, i-5):min(len(y_val), i+5)])
    
    if draw_ridge:
        axs1.plot(np.linspace(xrange[0], xrange[1], x_len), y_val_avg, 'k-', alpha=0.9)
        axs1.plot(np.linspace(xrange[0], xrange[1], x_len), im_max_val, 'k*', ms=1, alpha=0.5)

    # Plot a horizontal line at x=0 and a vertical line at y=0
    axs1.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    axs1.axvline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)

    if(draw_patch):
        patch = patches.Circle((0, 0), radius=(xrange[1] - xrange[0])/2., transform=axs1.transData,
                            fc='none', ec='k', lw=0.1)
        axs1.add_patch(patch)
        im1.set_clip_path(patch)

    axs1.set_xlabel(r'Y [GSM, $R_\oplus$]', fontsize=18)
    axs1.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=18)

    # Define the location of the colorbar, it's size relative to main figure and the padding
    # between the colorbar and the figure, the orientation the colorbar
    cax1 = divider1.append_axes("top", size="5%", pad=0.01)
    cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05,
                         pad=0.01)
    cbar1.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                         labelbottom=False, pad=0.01)
    cbar1.ax.xaxis.set_label_position('top')
    cbar1.ax.set_xlabel(f'{c_label} ({c_unit})', fontsize=18)


    # Write the timme range on the plot
    axs1.text(1.0, 0.5, f'Time range: {trange[0]} - {trange[1]}', horizontalalignment='left',
              verticalalignment='center', transform=axs1.transAxes, rotation=270, color='r')
    axs1.text(0.01, 0.99, f'Clock Angle: {np.round(imf_clok_angle, 2)}$^\circ$',
    horizontalalignment='left', verticalalignment='top', transform=axs1.transAxes, rotation=0,
    color='r')
    axs1.text(0.99, 0.99, f'Dipole tilt: {np.round(dipole_tilt_angle * 180/np.pi, 2)} $^\circ$',
    horizontalalignment='right', verticalalignment='top', transform=axs1.transAxes, rotation=0, color='r')
    # fig.show()

    if save_fig:
        try:
            fig_name = f'../figures/ridge_plot_vir_{fig_name}_{dr}dr_{m_p}mp_{t_range[0][:10]}_{t_range[0][-8:]}_{t_range[1][:10]}_{t_range[1][-8:]}.{fig_format}'
            plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=fig_format, dpi=300)
            print(f'Figure saved as {fig_name}')
        except  Exception as e:
            print(e)
            print(f'Figure not saved, folder does not exist. Create folder ../figures')
            #pass
        plt.close()
    return y_val


def draping_field_plot(x_coord=None, y_coord=None, by=None, bz=None, scale=None, save_fig=True,
                       fig_name="draping_field", figure_format="pdf", tick_fontsize=16,
                       label_fontsize=18):
    r"""
    Plots the draping field.

    Parameters
    ----------
    x_coord : ndarray
        The x-coordinates of the points.
    y_coord : ndarray
        The y-coordinates of the points.
    by : ndarray
        The y-component of the draping field.
    bz : ndarray
        The z-component of the draping field.
    scale : float, optional
        Number of data units per arrow length unit. Default is None.
    save_fig : bool, optional
        Whether to save the figure. Default is True.
    fig_name : str, optional
        The name of the figure. Default is 'draping_field'.
    figure_format : str, optional
        The format of the figure. Default is 'pdf'.
    Raises
    ------
    ValueError: If the x_coord, y_coord, by or bz are not numpy arrays.

    Returns
    -------
    fig : matplotlib figure
    """
    if x_coord is None or y_coord is None or by is None or bz is None:
        raise ValueError("No coordinates or field components given")

    fig, axs1 = plt.subplots(1, 1, figsize=(8, 6))
    im1 = axs1.quiver(x_coord, y_coord, by, bz, scale=scale, scale_units='inches', angles='uv', width=0.002)
    axs1.set_xlabel(r'Y [GSM, $R_\oplus$]', fontsize=label_fontsize)
    axs1.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=label_fontsize)
    patch = patches.Circle((0, 0), radius=15, transform=axs1.transData, fc='none', ec='none', lw=0.1)
    axs1.add_patch(patch)
    axs1.set_clip_path(patch)
    im1.set_clip_path(patch)

    axs1.set_xlim([-15, 15])
    axs1.set_ylim([-15, 15])
    axs1.set_aspect('equal')
    if fig_name == "magnetosheath":
        axs1.set_title(r'Magnetosheath Field (nT)', fontsize=label_fontsize)
    elif fig_name == "magnetosphere":
        axs1.set_title(r'Magnetosphere Field (nT)', fontsize=label_fontsize)
    else:
        axs1.set_title(r'Draping Field (nT)', fontsize=label_fontsize)
    
    axs1.tick_params(axis='both', which='major', labelsize=tick_fontsize)

    if save_fig:
        fig_name = f'../figures/draping_field_plot_{fig_name}.{figure_format}'
        plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=figure_format, dpi=300)
        print(f'Figure saved as {fig_name}')
    return fig


def model_run(*args):
    """
    Returns the value of the magnetic field at a given point in the model grid using three different
    models
    """
    j = args[0][0]
    k = args[0][1]
    y_max = args[0][2]
    z_max = args[0][3]
    model_type = args[0][4]

    print(j, k)
    y0 = int(j * dr) - y_max
    z0 = int(k * dr) - z_max
    rp = np.sqrt(y0**2 + z0**2)  # Projection of r into yz-plane

    for index in range(0, 100):

        theta = index * d_theta
        r = ro * (2/(1 + np.cos(theta))) ** alpha
        zp = r * np.sin(theta)  # not really in z direction, but a distance in yz plane
        x0 = r * np.cos(theta)

        if x0 == 0:
            signx = 1.0
        else:
            signx = np.sign(x0)

        if y0 == 0:
            signy = 1.0
        else:
            signy = np.sign(y0)

        if z0 == 0:
            signz = 1.0
        else:
            signz = np.sign(z0)

        if (rp <= zp):
            # print(index, rp, zp)
            # print(f'Value of theta = {theta}')

            y_coord = y0
            z_coord = z0
            x_shu = (r -m_p) * np.cos(theta)
            phi = np.arctan2(z0, y0)

            if (abs(y0) == 0 or abs(z0) == 0):
                if(abs(y0) == 0):
                    y_shu = 0
                    z_shu = (r -m_p) * np.sin(theta)
                elif (abs(z0) == 0):
                    z_shu = 0
                    y_shu = (r -m_p) * np.sin(theta)
            else:
                z_shu = np.sqrt((rp - 1.0)**2/(1 + np.tan(phi)**(-2)))
                y_shu = z_shu/np.tan(phi)

            rho_sh = rho * (1.509 * np.exp(x_shu/rmp) + .1285)
            #print(rmp)
            n_sh = rho_sh/m_proton

            y_shu = abs(y_shu)*signy
            z_shu = abs(z_shu)*signz

            #print(f'z_shu = {z_shu}, y_shu = {y_shu}')

            # Cooling JGR 2001 Model, equation 9 to 12
            # the distance from the focus to the magnetopause surface
            ll = 3 * rmp/2 - x0
            b_msx = - A * (- b_imf_x * (1 - rmp / (2 * ll)) + b_imf_y * (y0 / ll)
                           + b_imf_z * (z0 / ll))
            b_msy = A * (- b_imf_x * (y0 / (2 * ll)) + b_imf_y * (2 - y0**2/(
                           ll * rmp)) - b_imf_z * (y0 * z0 / (ll * rmp)))
            b_msz = A * (- b_imf_x * (z0 / (2 * ll)) - b_imf_y * (y0 * z0 / (ll
                         * rmp)) + b_imf_z * (2 - z0**2 / (ll * rmp)))
            try:
                if model_type == 't96':
                    bx_ext, by_ext, bz_ext = gp.t96.t96(param, ps, x_shu, y_shu, z_shu)
                elif model_type == 't01':
                    bx_ext, by_ext, bz_ext = gp.t01.t01(param, ps, x_shu, y_shu, z_shu)
            except exception as e:
                print(e)
                print(f'Skipped for {x_shu, y_shu, z_shu}')
                bx_ext = np.nan
                by_ext = np.nan
                bz_ext = np.nan
                print(f'{bx_ext, by_ext, bz_ext}, {j}, {k}')
                #pass

            bx_igrf, by_igrf, bz_igrf = gp.igrf_gsm(x_shu, y_shu, z_shu)

            #print(j, k, bx_ext, bx_igrf)
            # If both bx_ext and bx_igrf are nan, then print the value of j and k
            if np.isnan(bx_ext) and np.isnan(bx_igrf):
                print(f'{j}, {k}')
            bx = bx_ext + bx_igrf
            by = by_ext + by_igrf
            bz = bz_ext + bz_igrf
            #print(j, k, bx)
            #if (np.sqrt(y_shu**2 + z_shu**2) > 31):
            #    shear = np.nan
            #    rx_en = np.nan
            #    va_cs = np.nan
            #    bisec_msp = np.nan
            #    bisec_msh = np.nan
            #else:
            shear = get_shear([bx, by, bz], [b_msx, b_msy, b_msz], angle_unit="degrees")

            rx_en = get_rxben([bx, by, bz], [b_msx, b_msy, b_msz])
            va_cs = get_vcs([bx, by, bz], [b_msx, b_msy, b_msz], n_sh, 0.1)
            bisec_msp, bisec_msh = get_bis([bx, by, bz], [b_msx, b_msy, b_msz])
            break

    return j, k, bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh


#def rx_model_batch(
#    probe=None,
#    omni_level='hro',
#    maximum_shear=True,
#    mms_probe='mms1',
#    movie=None,
#    trange=None,
#    model_type='t96',
#    m_p=0.5,
#    dr=0.25
#    y_min= -20,
#    y_max= 20,
#    z_min= -20,
#    z_max= 20,
#    draw_patch=False,
#    save_data=False,
#    plot_type="shear"
#):
code_run = '1'
if(code_run):
#for i in range(1):
    r"""
    RX model from the IDL code.

    Parameters
    ----------
    probe : TYPE, optional
        DESCRIPTION. The default is probe.
    maximum_shear : TYPE, optional
        DESCRIPTION. The default is maximum_shear.
    movie : TYPE, optional.
        The default is movie.
    mmsprobe : TYPE, optional
        DESCRIPTION. The default is mmsprobe.
    times : TYPE, optional
        DESCRIPTION. The default is times.
    trange : TYPE, optional
        Range of times to use. The default is trange.
    model_type : TYPE, optional
        Type of Tsyganenko model. The default is 't96'.
    m_p : float, optional
        Thickness of the magnetopause. The default is 0.5.
    dr : float, optional
        Resolution of the model. The default is 0.25.
    y_min : float, optional
        Minimum y value. The default is -20.
    y_max : float, optional
        Maximum y value. The default is 20.
    z_min : float, optional
        Minimum z value. The default is -20.
    z_max : float, optional
        Maximum z value. The default is 20.
    draw_patch : bool, optional
        Draw the patch of a given radius around the figure. The default is False. If it is set to
        true then the radius of the patch is automatically computed as the maximum value of the
        x-axis.
    save_data : bool, optional
        If True, the data will be saved in an HDF file. The default is False.
    plot_type : String or array of strings, optional
        Type of the plot one wants to get once the code has finished running. The default is
        "shear". Other options are "rx_en", "va_cs", "bisec_msp", "bisec_msh" and "all" for all of
        them.

    Raises
    ------
    KeyError:
        If "plot_type" is not in the list of plot types.
    Returns
    -------
    None.
    """
    probe = None
    omni_level = 'hro'
    maximum_shear = True
    mms_probe = None
    movie = None
    #trange = ['2016-12-24 15:08:00', '2016-12-24 15:12:00']
    #trange = ['2016-12-07 05:11:00', '2016-12-07 05:21:00']
    #trange = ['2016-12-29 03:53:00', '2016-12-29 04:03:00']
    trange = ['2015-09-08 11:05:00', '2015-09-08 11:15:00']
    model_type = 't96'
    m_p = 0.5  # Magnetopause thichkness
    dr = 0.5  # Resolution of model run in R_E units 
    min_max_val = 20
    y_min = -min_max_val
    y_max = min_max_val
    z_min = -min_max_val
    z_max = min_max_val
    save_data = False
    plot_type = "ttt"
    draw_patch = True
    # Get time range as a datetime object
    trange_unix = [parser.parse(xx) for xx in trange]
    # For MMS add 5 minutes buffer on either side of the time range and give out 'trange_mms' in a
    # format that can be read by PySpedas
    trange_mms_dt = [trange_unix[0] - datetime.timedelta(minutes=5), 
                     trange_unix[1] + datetime.timedelta(minutes=5)]
    trange_mms = [xx.strftime('%Y-%m-%d %H-%M-%S') for xx in trange_mms_dt]

    # Download the OMNI data (default level of 'hro_1min') for the specified timerange.
    #if ('omni_vars' in locals()):
    #    pass
    #else:
    omni_vars = spd.omni.data(trange=trange, level=omni_level, time_clip=True)

    omni_time = ptt.get_data('BX_GSE')[0]
    # omni_time_unix = np.vectorize(datetime.utcfromtimestamp)(omni_time[:]) # converting omni_time
    # from unixtime to utc datetime object array in python

    omni_bx_gse = ptt.get_data('BX_GSE')[1]
    omni_by_gse = ptt.get_data('BY_GSE')[1]
    omni_bz_gse = ptt.get_data('BZ_GSE')[1]
    omni_by_gsm = ptt.get_data('BY_GSM')[1]
    omni_bz_gsm = ptt.get_data('BZ_GSM')[1]
    omni_np = ptt.get_data('proton_density')[1]
    omni_vx = ptt.get_data('Vx')[1]
    omni_vy = ptt.get_data('Vy')[1]
    omni_vz = ptt.get_data('Vz')[1]
    omni_sym_h = ptt.get_data('SYM_H')[1]

    if (probe is None):
        pass
    elif ~hasattr(probe, '__len__'):
        themis_stt_vars = spd.themis.state(probe=probe, trange=trange)
        themis_time = ptt.get_data(f'th{probe}_pos_gse')[0]
        themis_sc_pos = ptt.get_data(f'th{probe}_pos_gse')[1]
    elif hasattr(probe, '__len__'):
        themis_vars = []
        sc_pos = []
        for xx in probe:
            themis_vars.append(spd.themis.state(probe=xx, trange=trange))
            sc_pos.append(ptt.get_data(f'th{xx}_pos_gse')[1])

    if (mms_probe is not None):
        #if('mms_mec_vars' in locals()):
        #    pass
        #else:
        mms_mec_vars = spd.mms.mec(trange=trange_mms, data_rate='srvy', probe='1')
        mms_time = ptt.get_data('mms1_mec_r_gsm')[0]
        mms_sc_pos = ptt.get_data('mms1_mec_r_gsm')[1:3]
    else:
        pass

    if (movie is None):
        time_imf = np.nanmedian(omni_time)
        #print(time_imf, type(time_imf))
        b_imf_x = np.nanmedian(omni_bx_gse)
        b_imf_y = np.nanmedian(omni_by_gsm)
        b_imf_z = np.nanmedian(omni_bz_gsm)

        if (b_imf_z > 15 or b_imf_z < -18):
            warnings.warn(
            f"The given parameters produced the z-component of IMF field (b_imf_z) {b_imf_z} nT,"
            f"which is out of range in which model is valid (-18 nT < b_imf_z < 15 nT)"
            )

        time_imf_hrf = datetime.datetime.utcfromtimestamp(time_imf)
        np_imf = np.nanmedian(omni_np)
        vx_imf = np.nanmedian(omni_vx)
        vy_imf = np.nanmedian(omni_vy)
        vz_imf = np.nanmedian(omni_vz)
        sym_h_imf = np.nanmedian(omni_sym_h)
        v_imf = [vx_imf, vy_imf, vz_imf]
        b_imf = [b_imf_x, b_imf_y, b_imf_z]
        imf_clok_angle = np.arctan2(b_imf[1], b_imf[2]) * 180 / np.pi
        print("IMF parameters found:")
        print(tabulate(
            [["Time of observation (UTC)", time_imf_hrf],
             ["IMF Magnetic field [GSM] (nT)", b_imf],
             ["IMF Proton density (1/cm^-3)", np_imf],
             ["IMF Plasma velocity (km/sec)", v_imf],
             ["IMF clock angle (degrees)", imf_clok_angle],
             ["IMF Sym H", sym_h_imf]],
            headers=["Parameter", "Value"], tablefmt="fancy_grid", floatfmt=".2f",
            numalign="center"))

        # Check if the values are finite, if not then assign a default value to each of them
        if ~(np.isfinite(np_imf)):
            np_imf = 5
        if ~(np.isfinite(vx_imf)):
            vx_imf = -500
        if ~(np.isfinite(vy_imf)):
            vy_imf = 0
        if ~(np.isfinite(vz_imf)):
            vz_imf = 0
        if ~(np.isfinite(sym_h_imf)):
            sym_h_imf = -1

        m_proton = 1.672e-27  # Mass of proton in SI unit

        rho = np_imf * m_proton * 1.15

        #  Solar wind ram pressure in nPa, including roughly 4% Helium++ contribution
        p_dyn = 1.6726e-6 * 1.15 * np_imf * (vx_imf**2 + vy_imf**2 + vz_imf**2)

        if (p_dyn > 8.5 or p_dyn < 0.5):
            warnings.warn(
                f"The given parameters produced a dynamic pressure of {p_dyn} nPa which is out of"
                f" range in which model is valid (0.5 nPa < p_dyn < 8.5 nPa)",
            )
        param = [p_dyn, sym_h_imf, b_imf_y, b_imf_z, 0, 0, 0, 0, 0, 0]

        # Compute the dipole tilt angle
        ps = gp.recalc(time_imf)
        # ps = -0.25666021186831828

        n_arr_y = int((y_max - y_min) / dr) + 1
        n_arr_z = int((z_max - z_min) / dr) + 1

        bx = np.full((n_arr_y, n_arr_z), np.nan)
        by = np.full((n_arr_y, n_arr_z), np.nan)
        bz = np.full((n_arr_y, n_arr_z), np.nan)

        bx_ext = np.full((n_arr_y, n_arr_z), np.nan)
        by_ext = np.full((n_arr_y, n_arr_z), np.nan)
        bz_ext = np.full((n_arr_y, n_arr_z), np.nan)

        bx_igrf = np.full((n_arr_y, n_arr_z), np.nan)
        by_igrf = np.full((n_arr_y, n_arr_z), np.nan)
        bz_igrf = np.full((n_arr_y, n_arr_z), np.nan)

        b_msx = np.full((n_arr_y, n_arr_z), np.nan)
        b_msy = np.full((n_arr_y, n_arr_z), np.nan)
        b_msz = np.full((n_arr_y, n_arr_z), np.nan)

        x_shu = np.full((n_arr_y, n_arr_z), np.nan)
        y_shu = np.full((n_arr_y, n_arr_z), np.nan)
        z_shu = np.full((n_arr_y, n_arr_z), np.nan)

        rho_sh = np.full((n_arr_y, n_arr_z), np.nan)

        # rp = np.full((n_arr, n_arr), np.nan)

        # r = np.full((n_arr, n_arr), np.nan)
        # zp = np.full((n_arr, n_arr), np.nan)
        # x0 = np.full((n_arr, n_arr), np.nan)

        shear = np.full((n_arr_y, n_arr_z), np.nan)
        rx_en = np.full((n_arr_y, n_arr_z), np.nan)
        gd1 = np.full((n_arr_y, n_arr_z), np.nan)
        gd2 = np.full((n_arr_y, n_arr_z), np.nan)
        va_cs = np.full((n_arr_y, n_arr_z), np.nan)
        bisec_msp = np.full((n_arr_y, n_arr_z), np.nan)
        bisec_msh = np.full((n_arr_y, n_arr_z), np.nan)
        y_coord = np.full((n_arr_y, n_arr_z), np.nan)
        z_coord = np.full((n_arr_y, n_arr_z), np.nan)
        b_sh_ca = np.full((n_arr_y, n_arr_z), np.nan)
        b_sh_mag = np.full((n_arr_y, n_arr_z), np.nan)
        n_sh = np.full((n_arr_y, n_arr_z), np.nan)

        d_theta = np.pi/100

        # Shue et al.,1998, equation 10
        # if (b_imf_z >=0):
        #     ro = (11.44 + 0.013 * b_imf_z) * (p_dyn)**(-1.0/6.6)
        #     #ro = (10.22 + 1.29 * np.tanh(0.184 * (b_imf_z + 8.14))) * (p_dyn)**(-1.0/6.6)
        # else:
        #    ro = (11.44 + 0.14 * b_imf_z) * (p_dyn)**(-1.0/6.6)
        ro = (10.22 + 1.29 * np.tanh(0.184 * (b_imf_z + 8.14))) * (p_dyn)**(-1.0/6.6)

        # Shue et al.,1998, equation 11
        # alpha = (0.58 - 0.010 * b_imf_z) * (1 + 0.010 * p_dyn)
        alpha = (0.58 - 0.007 * b_imf_z) * (1 + 0.024 * np.log(p_dyn))
        rmp = ro * (2/(1 + np.cos(0.0))) ** alpha  # Stand off position of the magnetopause

        A = 2
        len_y = int((y_max - y_min)/dr) + 1
        len_z = int((z_max - z_min)/dr) + 1

        p = mp.Pool(processes=1)

        input = ((j, k, y_max, z_max, 't96') for j in range(len_y) for k in range(len_z))

        print("Running the model\n")
        res = p.map(model_run, input)
        print("Model run complete")

        p.close()
        p.join()

        for r in res:
            j = r[0]
            k = r[1]
            bx[j, k] = r[2]
            by[j, k] = r[3]
            bz[j, k] = r[4]

            shear[j, k] = r[5]
            rx_en[j, k] = r[6]
            va_cs[j, k] = r[7]
            bisec_msp[j, k] = r[8]
            bisec_msh[j, k] = r[9]

    if save_data:
        try:
            fn = f'../data/all_data_rx_model_{dr}re_{m_p}mp_{model_type}_{today_date}.h5'
            data_file = hf.File(fn, 'w')

            data_file.create_dataset('bx', data=bx)
            data_file.create_dataset('by', data=by)
            data_file.create_dataset('bz', data=bz)

            data_file.create_dataset('shear', data=shear)
            data_file.create_dataset('rx_en', data=rx_en)
            data_file.create_dataset('va_cs', data=va_cs)
            data_file.create_dataset('bisec_msp', data=bisec_msp)
            data_file.create_dataset('bisec_msh', data=bisec_msh)

            data_file.close()
            print(f'Date saved to file {fn}')
        except Exception as e:
            print(e)
            print(f'Data not saved to file {fn}. Please make sure that file name is correctly assigned and that the directory exists and you have write permissions')

    # Check if 'plot_type' has length attribute. If it has length attribute then plot the ridge plot
    # for each of the plot type in the list. If it does not have length attribute then plot the
    # ridge plot for the specified plot type.
    plot_type = 'shear'
    types_of_plot = ['shear', 'rx_en', 'va_cs', 'bisec_msp', 'bisec_msh', 'all']
    if isinstance(plot_type, list):
        for xx in plot_type:
            if xx not in types_of_plot:
                raise ValueError(
                    f'{xx} is not a valid plot type. Please choose from {types_of_plot}'
                    )
        if 'shear' in plot_type:
            print('Plotting shear')
            ridge_finder(image=shear, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
            sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='shear_vp', c_label='Shear',
            c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'rx_en' in plot_type:
            print('Plotting rx_en')
            ridge_finder(image=rx_en/np.nanmax(rx_en), t_range=trange, xrange=[y_min, y_max],
            yrange=[z_min, z_max], sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='rx_en_vp',
            c_label='Rx_en', c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'va_cs' in plot_type:
            print('Plotting va_cs')
            ridge_finder(image=va_cs, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
            sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='va_cs_vp', c_label='Va_cs',
            c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'bisec_msp' in plot_type:
            print('Plotting bisec_msp')
            ridge_finder(image=bisec_msp, t_range=trange, xrange=[y_min, y_max],
            yrange=[z_min, z_max], sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='bisec_msp_vp',
            c_label='Bisec_msp', c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'bisec_msh' in plot_type:
            print('Plotting bisec_msh')
            ridge_finder(image=bisec_msh, t_range=trange, xrange=[y_min, y_max],
            yrange=[z_min, z_max], sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='bisec_msh_vp',
            c_label='Bisec_msh', c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type == 'all':
        print('Plotting for all')
        ridge_finder(image=shear, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='shear_vp', c_label='Shear',
        c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=rx_en/np.nanmax(rx_en), t_range=trange, xrange=[y_min, y_max],
        yrange=[z_min, z_max], sigma=2.8, dr=dr, dipole_tilt_angle=ps,fig_name='rx-en_nPa_vp',
        c_label='Reconnection Energy', c_unit='nPa', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=va_cs, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=3., dr=dr, dipole_tilt_angle=ps,fig_name='va-cs_vp', c_label='Exhaust Velocity',
        c_unit='km/s', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=bisec_msp, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='bisec_msp_vp', c_label='Bisection Field',
        c_unit='nT', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=bisec_msh, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='bisec_msh_vp', c_label='Bisection Field',c_unit='nT', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='shear':
        print('Plotting shear')
        ridge_finder(image=shear, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=1, dr=dr, dipole_tilt_angle=ps, fig_name='shear_vp', c_label='Shear',
        c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='rx_en':
        print('Plotting rx_en')
        ridge_finder(image=rx_en/np.nanmax(rx_en), t_range=trange, xrange=[y_min, y_max],
        yrange=[z_min, z_max], sigma=2.8, dr=dr, dipole_tilt_angle=ps,fig_name='rx-en_nPa_vp',
        c_label='Reconnection Energy', c_unit='nPa', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='va_cs':
        print('Plotting va_cs')
        y_val =ridge_finder(image=va_cs, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=3., dr=dr, dipole_tilt_angle=ps,fig_name='va-cs_vp', c_label='Exhaust Velocity',
        c_unit='km/s', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='bisec_msp':
        print('Plotting bisec_msp')
        ridge_finder(image=bisec_msp, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='bisec_msp_vp', c_label='Bisection Field',c_unit='nT', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='bisec_msh':
        print('Plotting bisec_msh')
        ridge_finder(image=bisec_msh, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps,fig_name='bisec_msh_vp', c_label='Bisection Field',
        c_unit='nT', draw_patch=draw_patch, draw_ridge=True)
    else:
        raise KeyError('plot_type must be one or a list of: all, shear, rx-en, va-cs, bisec_msp, bisec_msh')

    #_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=b_msy, bz=b_msz, save_fig=True,
    #                      scale=40, fig_name="magnetosheath_vp")
    #_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=by, bz=bz, save_fig=True, scale=120,
    #                 fig_name="magnetosphere_vp")
print(f'Took {round(time.time() - start, 3)} seconds')