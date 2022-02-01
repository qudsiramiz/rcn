# functions to be used in rx_model_batch_parallel.py
import datetime
import multiprocessing as mp
import warnings

import geopack.geopack as gp
import h5py as hf
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pyspedas as spd
import pytplot as ptt
import scipy as sp
from dateutil import parser
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import frangi, hessian, meijering, sato
from tabulate import tabulate


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
    unit_vec_1 = b_vec_1 / np.linalg.norm(b_vec_1)
    unit_vec_2 = b_vec_2 / np.linalg.norm(b_vec_2)
    angle = np.arccos(np.dot(unit_vec_1, unit_vec_2))

    if (angle_unit == "radians"):
        return angle
    elif (angle_unit == "degrees"):
        return angle * 180 / np.pi
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

    #alpha = - 14.87 * np.pi / 180  # radians (From Hesse2013)
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

    # Reconnection energy (from Hesse2013)
    rx_en = rx_b_mag_1 ** 2 * rx_b_mag_2 ** 2

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

    # The bisector vector of thwo input vectors
    unit_vec_bisec = (unit_vec_1 + unit_vec_2) / np.linalg.norm(unit_vec_1 + unit_vec_2)

    # Cross product of the two input vectors with the bisector vector to get the reconnection
    # component of the magnetic field.
    rx_b_1 = np.cross(b_vec_1, unit_vec_bisec)
    rx_b_2 = np.cross(b_vec_2, unit_vec_bisec)

    # Magnitude of the reconnection component of the magnetic fields
    rx_mag_1 = np.linalg.norm(rx_b_1)
    rx_mag_2 = np.linalg.norm(rx_b_2)

    vcs = va_p1 * np.sqrt(rx_mag_1 * rx_mag_2 * (rx_mag_1 +
                          rx_mag_2) / (rx_mag_1 * n_2 + rx_mag_2 * n_1))

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
        The magnitude of bisected field line corresponding to the second input magnetic field
        vector.
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
    mms_probe_num='1',
    mms_sc_pos= [0, 0],
    dr=0.5,
    dipole_tilt_angle=None,
    imf_clock_angle=None,
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
    fig_format="png",
    c_label="none",
    c_unit="none",
    ):
    r"""
    Finds ridges in an image and plot the points with maximum ridge value on the given image.

    Parameters
    ----------
    image : ndarray
            The image to find ridges in. Default is None.
    t_range : list of str
            The time range to find the ridge in. Default is ['2016-12-24 15:08:00',
            '2016-12-24 15:12:00'].
    xrange : list of floats, optional
            The range of x-values for image. Default is [-15.1, 15].
    yrange : list of floats, optional
            The range of y-values for image. Default is [-15.1, 15].
    mms_probe_num : str, optional
            The probe number of the MMS spacecraft. Default is 1.
    mms_sc_pos : list of floats, optional
            The position of the spacecraft in the image. Default is [0, 0].
    dr : float, optional
            The step size for the grid. Default is 0.5.
    dipole_tilt_angle : float, optional
            The dipole tilt angle. Default is None.
    imf_clock_angle : float, optional
            The IMF clock angle. Default is None.
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

    # Draw the spacecraft position
    axs1.plot(mms_sc_pos[0], mms_sc_pos[1], 'k', marker=mms_probe_num, ms=10, alpha=1)

    # Write the timme range on the plot
    axs1.text(1.0, 0.5, f'Time range: {t_range[0]} - {t_range[1]}', horizontalalignment='left',
              verticalalignment='center', transform=axs1.transAxes, rotation=270, color='r')
    axs1.text(0.01, 0.99, f'Clock Angle: {np.round(imf_clock_angle, 2)}$^\circ$',
    horizontalalignment='left', verticalalignment='top', transform=axs1.transAxes, rotation=0,
    color='r')
    axs1.text(0.99, 0.99, f'Dipole tilt: {np.round(dipole_tilt_angle * 180/np.pi, 2)} $^\circ$',
              horizontalalignment='right', verticalalignment='top', transform=axs1.transAxes,
              rotation=0, color='r')

    if save_fig:
        try:
            fig_time_range = f"{parser.parse(t_range[0]).strftime('%Y-%m-%d_%H-%M-%S')}_{parser.parse(t_range[1]).strftime('%Y-%m-%d_%H-%M-%S')}"
            fig_name = f'../figures/{fig_name}/ridge_plot_{fig_name}_{fig_time_range}.{fig_format}'
            plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=fig_format, dpi=300)
            print(f'Figure saved as {fig_name}')
        except  Exception as e:
            print(e)
            print(f'Figure not saved, folder does not exist. Create folder ../figures')
            #pass
        plt.close()
    return y_val


def ridge_finder_multiple(
    image=[None, None, None, None],
    t_range=['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
    dt=5,
    xrange=[-15.1, 15],
    yrange=[-15.1, 15],
    mms_probe_num='1',
    mms_sc_pos= [0, 0],
    dr=0.5,
    dipole_tilt_angle=None,
    imf_clock_angle=None,
    sigma=[2.2, 2.2, 2.2, 2.2],
    mode="nearest",
    alpha=1.,
    vmin=[None, None, None, None],
    vmax=[None, None, None, None],
    cmap_list=["viridis", "viridis", "viridis", "viridis"],
    draw_patch=[True, True, True, True],
    draw_ridge=[True, True, True, True],
    save_fig=True,
    fig_name="new",
    fig_format="png",
    c_label=[None, None, None, None],
    c_unit=[None, None, None, None],
    wspace=0.1,
    hspace=0.1,
    fig_size=(10, 10),
    box_style=None,
    title_y_pos=0.95,
    interpolation='nearest',
    ):
    r"""
    Finds ridges in an image and plot the points with maximum ridge value on the given image.

    Parameters
    ----------
    image : list of numpy arrays
        List of images to be plotted.
    t_range : list of str
            The time range to find the ridge in. Default is ['2016-12-24 15:08:00',
            '2016-12-24 15:12:00'].
    dt : float, optional
        The time differential, in minutes, for observation if 't_range' has only one element.
        Default is 5 minutes.
    xrange : list of floats, optional
            The range of x-values for image. Default is [-15.1, 15].
    yrange : list of floats, optional
            The range of y-values for image. Default is [-15.1, 15].
    mms_probe_num : str, optional
            The probe number of the MMS spacecraft. Default is 1.
    mms_sc_pos : list of floats, optional
            The position of the spacecraft in the image. Default is [0, 0].
    dr : float, optional
            The step size for the grid. Default is 0.5.
    dipole_tilt_angle : float, optional
            The dipole tilt angle. Default is None.
    imf_clock_angle : float, optional
            The IMF clock angle. Default is None.
    sigma : list of floats, optional
            List of sigmas to be used for the ridge plot. Default is [2.2, 2.2, 2.2, 2.2].
    mode : str
            The mode of the filter. Can be 'nearest', 'reflect', 'constant', 'mirror', 'wrap' or
            'linear'. Default is 'nearest'.
    alpha : float
            The alpha value for the filter. Default is 0.5.
    vmin : list of floats, optional
            List of vmin values for the ridge plot. Default is [None, None, None, None].
    vmax : list of floats, optional
            List of vmax values for the ridge plot. Default is [None, None, None, None].
    cmap_list : list of str, optional
            List of colormaps to be used for the ridge plot. Default is ['viridis', 'viridis',
            'viridis', 'viridis'].
    draw_patch : list of bool, optional
            Whether to draw the circular patch. Default is [True, True, True, True].
    draw_ridge : list of bool, optional
            Whether to draw the ridge line. Default is [True, True, True, True].
    save_fig : bool, optional
            Whether to save the figure. Default is True.
    fig_name : str, optional
            The name of the figure. Default is "new".
    fig_format : str, optional
            The format of the figure. Default is "pdf".
    c_label : list of str, optional
            List of colorbar labels. Default is [None, None, None, None].
    c_unit : list of str, optional
            List of colorbar units. Default is [None, None, None, None].
    wspace : float, optional
            The width space between subplots. Default is 0.1.
    hspace : float, optional
            The height space between subplots. Default is 0.1.
    fig_size : tuple of floats, optional
            The size of the figure. Default is (10, 10).
    box_style : dict, optional
            The style of the box. Default is None.
    title_y_pos : float, optional
            The y-position of the title. Default is 0.95.
    interpolation : str, optional
            The interpolation method for imshow. Default is 'nearest'.
            Options are 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning',
            'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell'

    Raises
    ------
    ValueError: If the image is not a numpy array.

    Returns
    -------
    ridge_points : ndarray
    """
    if image is None:
        raise ValueError("No image given")

    if len(t_range) == 1:
        t_range_date = datetime.datetime.strptime(t_range[0], '%Y-%m-%d %H:%M:%S')
        t_range_date_min = t_range_date - datetime.timedelta(minutes=dt)
        t_range_date_max = t_range_date + datetime.timedelta(minutes=dt)
        t_range = [t_range_date_min.strftime('%Y-%m-%d %H:%M:%S'),
                  t_range_date_max.strftime('%Y-%m-%d %H:%M:%S')]

    fig = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor='w', edgecolor='gray')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=wspace, hspace=hspace)
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1])

    # Set the font size for the axes
    label_size = 20  # fontsize for x and y labels
    t_label_size = 18  # fontsize for tick label
    c_label_size = 18  # fontsize for colorbar label
    ct_tick_size = 14  # fontsize for colorbar tick labels
    l_label_size = 14  # fontsize for legend label

    #box_style = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    box_style = box_style
    y_vals = []
    for i in range(len(image)):
        image_rotated = np.transpose(image[i])

        if cmap_list is None:
            cmap_list = ["viridis", "viridis", "viridis", "viridis"]
        else:
            cmap_list = cmap_list
        if(vmin is not None and vmax is not None):
            norm = plt.Normalize(vmin=vmin[i], vmax=vmax[i])
        else:
            norm = plt.Normalize()

        kwargs = {'sigmas': [sigma[i]], 'black_ridges': False, 'mode': mode, 'alpha': 1}

        # Smoothen the image
        image_smooth = sp.ndimage.filters.gaussian_filter(image_rotated, sigma=[5, 5], mode=mode)
        result = meijering(image_smooth, **kwargs)  #frangi, hessian, meijering, sato

        x_len = image_rotated.shape[0]
        y_len = image_rotated.shape[1]

        y_val = np.full(y_len, np.nan)
        y_vals.append(y_val)
        im_max_val = np.full(y_len, np.nan)
        for xx in range(y_len):
            y_val[xx] = np.argmax(result[:, xx]) * dr + yrange[0]
            im_max_val[xx] = np.argmax(image_rotated[:, xx]) * dr + yrange[0]

        # plt.close('all')
        # TODO: Find a better way to do this
        if i==0:
            j = 0
            k = 0
        elif i==1:
            j = 0
            k = 1
        elif i==2:
            j = 1
            k = 0
        elif i==3:
            j = 1
            k = 1

        axs1 = plt.subplot(gs[j, k])
        im1 = axs1.imshow(image_smooth, extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                          origin='lower', cmap=cmap_list[i], norm=norm, interpolation=interpolation,
                          alpha=1)
        divider1 = make_axes_locatable(axs1)

        # Take rolling average of the y_val array
        y_val_avg = np.full(len(y_val), np.nan)
        for xx in range(len(y_val)):
            y_val_avg[xx] = np.nanmean(y_val[max(0, xx-25):min(len(y_val), xx+25)])
    
        if draw_ridge:
            axs1.plot(np.linspace(xrange[0], xrange[1], x_len), y_val_avg, 'k-', alpha=0.9)
            axs1.plot(np.linspace(xrange[0], xrange[1], x_len), im_max_val, 'k*', ms=1, alpha=0.5)

        # Plot a horizontal line at x=0 and a vertical line at y=0
        axs1.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
        axs1.axvline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)

        if(draw_patch):
            patch = patches.Circle((0, 0), radius=(xrange[1] - xrange[0])/2.,
                                    transform=axs1.transData, fc='none', ec='k', lw=0.1)
            axs1.add_patch(patch)
            im1.set_clip_path(patch)

        if i==0 or i==2:
            axs1.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=label_size)
        if i==2 or i==3:
            axs1.set_xlabel(r'Y [GSM, $R_\oplus$]', fontsize=label_size)
        if i==1 or i==3:
            axs1.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=label_size)
            axs1.yaxis.set_label_position("right")

        # Define the location of the colorbar, it's size relative to main figure and the padding
        # between the colorbar and the figure, the orientation the colorbar
        cax1 = divider1.append_axes("top", size="5%", pad=0.01)
        cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05,
                             pad=0.01)
        cbar1.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                             labelbottom=False, pad=0.01, labelsize=ct_tick_size)
        cbar1.ax.xaxis.set_label_position('top')
        cbar1.ax.set_xlabel(f'{c_label[i]} ({c_unit[i]})', fontsize=c_label_size)

        # Draw the spacecraft position
        axs1.plot(mms_sc_pos[0], mms_sc_pos[1], 'white', marker=mms_probe_num, ms=10, alpha=1)

        # Set tick label parameters
        if i==0 or i==2:
            axs1.tick_params(axis='both', direction='in', which='major', left=True, right=True,
                             top=True, bottom=True, labelleft=True, labelright=False,
                             labeltop=False, labelbottom=True, labelsize=t_label_size)
        else:
            axs1.tick_params(axis='both', direction='in', which='major', left=True, right=True,
                             top=True, bottom=True, labelleft=False, labelright=True,
                             labeltop=False, labelbottom=True, labelsize=t_label_size)
        # Write the timme range on the plot
        if i==2:
            axs1.text(-0.17, -0.1, f'Clock Angle: {np.round(imf_clock_angle, 2)}$^\circ$',
                  horizontalalignment='left', verticalalignment='top', transform=axs1.transAxes,
                  rotation=0, color='white', fontsize=l_label_size, bbox=box_style)
        elif i==3:
            axs1.text(1.17, -0.1,
                      f'Dipole tilt: {np.round(dipole_tilt_angle * 180/np.pi, 2)} $^\circ$',
                      horizontalalignment='right', verticalalignment='top',
                      transform=axs1.transAxes, rotation=0, color='white', fontsize=l_label_size,
                      bbox=box_style)
        # Set the title of the plot
        fig.suptitle(f'Time range: {t_range[0]} - {t_range[1]}', fontsize=label_size, color='r',
                     y=title_y_pos)

    # fig.show()

    if save_fig:
        try:
            fig_time_range = f"{parser.parse(t_range[0]).strftime('%Y-%m-%d_%H-%M-%S')}_{parser.parse(t_range[1]).strftime('%Y-%m-%d_%H-%M-%S')}"
            fig_name = f'../figures/{fig_name}/ridge_plot_{fig_time_range}.{fig_format}'
            plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=fig_format, dpi=300)
            print(f'Figure saved as {fig_name}')
        except  Exception as e:
            print(e)
            print(f'Figure not saved, folder does not exist. Create folder ../figures')
            #pass
        plt.close()
    return y_vals


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
    im1 = axs1.quiver(x_coord, y_coord, by, bz, scale=scale,
                      scale_units='inches', angles='uv', width=0.002)
    axs1.set_xlabel(r'Y [GSM, $R_\oplus$]', fontsize=label_fontsize)
    axs1.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=label_fontsize)
    patch = patches.Circle(
        (0, 0), radius=15, transform=axs1.transData, fc='none', ec='none', lw=0.1)
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
    dr = args[0][4]
    m_p = args[0][5]
    ro = args[0][6]
    alpha = args[0][7]
    rmp = args[0][8]
    sw_params = args[0][9]
    model_type = args[0][-1]

    y0 = int(j * dr) - y_max
    z0 = int(k * dr) - z_max
    rp = np.sqrt(y0**2 + z0**2)  # Projection of r into yz-plane

    d_theta = np.pi/100

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
            # print( j, k, theta, x_shu[j,k])

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

            rho_sh = sw_params['rho'] * (1.509 * np.exp(x_shu/rmp) + .1285)

            m_proton = 1.672e-27  # Mass of proton in SI unit
            n_sh = rho_sh/m_proton

            y_shu = abs(y_shu)*signy
            z_shu = abs(z_shu)*signz

            # Cooling JGR 2001 Model, equation 9 to 12
            # the distance from the focus to the magnetopause surface
            A = 2
            ll = 3 * rmp/2 - x0
            b_msx = - A * (- sw_params['b_imf'][0] * (1 - rmp / (2 * ll)) + sw_params['b_imf'][1]
                        * (y0 / ll) + sw_params['b_imf'][2] * (z0 / ll))
            b_msy = A * (- sw_params['b_imf'][0] * (y0 / (2 * ll)) + sw_params['b_imf'][1]
                      * (2 - y0**2/( ll * rmp)) - sw_params['b_imf'][2] * (y0 * z0 / (ll * rmp)))
            b_msz = A * (- sw_params['b_imf'][0] * (z0 / (2 * ll)) - sw_params['b_imf'][1]
                      * (y0 * z0 / (ll * rmp)) + sw_params['b_imf'][2] * (2 - z0**2 / (ll * rmp)))
            try:
                if model_type == 't96':
                    bx_ext, by_ext, bz_ext = gp.t96.t96(sw_params['param'], sw_params['ps'], x_shu,
                                                                                       y_shu, z_shu)
                elif model_type == 't01':
                    bx_ext, by_ext, bz_ext = gp.t01.t01(sw_params['param'], sw_params['ps'], x_shu,
                                                                                       y_shu, z_shu)
            except:
                    print(f'Skipped for {x_shu, y_shu, z_shu}')
                    pass

            bx_igrf, by_igrf, bz_igrf = gp.igrf_gsm(x_shu, y_shu, z_shu)

            #print(j, k, bx_ext, bx_igrf)
            bx = bx_ext + bx_igrf
            by = by_ext + by_igrf
            bz = bz_ext + bz_igrf

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


def get_sw_params(
    probe=None,
    omni_level="hro",
    time_clip=True,
    trange=None,
    mms_probe_num=None,
    verbose=False
    ):
    r"""
    Get the solar wind parameters from the OMNI database.

    Parameters
    ----------
    probe : str
        The probe to use. Default is 'None'.
    omni_level : str
        The omni data level to use. Options are 'hro' and 'hro2'. Default is 'hro'.
    time_clip : bool
        If True, the data will be clipped to the time range specified by trange. Default is True.
    trange : list or an array of length 2
        The time range to use. Should in the format [start, end], where start and end times should
        be a string in the format 'YYYY-MM-DD HH:MM:SS'.
    mms_probe_num : str
        The MMS probe to use. Options are '1', '2', '3' and '4'. Default is None.
    verbose : bool
        If True, print out a few messages and the solar wind parameters. Default is False.

    Raises
    ------
    ValueError: If the probe is not one of the options.
    ValueError: If the trange is not in the correct format.

    Returns
    -------
    sw_params : dict
        The solar wind parameters.
    """

    if trange is None:
        raise ValueError("trange must be specified as a list of start and end times in the format 'YYYY-MM-DD HH:MM:SS'.")

    # Check if trange is either a list or an array of length 2
    if not isinstance(trange, (list, np.ndarray)) or len(trange) != 2:
        raise ValueError(
            "trange must be specified as a list or array of length 2 in the format 'YYYY-MM-DD HH:MM:SS.")

    # Download the OMNI data (default level of 'hro_1min') for the specified timerange.
    omni_varnames = ['BX_GSE', 'BY_GSM', 'BZ_GSM', 'proton_density', 'Vx', 'Vy', 'Vz', 'SYM_H']
    omni_vars = spd.omni.data(trange=trange, varnames=omni_varnames, level=omni_level,
                             time_clip=time_clip)

    omni_time = ptt.get_data(omni_vars[0])[0]

    omni_bx_gse = ptt.get_data(omni_vars[0])[1]
    omni_by_gsm = ptt.get_data(omni_vars[1])[1]
    omni_bz_gsm = ptt.get_data(omni_vars[2])[1]
    omni_np = ptt.get_data(omni_vars[3])[1]
    omni_vx = ptt.get_data(omni_vars[4])[1]
    omni_vy = ptt.get_data(omni_vars[5])[1]
    omni_vz = ptt.get_data(omni_vars[6])[1]
    omni_sym_h = ptt.get_data(omni_vars[7])[1]

    # Get mms postion in GSM coordinates for the specified time range

    if (mms_probe_num is not None):
        mms_varnames = [f'mms{mms_probe_num}_mec_r_gsm']
        mms_vars = spd.mms.mec(trange=trange, varnames=mms_varnames, probe=mms_probe_num,
                               data_rate='srvy', level='l2', time_clip=time_clip,
                               latest_version=True)
        mms_time = ptt.get_data(mms_vars[0])[0]
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        mms_sc_pos = ptt.get_data(mms_vars[0])[1:3][0]/r_e
    else:
        mms_time = None
        mms_sc_pos = None
        pass

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
    imf_clock_angle = np.arctan2(b_imf[1], b_imf[2]) * 180 / np.pi
    mean_mms_sc_pos = np.round([np.nanmedian(mms_sc_pos[:,0]), np.nanmedian(mms_sc_pos[:,1]),
                       np.nanmedian(mms_sc_pos[:,2])], 2)
    print("IMF parameters found:")
    if (verbose):
        print(tabulate(
            [["Time of observation (UTC)", time_imf_hrf],
             ["IMF Magnetic field [GSM] (nT)", b_imf],
             ["IMF Proton density (1/cm^-3)", np_imf],
             ["IMF Plasma velocity (km/sec)", v_imf],
             ["IMF clock angle (degrees)", imf_clock_angle],
             ["IMF Sym H", sym_h_imf],
             ["MMS position (GSM) (R_E)", mean_mms_sc_pos]],
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

    # Make a dictionary of all the solar wind parameters
    sw_dict = {}
    sw_dict['time'] = time_imf
    sw_dict['b_imf'] = b_imf
    sw_dict['rho'] = rho
    sw_dict['ps'] = ps
    sw_dict['p_dyn'] = p_dyn
    sw_dict['sym_h'] = sym_h_imf
    sw_dict['imf_clock_angle'] = imf_clock_angle
    sw_dict['param'] = param
    sw_dict['mms_time'] = mms_time
    sw_dict['mms_sc_pos'] = mms_sc_pos

    return sw_dict


def rx_model(
    probe=None,
    trange=['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
    dt=5,
    omni_level = 'hro',
    mms_probe_num = '3',
    model_type = 't96',
    m_p = 0.5,
    dr = 0.5,
    min_max_val = 15,
    y_min = None,
    y_max = None,
    z_min = None,
    z_max = None,
    save_data = False
):
    """
    This function computes the magnetosheath and magnetospheric magnetic fields using the T96 model
    and  sowlarind paramters and returns the value of those fields as well as the value of other
    parameters associated with dayside reconnection, such as shear angle, reconnection energy,
    exhaust velocity, and the bisection field based on bisection model, for both the magnetospheric
    and magnetosheath bisection fields.

    Parameters
    ----------
    probe : str
        The probe to be used for the computation.
        Options: 'mms1', 'mms2', 'mms3', 'mms4'
    trange : list or an array of length 2
        The time range to use. Should in the format [start, end], where start and end times should
        be a string in the format 'YYYY-MM-DD HH:MM:SS'.
    dt : float
        The time interval to use in case the 'trange' has only one element.
    omni_level : str, optional
        The omni-data level to be used for the computation. Default is 'hro'.
        Options: 'hro', 'hro1'
    mms_probe_num : str. optional
        The MMS probe number to be used for the computation. Default is '3'.
        Options: '1', '2', '3', '4'
    model_type : str
        The model type to be used for the computation. Default is 't96'.
        Options: 't96', 't01'
    m_p : float
        Thickness of the magnetosphere in earth radius units . Default is 0.5 R_E.
    dr : float
        Finess of the grid for model computation. Default is 0.5 R_E.
    min_max_val : float
        The minimum and maximum values of the y and z-axis to be used for the computation. Default
        is 15, meaning the model will be computed for -15 < y < 15 and -15 < z < 15.
    y_min : float, optional
        The minimum value of the y-axis to be used for the computation. Default is None, in which
        case y_min is set to -min_max_val.
    y_max : float, optional
        The maximum value of the y-axis to be used for the computation. Default is None, in which
        case y_max is set to min_max_val.
    z_min : float, optional
        The minimum value of the z-axis to be used for the computation. Default is None, in which
        case z_min is set to -min_max_val.
    z_max : float, optional
        The maximum value of the z-axis to be used for the computation. Default is None, in which
        case z_max is set to min_max_val.
    save_data : bool, optional
        Whether to save the data or not. Default is False. If True, the data will be saved in a
        "HDF5" file.
    """

    if y_min is None:
        y_min = -min_max_val
    if y_max is None:
        y_max = min_max_val
    if z_min is None:
        z_min = -min_max_val
    if z_max is None:
        z_max = min_max_val
    
    # If trange has only one element, then make it a list of length 2 with 5 min of padding
    if len(trange) == 1:
        trange_date = datetime.datetime.strptime(trange[0], '%Y-%m-%d %H:%M:%S')
        trange_date_min = trange_date - datetime.timedelta(minutes=dt)
        trange_date_max = trange_date + datetime.timedelta(minutes=dt)
        trange = [trange_date_min.strftime('%Y-%m-%d %H:%M:%S'),
                  trange_date_max.strftime('%Y-%m-%d %H:%M:%S')]

    # Get the solar wind parameters for the model
    sw_params = get_sw_params(probe=probe, omni_level=omni_level, trange=trange,
                              mms_probe_num=mms_probe_num, verbose=True)

    n_arr_y = int((y_max - y_min) / dr) + 1
    n_arr_z = int((z_max - z_min) / dr) + 1
    bx = np.full((n_arr_y, n_arr_z), np.nan)
    by = np.full((n_arr_y, n_arr_z), np.nan)
    bz = np.full((n_arr_y, n_arr_z), np.nan)
    b_msx = np.full((n_arr_y, n_arr_z), np.nan)
    b_msy = np.full((n_arr_y, n_arr_z), np.nan)
    b_msz = np.full((n_arr_y, n_arr_z), np.nan)
    shear = np.full((n_arr_y, n_arr_z), np.nan)
    rx_en = np.full((n_arr_y, n_arr_z), np.nan)
    va_cs = np.full((n_arr_y, n_arr_z), np.nan)
    bisec_msp = np.full((n_arr_y, n_arr_z), np.nan)
    bisec_msh = np.full((n_arr_y, n_arr_z), np.nan)

    # Shue et al.,1998, equation 10
    ro = (10.22 + 1.29 * np.tanh(0.184 * (sw_params['b_imf'][2] + 8.14))) * (
                                          sw_params['p_dyn'])**(-1.0/6.6)

    # Shue et al.,1998, equation 11
    alpha = (0.58 - 0.007 * sw_params['b_imf'][2]) * (1 + 0.024 * np.log(sw_params['p_dyn']))
    rmp = ro * (2/(1 + np.cos(0.0))) ** alpha  # Stand off position of the magnetopause

    len_y = int((y_max - y_min)/dr) + 1
    len_z = int((z_max - z_min)/dr) + 1

    p = mp.Pool()

    input = ((j, k, y_max, z_max, dr, m_p, ro, alpha, rmp, sw_params, model_type)
             for j in range(len_y) for k in range(len_z))

    print("Running the model \n")
    res = p.map(model_run, input)
    print("Model run complete \n")

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
            today_date = datetime.datetime.today().strftime('%Y-%m-%d')
            fn = f'../data/all_data_rx_model_{dr}re_{m_p}mp_{model_type}_{today_date}.h5'
            data_file = hf.File(fn, 'w')

            data_file.create_dataset('bx', data=bx)
            data_file.create_dataset('by', data=by)
            data_file.create_dataset('bz', data=bz)

            data_file.create_dataset('b_msx', data=b_msx)
            data_file.create_dataset('b_msy', data=b_msy)
            data_file.create_dataset('b_msz', data=b_msz)

            data_file.create_dataset('shear', data=shear)
            data_file.create_dataset('rx_en', data=rx_en)
            data_file.create_dataset('va_cs', data=va_cs)
            data_file.create_dataset('bisec_msp', data=bisec_msp)
            data_file.create_dataset('bisec_msh', data=bisec_msh)

            data_file.close()
            print(f'Date saved to file {fn} \n')
        except Exception as e:
            print(e)
            print(
                f'Data not saved to file {fn}. Please make sure that file name is correctly assigned and that the directory exists and you have write permissions')

    return bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh, sw_params
