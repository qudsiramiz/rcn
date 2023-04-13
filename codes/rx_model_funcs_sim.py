# functions to be used in rx_model_batch_parallel.py
import datetime
import multiprocessing as mp
import os

import geopack.geopack as gp
import h5py as hf
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from matplotlib.pyplot import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import frangi


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

    # alpha = - 14.87 * np.pi / 180  # radians (From Hesse2013)
    b_vec_1 = np.array(b_vec_1)
    b_vec_2 = np.array(b_vec_2)
    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1 / mag_vec_1
    unit_vec_2 = b_vec_2 / mag_vec_2

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

    unit_vec_1 = b_vec_1 / mag_vec_1
    unit_vec_2 = b_vec_2 / mag_vec_2

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

    unit_vec_1 = b_vec_1 / mag_vec_1
    unit_vec_2 = b_vec_2 / mag_vec_2

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
    angle = np.arctan(b_vec[1] / b_vec[2])

    if (angle_unit == "radians"):
        return angle
    elif (angle_unit == "degrees"):
        return angle * 180 / np.pi
    else:
        raise KeyError("angle_unit must be radians or degrees")


def ridge_finder_multiple(
    image=[None, None, None, None],
    convolution_order=[1, 1, 1, 1],
    t_range=['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
    dt=5,
    b_imf=[-5, 0, 0],
    b_msh=[-5, 0, 0],
    v_msh=[-200, 50, 50],
    xrange=[-15.1, 15],
    yrange=[-15.1, 15],
    mms_probe_num='1',
    mms_sc_pos=[0, 0],
    dr=0.5,
    dipole_tilt_angle=None,
    p_dyn=None,
    imf_clock_angle=None,
    sigma=[2.2, 2.2, 2.2, 2.2],
    mode="nearest",
    alpha=1.,
    vmin=[None, None, None, None],
    vmax=[None, None, None, None],
    cmap_list=["viridis", "viridis", "viridis", "viridis"],
    draw_patch=[True, True, True, True],
    draw_ridge=[False, False, False, False],
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
    tsy_model="t96",
    dark_mode=True,
    rc_file_name="rc_file.csv",
    rc_folder="../data",
    save_rc_file=False,
    walen1=False,
    walen2=False,
    jet_detection=False,
    fig_version='v6',
    r_W=None,
    theta_W=None,
    jet_time=None,
    np_median_msp=None,
    np_median_msh=None,
    df_jet_reversal=None,
):
    r"""
    Finds ridges in an image and plot the points with maximum ridge value on the given image.

    Parameters
    ----------
    image : list of numpy arrays
        List of images to be plotted.
    convolution_order : list of ints
        List of the order of the convolution to be used while smoothing each image. Values must be
        non-negative integers. Default is [1, 1, 1, 1].
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
    dark_mode : bool, optional
        Sets the dark mode for the plot and adjusts the color of labels and tickmarks accordingly.
        Default is True.

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
        # Check if t_range is a datetime object
        if isinstance(t_range[0], datetime.datetime):
            t_range_date = t_range[0]
        else:
            t_range_date = datetime.datetime.strptime(t_range[0], '%Y-%m-%d %H:%M:%S')
        t_range_date_min = t_range_date - datetime.timedelta(minutes=dt)
        t_range_date_max = t_range_date + datetime.timedelta(minutes=dt)
        t_range = [t_range_date_min.strftime('%Y-%m-%d %H:%M:%S'),
                   t_range_date_max.strftime('%Y-%m-%d %H:%M:%S')]

    if dark_mode:
        plt.style.use('dark_background')
        # tick_color = 'w'  # color of the tick lines
        mtick_color = 'w'  # color of the minor tick lines
        label_color = 'w'  # color of the tick labels
        clabel_color = 'w'  # color of the colorbar label
    else:
        plt.style.use('default')
        # tick_color = 'k'  # color of the tick lines
        mtick_color = 'k'  # color of the minor tick lines
        label_color = 'k'  # color of the tick labels
        clabel_color = 'k'  # color of the colorbar label

    # Set the fontstyle to Times New Roman
    font = {'family': 'serif', 'weight': 'normal', 'size': 10}
    plt.rc('font', **font)
    plt.rc('text', usetex=True)
    fig = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor='w', edgecolor='k')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=wspace, hspace=hspace)
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1])

    # Set the font size for the axes
    label_size = 18  # fontsize for x and y labels
    t_label_size = 18  # fontsize for tick label
    c_label_size = 18  # fontsize for colorbar label
    ct_tick_size = 14  # fontsize for colorbar tick labels
    l_label_size = 14  # fontsize for legend label

    tick_len = 10  # length of the tick lines
    mtick_len = 7  # length of the minor tick lines
    tick_width = 1  # tick width in points
    mtick_width = 0.7  # minor tick width in points

    # box_style = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    if dark_mode:
        box_style = box_style
    else:
        box_style = dict(boxstyle="round", color="w", alpha=0.8, linewidth=1)
    y_vals = []
    x_intr_vals_list = []
    y_intr_vals_list = []
    for i in range(len(image)):
        image_rotated = np.transpose(image[i])

        # Create the masked image from result for all the new processings
        # Find the number of rows in the original image
        n_rows, n_cols = image_rotated.shape

        # Make a grid of the data based on mumber of rows and columns
        X, Y = np.ogrid[:n_rows, :n_cols]

        # Find the central row and column
        c_row = int(n_rows/2)
        c_col = int(n_cols/2)
        # Find the distance of each pixel from the central pixel in terms of pixels
        dist_pxl = np.sqrt((X - c_row) ** 2 + (Y - c_col) ** 2)
        mask_image = dist_pxl > 15 / dr

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
        image_smooth = sp.ndimage.gaussian_filter(image_rotated, order=convolution_order[i],
                                                  sigma=[5, 5], mode=mode)
        image_smooth_p = sp.ndimage.gaussian_filter(image_rotated, order=0, sigma=[5, 5],
                                                    mode=mode)
        result = frangi(image_smooth, **kwargs)  # frangi, hessian, meijering, sato

        m_result = result.copy()
        m_result[mask_image] = np.nan
        new_image_rotated = image_rotated.copy()
        new_image_rotated[mask_image] = np.nan

        x_len = image_rotated.shape[0]
        y_len = image_rotated.shape[1]

        y_val = np.full(y_len, np.nan)
        y_vals.append(y_val)
        im_max_val = np.full(y_len, np.nan)
        for xx in range(y_len):
            try:
                y_val[xx] = np.nanargmax(m_result[:, xx]) * dr + yrange[0]
                im_max_val[xx] = np.nanargmax(new_image_rotated[:, xx]) * dr + yrange[0]
            except Exception:
                pass

        axs1 = plt.subplot(gs[0, i])
        im1 = axs1.imshow(image_smooth_p, extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                          origin='lower', cmap=cmap_list[i], norm=norm, interpolation=interpolation,
                          alpha=1)
        divider1 = make_axes_locatable(axs1)
        # Draw a circle of radius 10 around the center of the image
        axs1.add_patch(plt.Circle((0, 0), radius=15, color='gray', fill=False, lw=0.5))

        # Take rolling average of the y_val array
        y_val_avg = np.full(len(y_val), np.nan)
        im_max_val_avg = np.full(len(y_val), np.nan)

        r_a_l = 5
        for xx in range(len(y_val)):
            y_val_avg[xx] = np.nanmean(y_val[max(0, xx - r_a_l):min(len(y_val), xx + r_a_l)])
            im_max_val_avg[xx] = np.nanmean(im_max_val[max(0, xx - r_a_l):min(len(y_val),
                                            xx + r_a_l)])

        axs1.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
        axs1.axvline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)

        if draw_ridge:
            # axs1.plot(np.linspace(xrange[0], xrange[1], x_len), y_val_avg, color='aqua', ls='-',
            #          alpha=0.9)
            x_intr_vals = np.linspace(xrange[0], xrange[1], x_len)
            y_intr_vals = im_max_val_avg
            # If the square root of the sum of squares of x_intr_vals and y_intr_vals is greater
            # than 15, then mask those values
            # r_intr_vals = np.sqrt(x_intr_vals ** 2 + y_intr_vals ** 2)
            # mask = r_intr_vals > 15
            # Mask the values of x_intr_vals and y_intr_vals
            # x_intr_vals[mask] = np.nan
            # y_intr_vals[mask] = np.nan
            # if z component of b_imf is negative, then the ridge is on the left side of the
            # image
            #if b_imf[2] <= 0:
            #    axs1.plot(x_intr_vals, y_intr_vals, color='aqua', ls='-', alpha=0.9)
        # Plot a horizontal line at x=0 and a vertical line at y=0
        if (draw_patch):
            patch = patches.Circle((0, 0), radius=xrange[1], transform=axs1.transData, fc='none',
                                   ec='k', lw=0.5)
            im1.set_clip_path(patch)
        axs1.add_patch(patch)
        if i == 0 or i == 3:
            axs1.set_ylabel(r'Z [GSM, $R_{\rm E}$]', fontsize=label_size, color=label_color)
        if i == 3:
            axs1.yaxis.set_label_position("right")

        axs1.set_xlabel(r'Y [GSM, $R_{\rm E}$]', fontsize=label_size, color=label_color)
        if dark_mode:
            text_color = 'white'
        else:
            text_color = 'black'
        if i == 0:
            # axs1.text(-0.3, 1.16, f'Model: {tsy_model}', horizontalalignment='left',
            #           verticalalignment='bottom', transform=axs1.transAxes, rotation=0,
            #           color=text_color, fontsize=l_label_size, bbox=box_style)
            axs1.text(-0.3, 1.16, f'Clock Angle: {np.round(imf_clock_angle, 2)}$^\\circ$',
                      horizontalalignment='left', verticalalignment='bottom', transform=axs1.transAxes,
                      rotation=0, color=text_color, fontsize=l_label_size, bbox=box_style)

        if i == 3:
            axs1.text(1.3, 1.16, f'$B_{{\\rm {{imf}}}}$ = [{b_imf[0]}, {b_imf[1]}, {b_imf[2]}]',
                      horizontalalignment='right',
                      verticalalignment='bottom', transform=axs1.transAxes, rotation=0,
                      color=text_color, fontsize=l_label_size, bbox=box_style)
        # elif i == 3:
        #     axs1.text(1.3, -0.15,
        #               f'Dipole tilt: {np.round(dipole_tilt_angle * 180/np.pi, 2)} ${{\\hspace{{-.2em}}}}^\\circ$',
        #               horizontalalignment='right', verticalalignment='top',
        #               transform=axs1.transAxes, rotation=0, color=text_color, fontsize=l_label_size,
        #               bbox=box_style)

        # Define the location of the colorbar, it's size relative to main figure and the padding
        # between the colorbar and the figure, the orientation the colorbar
        cax1 = divider1.append_axes("top", size="5%", pad=0.01)
        cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05,
                             pad=0.01)
        cbar1.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                             labelbottom=False, pad=0.01, labelsize=ct_tick_size,
                             labelcolor=label_color)
        # Get the location of all ticks on the colorbar
        cbar_ticks = cbar1.ax.get_xticks()
        # Remove the first tick
        cbar_ticks = cbar_ticks[1:]
        # Set the ticks to the new tick values
        cbar1.ax.set_xticks(cbar_ticks)

        cbar1.ax.xaxis.set_label_position('top')

        cbar1.ax.set_xlabel(f'{c_label[i]}', fontsize=c_label_size,
                            color=clabel_color)

        # Set tick label parameters
        if i == 0:
            axs1.tick_params(axis='both', direction='in', which='major', left=True, right=True,
                             top=True, bottom=True, labelleft=True, labelright=False,
                             labeltop=False, labelbottom=True, labelsize=t_label_size,
                             length=tick_len, width=tick_width, labelcolor=label_color)
        elif i == 1 or i == 2:
            axs1.tick_params(axis='both', direction='in', which='major', left=True, right=True,
                             top=True, bottom=True, labelleft=False, labelright=False,
                             labeltop=False, labelbottom=True, labelsize=t_label_size,
                             length=tick_len, width=tick_width, labelcolor=label_color)
        else:
            axs1.tick_params(axis='both', direction='in', which='major', left=True, right=True,
                             top=True, bottom=True, labelleft=False, labelright=True,
                             labeltop=False, labelbottom=True, labelsize=t_label_size,
                             length=tick_len, width=tick_width, labelcolor=label_color)

        if i == 0:
            # Add a label '(a)' to the plot to indicate the panel number
            axs1.text(0.05, 0.15, '(a)', horizontalalignment='left', verticalalignment='top',
                      transform=axs1.transAxes, rotation=0, color=text_color,
                      fontsize=1.2 * l_label_size)
        elif i == 1:
            # Add a label '(b)' to the plot to indicate the panel number
            axs1.text(0.05, 0.15, '(b)', horizontalalignment='left', verticalalignment='top',
                      transform=axs1.transAxes, rotation=0, color=text_color,
                      fontsize=1.2 * l_label_size)
        elif i == 2:
            # Add a label '(c)' to the plot to indicate the panel number
            axs1.text(0.05, 0.15, '(c)', horizontalalignment='left', verticalalignment='top',
                      transform=axs1.transAxes, rotation=0, color=text_color,
                      fontsize=1.2 * l_label_size)
        elif i == 3:
            # Add a label '(d)' to the plot to indicate the panel number
            axs1.text(0.05, 0.15, '(d)', horizontalalignment='left', verticalalignment='top',
                      transform=axs1.transAxes, rotation=0, color=text_color,
                      fontsize=1.2 * l_label_size)

        # Show minor ticks
        axs1.minorticks_on()
        axs1.tick_params(axis='both', which='minor', direction='in', length=mtick_len, left=True,
                         right=True, top=True, bottom=True, color=mtick_color, width=mtick_width)
        # Set the number of ticks on the x-axis
        axs1.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
        # Set the number of ticks on the y-axis
        axs1.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))

        # Setting the tickmarks labels in such a way that they don't overlap
        plt.setp(axs1.get_xticklabels(), rotation=0, ha='right', va='top', visible=True)
        plt.setp(axs1.get_yticklabels(), rotation=0, va='center', visible=True)
        # Set the title of the plot
        # fig.suptitle(f'$B_{{\\rm {{imf}}}}$ = {b_imf}',
        #              fontsize=label_size, color=text_color, y=title_y_pos, alpha=0.65)

    #plt.show()
    if save_fig:
        try:
            # TODO: Add folder name as one of the path and make sure that the code creates the
            # folder. Gives out error if the folder can't be created.
            fig_folder = f"../figures/test/{tsy_model}/{interpolation}" +\
                         f"_interpolation_mms{mms_probe_num}/{fig_version}"
            check_folder = os.path.isdir(fig_folder)
            # If folder doesn't exist, then create it.
            if not check_folder:
                os.makedirs(fig_folder)
                print("created folder : ", fig_folder)
            else:
                print(f"folder already exists: {fig_folder}\n")

            # fig_folder = "../figures/test"
            fig_name = f'{fig_folder}/ridge_plot_{int(b_imf[0])}_{int(b_imf[1])}_{int(b_imf[2])}.{fig_format}'
            plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=fig_format, dpi=200)
            print(f'Figure saved as {fig_name}')
        except Exception as e:
            print(e)
            print('Figure not saved, folder does not exist. Create folder ../figures')
            # pass
        # plt.close()
    plt.close()
    return y_vals, x_intr_vals_list, y_intr_vals_list


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

    d_theta = np.pi / 100

    for index in range(0, 100):

        theta = index * d_theta
        r = ro * (2 / (1 + np.cos(theta))) ** alpha
        zp = r * np.sin(theta)  # not really in z direction, but a distance in yz plane
        x0 = r * np.cos(theta)

        # if x0 == 0:
        #     signx = 1.0
        # else:
        #     signx = np.sign(x0)

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

            # y_coord = y0
            # z_coord = z0
            x_shu = (r - m_p) * np.cos(theta)
            phi = np.arctan2(z0, y0)
            # print( j, k, theta, x_shu[j,k])

            if (abs(y0) == 0 or abs(z0) == 0):
                if (abs(y0) == 0):
                    y_shu = 0
                    z_shu = (r - m_p) * np.sin(theta)
                elif (abs(z0) == 0):
                    z_shu = 0
                    y_shu = (r - m_p) * np.sin(theta)
            else:
                z_shu = np.sqrt((rp - 1.0)**2 / (1 + np.tan(phi)**(-2)))
                y_shu = z_shu / np.tan(phi)

            m_proton = 1.672e-27  # Mass of proton in SI unit
            n_sh = sw_params['rho'] * (1.509 * np.exp(x_shu / rmp) + .1285) / m_proton

            y_shu = abs(y_shu) * signy
            z_shu = abs(z_shu) * signz

            # Cooling JGR 2001 Model, equation 9 to 12
            # the distance from the focus to the magnetopause surface
            A = 2
            ll = 3 * rmp / 2 - x0
            b_msx = - A * (- sw_params['b_imf'][0] * (1 - rmp / (2 * ll)) + sw_params['b_imf'][1]
                           * (y0 / ll) + sw_params['b_imf'][2] * (z0 / ll))
            b_msy = A * (- sw_params['b_imf'][0] * (y0 / (2 * ll)) + sw_params['b_imf'][1]
                         * (2 - y0**2 / (ll * rmp)) - sw_params['b_imf'][2] * (y0 * z0 / (ll *
                                                                                          rmp)))
            b_msz = A * (- sw_params['b_imf'][0] * (z0 / (2 * ll)) - sw_params['b_imf'][1]
                         * (y0 * z0 / (ll * rmp)) + sw_params['b_imf'][2] * (2 - z0**2 / (ll *
                                                                                          rmp)))
            try:
                if model_type == 't96':
                    bx_ext, by_ext, bz_ext = gp.t96.t96(sw_params['param'], sw_params['ps'], x_shu,
                                                        y_shu, z_shu)
                elif model_type == 't01':
                    bx_ext, by_ext, bz_ext = gp.t01.t01(sw_params['param'], sw_params['ps'], x_shu,
                                                        y_shu, z_shu)
            except Exception:
                print(f'Skipped for {x_shu, y_shu, z_shu}')
                pass

            bx_igrf, by_igrf, bz_igrf = gp.igrf_gsm(x_shu, y_shu, z_shu)

            # print(j, k, bx_ext, bx_igrf)
            bx = bx_ext + bx_igrf
            by = by_ext + by_igrf
            bz = bz_ext + bz_igrf

            shear = get_shear([bx, by, bz], [b_msx, b_msy, b_msz], angle_unit="degrees")

            rx_en = get_rxben([bx, by, bz], [b_msx, b_msy, b_msz])
            va_cs = get_vcs([bx, by, bz], [b_msx, b_msy, b_msz], n_sh, 0.1)
            bisec_msp, bisec_msh = get_bis([bx, by, bz], [b_msx, b_msy, b_msz])
            break

    return (j, k, bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh, x_shu, y_shu, z_shu, b_msx,
            b_msy, b_msz)


def rx_model(
    probe=None,
    trange=['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
    sw_params=None,
    dt=5,
    omni_level='hro',
    mms_probe_num='3',
    model_type='t96',
    m_p=0.5,
    dr=0.5,
    min_max_val=15,
    y_min=None,
    y_max=None,
    z_min=None,
    z_max=None,
    save_data=False,
    nprocesses=None
):
    """
    This function computes the magnetosheath and magnetospheric magnetic fields using the T96 model
    and  solar wind paramters and returns the value of those fields as well as the value of other
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
    nprocesses : int, optional
        The number of processes to use for the computation. Default is None, in which case the
        number of processes will be set to the number of cores in the system.
    """

    if y_min is None:
        y_min = -min_max_val
    if y_max is None:
        y_max = min_max_val
    if z_min is None:
        z_min = -min_max_val
    if z_max is None:
        z_max = min_max_val

    # If trange has only one element, then make it a list of length 2 with 'dt' minutes of padding
    if len(trange) == 1:
        # Check if trange is a datetime object
        if isinstance(trange[0], datetime.datetime):
            trange_date = trange[0]
        else:
            trange_date = datetime.datetime.strptime(trange[0], '%Y-%m-%d %H:%M:%S')
        trange_date_min = trange_date - datetime.timedelta(minutes=dt)
        trange_date_max = trange_date + datetime.timedelta(minutes=dt)
        trange = [trange_date_min.strftime('%Y-%m-%d %H:%M:%S'),
                  trange_date_max.strftime('%Y-%m-%d %H:%M:%S')]
        trange = [trange[0] + 'Z', trange[1] + 'Z']

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

    x_shu = np.full((n_arr_y, n_arr_z), np.nan)
    y_shu = np.full((n_arr_y, n_arr_z), np.nan)
    z_shu = np.full((n_arr_y, n_arr_z), np.nan)

    b_msx = np.full((n_arr_y, n_arr_z), np.nan)
    b_msy = np.full((n_arr_y, n_arr_z), np.nan)
    b_msz = np.full((n_arr_y, n_arr_z), np.nan)

    # Shue et al.,1998, equation 10
    ro = (10.22 + 1.29 * np.tanh(0.184 * (sw_params['b_imf'][2] + 8.14))) * (
        sw_params['p_dyn'])**(-1.0 / 6.6)

    # Shue et al.,1998, equation 11
    alpha = (0.58 - 0.007 * sw_params['b_imf'][2]) * (1 + 0.024 * np.log(sw_params['p_dyn']))
    rmp = ro * (2 / (1 + np.cos(0.0))) ** alpha  # Stand off position of the magnetopause

    len_y = int((y_max - y_min) / dr) + 1
    len_z = int((z_max - z_min) / dr) + 1

    if nprocesses is None:
        p = mp.Pool()
    else:
        p = mp.Pool(processes=nprocesses)

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

        x_shu[j, k] = r[10]
        y_shu[j, k] = r[11]
        z_shu[j, k] = r[12]

        b_msx[j, k] = r[13]
        b_msy[j, k] = r[14]
        b_msz[j, k] = r[15]

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

            data_file.create_dataset('x_shu', data=x_shu)
            data_file.create_dataset('y_shu', data=y_shu)
            data_file.create_dataset('z_shu', data=z_shu)

            data_file.close()
            print(f'Date saved to file {fn} \n')
        except Exception as e:
            print(e)
            print(f'Data not saved to file {fn}. Please make sure that file name is correctly' +
                  ' assigned and that the directory exists and you have write permissions')

    return (bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh, sw_params, x_shu, y_shu, z_shu,
            b_msx, b_msy, b_msz)


def line_fnc(
    r0=np.array([0, 0, 0]),
    b_msh=np.array([0, 0, 0]),
    r=0,
):
    """
    Function to compute the equation of a line along the direction of the magnetic field.

    Parameters
    ----------
    r0 : array
        Starting position of the line
    b_msh : array, optional
        Magnetic field from magnetosheath at the magnetopause
    r : float, optional
        Position of the point where the line is computed
    """
    # Compute the direction of the magnetic field line
    b_msh_norm = np.linalg.norm(b_msh)
    b_msh_dir = b_msh / b_msh_norm

    # return the line function
    return r0 + r * b_msh_dir


def line_fnc_der(x, y):
    """
    Function to give a line interpolation function

    Parameters
    ----------
    x : array
        x coordinates of the points
    y : array
        y coordinates of the points

    Returns
    -------
    line_intrp : function
        Function to interpolate the line
    """

    nans, x_nan = nan_helper(x)
    x[nans] = np.interp(x_nan(nans), x_nan(~nans), x[~nans])

    nans, y_nan = nan_helper(y)
    y[nans] = np.interp(y_nan(nans), y_nan(~nans), y[~nans])

    line_intrp = sp.interpolate.CubicSpline(x, y)
    return line_intrp


def target_fnc(r, r0, b_msh, line_fnc, line_intrp):
    p_line = line_fnc(r0=r0, b_msh=b_msh, r=r)
    z_surface = line_intrp(p_line[1])
    return np.sum((p_line[2] - z_surface)**2)


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]
