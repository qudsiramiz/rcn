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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import meijering  # sato, frangi, hessian

# Set the fontstyle to Times New Roman
font = {'family': 'sans-serif', 'weight': 'normal', 'size': 10}
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
        Angle between the two vectors in radians by default
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
    # TODO: Update the documentation of this function

    b_vec_1 = np.array(b_vec_1)
    b_vec_2 = np.array(b_vec_2)
    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    bisector = mag_vec_2 * b_vec_1 + mag_vec_1 * b_vec_2
    #u_bisect = bisector/np.linalg.norm(bisector)
    u_bisect = (unit_vec_1 + unit_vec_2)/np.sqrt(2)
    rx_bmag1 = np.dot(u_bisect, b_vec_1)
    rx_bmag2 = np.dot(u_bisect, b_vec_2)
    b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)

    #rx_en = 0.5 * (mag_vec_1 * mag_vec_2) * (1 + b1_b2_dotp)
    rx_en = (rx_bmag1 - rx_bmag2)**2 * 3.98e-4  #nPa
    #rx_en = (rx_bmag1**2 + rx_bmag2**2) * 1.03  # MJ/RE^3

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


def get_bis(b_vec_1, b_vec_2, angle_unit="radians"):
    r"""
    Get the shear angle between two magnetic field lines.

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        First input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Second input magnetic field vector.
    angle_unit : str, optional
        Preferred unit of angle returned by the code. Default is "radians".
    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
     Returns bisect stuff
    """
    b_vec_1 = np.array(b_vec_1)
    b_vec_2 = np.array(b_vec_2)

    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)

    angle = np.arccos(b1_b2_dotp)

    bisector = mag_vec_1*b_vec_1 + mag_vec_2*b_vec_2
    u_bisect = bisector/np.linalg.norm(bisector)

    bis_theta = np.arccos(np.dot(u_bisect, unit_vec_1))
    rx_mag_1 = np.dot(u_bisect, b_vec_1)
    rx_mag_2 = np.dot(u_bisect, b_vec_2)

    bis_field = mag_vec_2 * np.sin(bis_theta)
    # bis_field = rx_mag_1**2 + rx_mag_2**2

    #if (angle_unit == "radians"):
    #    return angle
    #elif (angle_unit == "degrees"):
    #    return angle * 180/np.pi
    #else:
    #    raise KeyError("angle_unit must be radians or degrees")
    return bis_field


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
    xrange=[-15.1, 15],
    yrange=[-15.1, 15],
    dr=0.5,
    sigma=2.2,
    mode="nearest",
    alpha=0.5,
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
    save_fig : bool, optional
            Whether to save the figure. Default is True.
    fig_name : str, optional
            The name of the figure. Default is "new".
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
    image = np.transpose(image)

    cmap = plt.cm.viridis

    kwargs = {'sigmas': [sigma], 'black_ridges': False, 'mode': mode, 'alpha': alpha}

    result = meijering(image, **kwargs)

    x_len = image.shape[0]
    y_len = image.shape[1]
    y_val = np.full(y_len, np.nan)
    im_max_val = np.full(y_len, np.nan)
    for i in range(y_len):
        y_val[i] = np.argmax(result[:, i]) * dr + yrange[0]
        im_max_val[i] = np.argmax(image[:, i]) * dr + yrange[0]

    # plt.close('all')
    fig, axs1 = plt.subplots(1, 1, figsize=(8, 6))

    im1 = axs1.imshow(abs(image), extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                      origin='lower', cmap=cmap)
    divider1 = make_axes_locatable(axs1)

    axs1.plot(np.linspace(xrange[0], xrange[1], x_len), y_val, 'k-', ms=2)

    axs1.plot(np.linspace(xrange[0], xrange[1], x_len), im_max_val, 'r*', ms=2)

    patch = patches.Circle((0, 0), radius=15, transform=axs1.transData, fc='none', ec='k', lw=0.1)
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

    # fig.show()

    if save_fig:
        fig_name = f'../figures/ridge_plot_vir_{fig_name}_{dr}dr_{mp}mp.{fig_format}'
        plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=fig_format, dpi=300)
        print(f'Figure saved as {fig_name}')
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


#def rx_model_batch(
#    probe=None,
#    omni_level='hro',
#    maximum_shear=True,
#    mms_probe='mms1',
#    movie=None,
#    trange=None,
#    model_type='t96',
#    mp=0.5,
#    dr=0.25
#    save_data=False,
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
    mp : float, optional
        Thickness of the magnetopause. The default is 0.5.
    dr : float, optional
        Resolution of the model. The default is 0.25.

    Returns
    -------
    None.
    """
    probe = None
    omni_level = 'hro'
    maximum_shear = True
    mms_probe = None
    movie = None
    trange = ['2016-12-24 15:08:00', '2016-12-24 15:12:00']
    model_type = 't96'
    mp = 0.5  # Magnetopause thichkness
    dr = 0.5  # Resolution of model run in R_E units
    save_data = True

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
    omni_vars = spd.omni.data(trange=trange, level=omni_level)

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
        mms_time = ptt.get_data('mms1_mec_r_gse')[0]
        mms_sc_pos = ptt.get_data('mms1_mec_r_gse')[1:3]
    else:
        pass

    if (movie is None):
        time_imf = np.nanmedian(omni_time)
        b_imf_x = np.nanmedian(omni_bx_gse)
        b_imf_y = np.nanmedian(omni_by_gsm)
        b_imf_z = np.nanmedian(omni_bz_gsm)

        if (b_imf_z > 15 or b_imf_z < -18):
            warnings.warn(
            f"The given parameters produced the z-component of IMF field (b_imf_z) {b_imf_z} nT,"
            f"which is out of range in which model is valid (-18 nT < b_imf_z < 15 nT)"
            )

        np_imf = np.nanmedian(omni_np)
        vx_imf = np.nanmedian(omni_vx)
        vy_imf = np.nanmedian(omni_vy)
        vz_imf = np.nanmedian(omni_vz)
        sym_h_imf = np.nanmedian(omni_sym_h)

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

        m_p = 1.672e-27  # Mass of proton in SI unit

        rho = np_imf * m_p * 1.15

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

        n_arr = int(30/dr) + 1

        bx = np.full((n_arr, n_arr), np.nan)
        by = np.full((n_arr, n_arr), np.nan)
        bz = np.full((n_arr, n_arr), np.nan)

        bx_ext = np.full((n_arr, n_arr), np.nan)
        by_ext = np.full((n_arr, n_arr), np.nan)
        bz_ext = np.full((n_arr, n_arr), np.nan)

        bx_igrf = np.full((n_arr, n_arr), np.nan)
        by_igrf = np.full((n_arr, n_arr), np.nan)
        bz_igrf = np.full((n_arr, n_arr), np.nan)

        b_msx = np.full((n_arr, n_arr), np.nan)
        b_msy = np.full((n_arr, n_arr), np.nan)
        b_msz = np.full((n_arr, n_arr), np.nan)

        x_shu = np.full((n_arr, n_arr), np.nan)
        y_shu = np.full((n_arr, n_arr), np.nan)
        z_shu = np.full((n_arr, n_arr), np.nan)

        rho_sh = np.full((n_arr, n_arr), np.nan)

        # rp = np.full((n_arr, n_arr), np.nan)

        # r = np.full((n_arr, n_arr), np.nan)
        # zp = np.full((n_arr, n_arr), np.nan)
        # x0 = np.full((n_arr, n_arr), np.nan)

        shear = np.full((n_arr, n_arr), np.nan)
        rx_en = np.full((n_arr, n_arr), np.nan)
        va_cs = np.full((n_arr, n_arr), np.nan)
        bisec = np.full((n_arr, n_arr), np.nan)
        y_coord = np.full((n_arr, n_arr), np.nan)
        z_coord = np.full((n_arr, n_arr), np.nan)
        b_sh_ca = np.full((n_arr, n_arr), np.nan)
        b_sh_mag = np.full((n_arr, n_arr), np.nan)
        n_sh = np.full((n_arr, n_arr), np.nan)

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
        len_y = int(30/dr) + 1
        len_z = int(30/dr) + 1
        count = 0
        for j in range(0, len_y):
            y0 = 15 - int(j * dr)
            for k in range(0, len_z):
                z0 = 15 - int(k * dr)
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

                        y_coord[j, k] = y0
                        z_coord[j, k] = z0
                        x_shu[j, k] = (r - mp) * np.cos(theta)
                        phi = np.arctan2(z0, y0)
                        # print( j, k, theta, x_shu[j,k])

                        if (abs(y0) == 0 or abs(z0) == 0):
                            if(abs(y0) == 0):
                                y_shu[j, k] = 0
                                z_shu[j, k] = (r - mp) * np.sin(theta)
                            elif (abs(z0) == 0):
                                z_shu[j, k] = 0
                                y_shu[j, k] = (r - mp) * np.sin(theta)
                        else:
                            z_shu[j, k] = np.sqrt((rp - 1.0)**2/(1 + np.tan(phi)**(-2)))
                            y_shu[j, k] = z_shu[j, k]/np.tan(phi)

                        rho_sh[j, k] = rho * (1.509 * np.exp(x_shu[j, k]/rmp) + .1285)
                        n_sh[j, k] = rho_sh[j, k]/m_p

                        y_shu[j, k] = abs(y_shu[j, k])*signy
                        z_shu[j, k] = abs(z_shu[j, k])*signz

                        # Cooling JGR 2001 Model, equation 9 to 12
                        # the distance from the focus to the magnetopause surface
                        ll = 3 * rmp/2 - x0
                        b_msx[j, k] = - A * (- b_imf_x * (1 - rmp / (2 * ll)) + b_imf_y * (y0 / ll)
                                      + b_imf_z * (z0 / ll))
                        b_msy[j, k] = A * (- b_imf_x * (y0 / (2 * ll)) + b_imf_y * (2 - y0**2/(
                                      ll * rmp)) - b_imf_z * (y0 * z0 / (ll * rmp)))
                        b_msz[j, k] = A * (- b_imf_x * (z0 / (2 * ll)) - b_imf_y * (y0 * z0 / (ll
                                      * rmp)) + b_imf_z * (2 - z0**2 / (ll * rmp)))

                        # TODO: Implement Geopack T96!!!
                        # Compute the external magnetic field from the T95 model for a given
                        # position and time, in GSM coordinate
                        try:
                            if(model_type == 't96'):
                                bx_ext[j, k], by_ext[j, k], bz_ext[j, k] = gp.t96.t96(
                                                                                      param,
                                                                                      ps,
                                                                                      x_shu[j, k],
                                                                                      y_shu[j, k],
                                                                                      z_shu[j, k]
                                                                                      )
                            elif(model_type == 't01'):
                                bx_ext[j, k], by_ext[j, k], bz_ext[j, k] = gp.t01.t01(
                                                                                      param,
                                                                                      ps,
                                                                                      x_shu[j, k],
                                                                                      y_shu[j, k],
                                                                                      z_shu[j, k]
                                                                                      )
                            else:
                                raise ValueError("Model type must be set to 't96' or 't01'.")
                        except:
                            if((y_shu[j,k] < 15 and y_shu[j, k] > -15) or
                               (z_shu[j,k] < 15 and z_shu[j, k] > -15)):
                                print(f'Skipped for {x_shu[j, k], y_shu[j, k], z_shu[j, k]}')
                                count += 1
                            else:
                                print(f'~Skipped for {x_shu[j, k], y_shu[j, k], z_shu[j, k]}')

                        # Compute the internal magnetic field from the IGRF model for a given
                        # position in GSM coordinate
                        bx_igrf[j, k], by_igrf[j, k], bz_igrf[j, k] = gp.igrf_gsm(x_shu[j, k],
                                                                                  y_shu[j, k],
                                                                                  z_shu[j, k])

                        bx[j, k] = bx_ext[j, k] + bx_igrf[j, k]
                        by[j, k] = by_ext[j, k] + by_igrf[j, k]
                        bz[j, k] = bz_ext[j, k] + bz_igrf[j, k]

                        if (np.sqrt(y_shu[j, k]**2 + z_shu[j, k]**2) > 31):
                            shear[j, k] = np.nan
                        else:
                            shear[j, k] = get_shear([bx[j, k], by[j, k], bz[j, k]], [b_msx[j, k],
                                                   b_msy[j, k], b_msz[j, k]], angle_unit="degrees")

                        rx_en[j, k] = get_rxben([bx[j, k], by[j, k], bz[j, k]], [b_msx[j, k],
                                                 b_msy[j, k], b_msz[j, k]])
                        va_cs[j, k] = get_vcs([bx[j, k], by[j, k], bz[j, k]], [b_msx[j, k],
                                              b_msy[j, k], b_msz[j, k]], n_sh[j, k], 0.1)
                        bisec[j, k] = get_bis([bx[j, k], by[j, k], bz[j, k]],
                                              [b_msx[j, k], b_msy[j, k], b_msz[j, k]])
                        b_sh_ca[j, k] = get_ca([b_msx[j, k], b_msy[j, k], b_msz[j, k]])
                        b_sh_mag[j, k] = np.linalg.norm([b_msx[j, k], b_msy[j, k], b_msz[j, k]])
                        break

    if save_data:

        fn = f'../data/all_data_rx_model_{dr}re_{mp}mp_{model_type}_{today_date}.h5'
        data_file = hf.File(fn, 'w')

        data_file.create_dataset('bx', data=bx)
        data_file.create_dataset('by', data=by)
        data_file.create_dataset('bz', data=bz)

        data_file.create_dataset('bx_ext', data=bx_ext)
        data_file.create_dataset('by_ext', data=by_ext)
        data_file.create_dataset('bz_ext', data=bz_ext)
                                                  
        data_file.create_dataset('bx_igrf', data=bx_igrf)
        data_file.create_dataset('by_igrf', data=by_igrf)
        data_file.create_dataset('bz_igrf', data=bz_igrf)
                                                  
        data_file.create_dataset('b_msx', data=b_msx)
        data_file.create_dataset('b_msy', data=b_msy)
        data_file.create_dataset('b_msz', data=b_msz)
                                                      
        data_file.create_dataset('x_shu', data=x_shu)
        data_file.create_dataset('y_shu', data=y_shu)
        data_file.create_dataset('z_shu', data=z_shu)
                                                      
        data_file.create_dataset('rho_sh', data=rho_sh)
        #data_file.create_dataset('rp', data=rp)
                                                  
        #data_file.create_dataset('r', data=r)
        #data_file.create_dataset('zp', data=zp)
        #data_file.create_dataset('x0', data=x0)
                                                      
        data_file.create_dataset('shear', data=shear)
        data_file.create_dataset('rx_en', data=rx_en)
        data_file.create_dataset('va_cs', data=va_cs)
        data_file.create_dataset('bisec', data=bisec)
        data_file.create_dataset('y_coord', data=y_coord)
        data_file.create_dataset('z_coord', data=z_coord)
        data_file.create_dataset('b_sh_ca', data=b_sh_ca)
        data_file.create_dataset('b_sh_mag', data=b_sh_mag)
        data_file.create_dataset('n_sh', data=n_sh)

        data_file.close()
        print(f'Date saved to file {fn}')

    ridge_finder(image=shear, sigma=2.2, dr=dr, fig_name='shear', c_label='Shear', c_unit=r'${}^\circ$')
    ridge_finder(image=rx_en, sigma=2.8, dr=dr, fig_name='rx-en_nPa_v2', c_label='Reconnection Energy', c_unit='nPa')
    ridge_finder(image=va_cs, sigma=3., dr=dr, fig_name='va-cs', c_label='Exhaust Velocity', c_unit='km/s')
    ridge_finder(image=bisec, sigma=2.2, dr=dr, fig_name='bisec', c_label='Bisection Field', c_unit='nT')
    _ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=b_msy, bz=b_msz, save_fig=True, scale=40,
                           fig_name="magnetosheath")
    _ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=by, bz=bz, save_fig=True, scale=120,
                       fig_name="magnetosphere")


ridge_finder(image=shear, sigma=2.2, dr=dr, fig_name='shear', c_label='Shear', c_unit=r'${}^\circ$')
ridge_finder(image=rx_en, sigma=2.2, dr=dr, fig_name='rx-en', c_label='Reconnection Energy', c_unit='nPa')
ridge_finder(image=va_cs, sigma=2.2, dr=dr, fig_name='va-cs', c_label='Exhaust Velocity', c_unit='km/s')
ridge_finder(image=bisec, sigma=3, dr=dr, fig_name='bisec_sh', c_label='Bisection Field', c_unit='nT')
_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=b_msy, bz=b_msz, save_fig=True, scale=90,
                           fig_name="magnetosheath")
_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=by, bz=bz, save_fig=True, scale=270,
                           fig_name="magnetosphere")

print(f'Took {round(time.time() - start, 3)} seconds')
