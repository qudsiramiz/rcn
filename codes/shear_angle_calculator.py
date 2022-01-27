import calendar
import datetime
import logging
import os
import traceback
import warnings

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

def get_shear(b_vec_1, b_vec_2, angle_units="radians"):
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

    if (angle_units == "radians"):
        return angle
    elif (angle_units == "degrees"):
        return angle * 180/np.pi
    else:
        raise KeyError("angle_unit must be radians or degrees")


def shear_angle_calculator(
    b_imf=None,
    np_imf=None,
    v_imf=None,
    dmp=0.5,
    dr=0.5,
    model_type="t96",
    angle_units="radians",
    use_real_data=False,
    time_observation=None,
    dt=30,
    save_data=False,
    data_file="shear_data",
    plot_figure=False,
    save_figure=False,
    figure_file="shear_angle_calculator",
    figure_format="png",
    verbose=True,
):
    r"""
    Calculate the shear angle between the IMF and the Magnetosheath magnetic field. The code also
    saves the data to a hdf5 file and plots the figure.
    For this code to work, the following modules must be available/installed:
        - geopack (https://github.com/tsssss/geopack)
        - pyspedas (https://github.com/spedas/pyspedas)

    Parameters
    ----------
    b_imf : array of shape 1x3
        Interplanetary magnetic field vector. If not given, and "use_real_data" is set to False,
        then the code raises a ValueError. In order to use the real data, set "use_real_data" to
        True.

    np_imf : array of shape 1x1
        Interplanetary proton density. If not given, and "use_real_data" is set to False, then the
        code assumes a default value of 5.0 /cm^-3.

    v_imf : array of shape 1x3
        Interplanetary plasma bulk velocity vector. If not given, and "use_real_data" is set to
        False, then the code assumes a default value of [-500, 0, 0] km/sec.

    dmp : float
        Thickness of the magnetopause in earth radii units. Default is 0.5.

    dr : float
        Grid size in earth radii units. Default is 0.5.

    model_type : str
        Model type to use. Default is "t96". Other option is "t01". Needs "geopack" to be installed.

    angle_units : str
        Units of the angle returned by the code. Default is "radians".

    use_real_data : bool
        If set to True, then the code will use the real data from the CDAWEB website. If set to
        False, then the code will use the default values for the parameters. For this to work, you
        must have the following modules installed:
        - pyspedas (https://github.com/spedas/pyspedas)

    time_observation : str
        Time of observation. If not given, then the code assumes the time of observation is the
        current time. The time of observation must be in the format "YYYY-MM-DDTHH:MM:SS".

    dt : float
        Time duration, in minutes, for which the IMF data is to be considered. Default is 30
        mcentered around the time of observation.

    save_data : bool
        If set to True, then the data will be saved to a hdf5 file. Default is False. Needs the h5py
        package to be installed.

    data_file : str
        Name of the hdf5 file to which data is to be saved. Default is "shear_data".

    plot_figure : bool
        If set to True, then the figure will be plotted and shown. Default is False.

    save_figure : bool
        If set to True, then the figure will be saved. Default is False.

    figure_file : str
        Name of the figure file. Default is "shear_angle_calculator".

    figure_format : str
        Format of the figure file. Default is "png". Other option can be "pdf".

    verbose : bool If set to True, then the code will print out the progress of the code at several
        points. Default is True. For this to work, you must have the following modules installed:
         - tabulate (https://pypi.org/project/tabulate/)
    """
    # Check if the required modules are installed
    try:
        import geopack.geopack as gp
    except ImportError:
        raise ImportError("geopack is not installed. Please install it using the command: pip install geopack or directly from the source at GitHub (https://github.com/tsssss/geopack).")

    if use_real_data:
        try:
            import pyspedas as spd
            import pytplot as ptt
        except ImportError:
            raise ImportError("pyspedas is not installed. Please install it using the command: pip install pyspedas")

    if save_data:
        try:
            import h5py as hf
        except ImportError:
            raise ImportError("h5py is not installed. Please install it using the command: pip install h5py")

    if verbose:
        try:
            from tabulate import tabulate
        except ImportError:
            raise ImportError("The required module 'tabulate' is not installed. Please install it using the command pip install tabulate.")

    if (b_imf is None and use_real_data is False):
        raise ValueError("Interplanetary magnetic field b_imf is not defined. If you do not wish to provide IMF magnetic field, set use_real_data to True")
    if use_real_data:
        if verbose:
            print("Attempting to download real data from the CDSAWEB website using PysSpedas \n")
        if time_observation is None:
            time_observation = "2020-02-16 17:00:00"
            time_observation = datetime.datetime.strptime(time_observation, "%Y-%m-%d %H:%M:%S")
            if verbose:
                print(f"Time of observation is not given. Defaulting to: {time_observation} UTC \n")
            if dt is None:
                print("dt, observation time range, is not defined. Setting dt to 30 minutes \n")
                dt = 30
            time_range = [(time_observation - datetime.timedelta(minutes=dt)).strftime(
                "%Y-%m-%d %H:%M:%S"), (time_observation + datetime.timedelta(minutes=dt)).strftime(
                "%Y-%m-%d %H:%M:%S")]
            if verbose:
                print(f"Downloading data from {time_range[0]} to {time_range[1]} \n")

        else:
            time_observation = datetime.datetime.strptime(time_observation, "%Y-%m-%d %H:%M:%S")
            if dt is None:
                print("dt, observation time range, is not defined. Setting dt to 30 minutes \n")
                dt = 30
            time_range = [(time_observation - datetime.timedelta(minutes=dt)).strftime(
                "%Y-%m-%d %H:%M:%S"), (time_observation + datetime.timedelta(minutes=dt)).strftime(
                "%Y-%m-%d %H:%M:%S")]
            if verbose:
                print(f"Downloading data from {time_range[0]} to {time_range[1]} \n")

        spd.omni.data(trange=time_range, level="hro", time_clip=False)

        omni_time = ptt.get_data('BX_GSE')[0]
        # omni_time_unix = np.vectorize(datetime.utcfromtimestamp)(omni_time[:]) # converting omni_time
        # from unixtime to utc datetime object array in python

        omni_bx_gse = ptt.get_data('BX_GSE')[1]  # NOTE: x-direction of GSE and GSM are same
        omni_by_gsm = ptt.get_data('BY_GSM')[1]
        omni_bz_gsm = ptt.get_data('BZ_GSM')[1]
        omni_np = ptt.get_data('proton_density')[1]
        omni_vx = ptt.get_data('Vx')[1]
        omni_vy = ptt.get_data('Vy')[1]
        omni_vz = ptt.get_data('Vz')[1]
        omni_sym_h = ptt.get_data('SYM_H')[1]

        if verbose:
            print("Data downloaded, processing data to find IMF parameters \n")
        time_imf = np.nanmedian(omni_time)
        b_imf = np.array([np.nanmedian(omni_bx_gse),
                          np.nanmedian(omni_by_gsm),
                          np.nanmedian(omni_bz_gsm)])

        if (b_imf[2] > 15 or b_imf[2] < -18):
            warnings.warn(
                f"The given parameters produced the z-component of IMF field (b_imf_z) {b_imf[2]} nT,"
                f"which is out of range in which model is valid (-18 nT < b_imf_z < 15 nT)"
            )

        np_imf = np.nanmedian(omni_np)
        v_imf = np.array([np.nanmedian(omni_vx), np.nanmedian(omni_vy), np.nanmedian(omni_vz)])
        sym_h_imf = np.nanmedian(omni_sym_h)
        clock_angle = np.arctan2(b_imf[1], b_imf[2]) * 180 / np.pi

        if verbose:
            print("Computed IMF parameters:")
            print(tabulate(
                [["Time of observation (UTC)", time_observation],
                 ["IMF Magnetic field (nT)", np.round(b_imf, 3)],
                 ["IMF Proton density (1/cm^-3)", np.round(np_imf, 3)],
                 ["IMF Plasma velocity (km/sec)", np.round(v_imf, 3)],
                 ["IMF Sym H", np.round(sym_h_imf, 3)]],
                headers=["Parameter", "Value"], tablefmt="fancy_grid", floatfmt=".2f",
                numalign="center"))

    if use_real_data is False:
        if verbose:
            print("Using default solar wind parameter values. If you want to use real time " +
                  "data, please use the function with the argument 'use_real_data=True' \n")
        v_imf = np.array([-500, 0.0, 0.0])
        np_imf = 5.0
        sym_h_imf = -30
        clock_angle = np.arctan2(b_imf[1], b_imf[2]) * 180 / np.pi

    m_p = 1.672e-27  # Mass of proton in SI unit

    rho = np_imf * m_p * 1.15  # 1.15 instead of 1.0 to account for the alpha particles

    #  Solar wind ram pressure in nPa, including roughly 4% Helium++ contribution
    p_dyn = 1.6726e-6 * 1.15 * np_imf * np.linalg.norm(v_imf) ** 2

    if (p_dyn > 8.5 or p_dyn < 0.5):
        warnings.warn(
            f"The given parameters produced a dynamic pressure of {p_dyn} nPa which is out of"
            f" range in which model is valid (0.5 nPa < p_dyn < 8.5 nPa)",
        )
    param = [p_dyn, sym_h_imf, b_imf[1], b_imf[2], 0, 0, 0, 0, 0, 0]

    if verbose:
        print("Input parameters for the model:")
        print(tabulate(
            [["Solar wind dynamic pressure (nPa)", np.round(p_dyn, 3)],
             ["IMF Sym H", np.round(sym_h_imf, 3)],
             ["B_IMF_Y (nT)", np.round(b_imf[1], 3)],
             ["B_IMF_Z (nT)", np.round(b_imf[2], 3)],
             ["Clock Angle (degrees)", np.round(clock_angle, 3)]],
            headers=["Parameter", "Value"], tablefmt="fancy_grid", floatfmt=".2f",
            numalign="center"))

    # Compute the dipole tilt angle
    if use_real_data:
        time_dipole = calendar.timegm(time_observation.utctimetuple())
    else:
        time_dipole = calendar.timegm(datetime.datetime.utcnow().utctimetuple())
        #time_observation = datetime.datetime.strptime(time_observation, "%Y-%m-%d %H:%M:%S")
        #time_dipole = calendar.timegm(time_observation.utctimetuple())

    # Compute the dipole tilt angle
    dipole_tilt_angle = gp.recalc(time_dipole)

    if verbose:
        if angle_units == "radians":
            print(tabulate([["Dipole tilt angle units", angle_units],
                            ["Dipole tilt angle", np.round(dipole_tilt_angle, 3)]],
                           headers=["Parameter", "Value"], tablefmt="fancy_grid", floatfmt=".3f",
                           numalign="center"))
        elif angle_units == "degrees":
            print(tabulate([["Dipole tilt angle units", angle_units],
                            ["Dipole tilt angle", np.round(dipole_tilt_angle * 180 / np.pi, 3)]],
                           headers=["Parameter", "Value"], tablefmt="fancy_grid", floatfmt=".3f",
                           numalign="center"))

    if verbose:
        print("Computing Earth's magnetic field \n")

    if dr is None:
        dr = 0.5
    if dmp is None:
        dmp = 0.5        
    n_arr = int(30 / dr) + 1

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

    shear = np.full((n_arr, n_arr), np.nan)
    y_coord = np.full((n_arr, n_arr), np.nan)
    z_coord = np.full((n_arr, n_arr), np.nan)
    n_sh = np.full((n_arr, n_arr), np.nan)

    d_theta = np.pi / 100

    # Shue et al.,1998, equation 10
    ro = (10.22 + 1.29 * np.tanh(0.184 * (b_imf[2] + 8.14))) * (p_dyn)**(-1.0 / 6.6)

    # Shue et al.,1998, equation 11
    # alpha = (0.58 - 0.010 * b_imf_z) * (1 + 0.010 * p_dyn)
    alpha = (0.58 - 0.007 * b_imf[2]) * (1 + 0.024 * np.log(p_dyn))
    # Stand off position of the magnetopause
    rmp = ro * (2 / (1 + np.cos(0.0))) ** alpha

    A = 2
    len_y = int(30 / dr) + 1
    len_z = int(30 / dr) + 1
    count = 0

    for j in range(0, len_y):
        y0 = 15 - int(j * dr)
        for k in range(0, len_z):
            z0 = 15 - int(k * dr)
            rp = np.sqrt(y0**2 + z0**2)  # Projection of r into yz-plane

            for index in range(0, 100):

                theta = index * d_theta
                r = ro * (2 / (1 + np.cos(theta))) ** alpha
                # not really in z direction, but a distance in yz plane
                zp = r * np.sin(theta)
                x0 = r * np.cos(theta)

                #if x0 == 0:
                #    signx = 1.0
                #else:
                #    signx = np.sign(x0)

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
                    x_shu[j, k] = (r - dmp) * np.cos(theta)
                    phi = np.arctan2(z0, y0)
                    # print( j, k, theta, x_shu[j,k])

                    if (abs(y0) == 0 or abs(z0) == 0):
                        if(abs(y0) == 0):
                            y_shu[j, k] = 0
                            z_shu[j, k] = (r - dmp) * np.sin(theta)
                        elif (abs(z0) == 0):
                            z_shu[j, k] = 0
                            y_shu[j, k] = (r - dmp) * np.sin(theta)
                    else:
                        z_shu[j, k] = np.sqrt(
                            (rp - 1.0)**2 / (1 + np.tan(phi)**(-2)))
                        y_shu[j, k] = z_shu[j, k] / np.tan(phi)

                    rho_sh[j, k] = rho * (1.509 * np.exp(x_shu[j, k] / rmp) + .1285)
                    n_sh[j, k] = rho_sh[j, k] / m_p
                    y_shu[j, k] = abs(y_shu[j, k]) * signy
                    z_shu[j, k] = abs(z_shu[j, k]) * signz

                    # Cooling JGR 2001 Model, equation 9 to 12
                    # the distance from the focus to the magnetopause surface
                    ll = 3 * rmp/2 - x0
                    b_msx[j, k] = - A * (- b_imf[0] * (1 - rmp / (2 * ll)) + b_imf[1] * (y0 / ll)
                                    + b_imf[2] * (z0 / ll))
                    b_msy[j, k] = A * (- b_imf[0] * (y0 / (2 * ll)) + b_imf[1] * (2 - y0**2/(
                                  ll * rmp)) - b_imf[2]* (y0 * z0 / (ll * rmp)))
                    b_msz[j, k] = A * (- b_imf[0] * (z0 / (2 * ll)) - b_imf[1] * (y0 * z0 / (ll
                                  * rmp)) + b_imf[2] * (2 - z0**2 / (ll * rmp)))

                    # Compute the external magnetic field from the T95 model for a given
                    # position and time, in GSM coordinate
                    try:
                        if(model_type == 't96'):
                            bx_ext[j, k], by_ext[j, k], bz_ext[j, k] = gp.t96.t96(
                                param,
                                dipole_tilt_angle,
                                x_shu[j, k],
                                y_shu[j, k],
                                z_shu[j, k]
                            )
                        elif(model_type == 't01'):
                            bx_ext[j, k], by_ext[j, k], bz_ext[j, k] = gp.t01.t01(
                                param,
                                dipole_tilt_angle,
                                x_shu[j, k],
                                y_shu[j, k],
                                z_shu[j, k]
                            )
                        else:
                            raise ValueError(
                                "Model type must be set to 't96' or 't01'.")
                    except Exception:
                        logging.error(traceback.format_exc())
                        if((y_shu[j, k] < 15 and y_shu[j, k] > -15) or
                           (z_shu[j, k] < 15 and z_shu[j, k] > -15)):
                            if verbose:
                                print(f'Skipped for {x_shu[j, k], y_shu[j, k], z_shu[j, k]}')
                            count += 1
                        else:
                            if verbose:
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
                        shear[j, k] = get_shear([bx[j, k], by[j, k], bz[j, k]],
                                                [b_msx[j, k], b_msy[j, k], b_msz[j, k]],
                                                angle_units=angle_units)
                    break
    if verbose:
        print("Earth's magnetic field computed at the location of the magnetopause and shear"+\
            " calculated \n")

    if (save_data):
        # Check if the data folder exists, if not then create it.
        if not os.path.exists("data_folder"):
            os.makedirs("data_folder")
            if verbose:
                print("Created data_folder folder")

        # Save shear data to an hdf5 file
        with hf.File(f'data_folder/{data_file}_{model_type}_{dr}dr.hdf5', 'w') as f:
            f.create_dataset('x_shu', data=x_shu)
            f.create_dataset('y_shu', data=y_shu)
            f.create_dataset('z_shu', data=z_shu)
            f.create_dataset('bx_ext', data=bx_ext)
            f.create_dataset('by_ext', data=by_ext)
            f.create_dataset('bz_ext', data=bz_ext)
            f.create_dataset('bx_igrf', data=bx_igrf)
            f.create_dataset('by_igrf', data=by_igrf)
            f.create_dataset('bz_igrf', data=bz_igrf)
            f.create_dataset('bx', data=bx)
            f.create_dataset('by', data=by)
            f.create_dataset('bz', data=bz)
            f.create_dataset('shear', data=shear)
            f.create_dataset('n_sh', data=n_sh)
            f.create_dataset('rho_sh', data=rho_sh)
            f.create_dataset('b_msx', data=b_msx)
            f.create_dataset('b_msy', data=b_msy)
            f.create_dataset('b_msz', data=b_msz)
            f.close()
        print(f"Data saved to {data_file}_{model_type}_{dr}dr.hdf5")

    fig, axs = plt.subplots(1, 1, figsize=(8, 6))

    # NOTE: This is a hack to ensure that the output of shear angle, reconnection energy etc. agrees
    # with what has been reported in literature. Plot for shear angle seems to agree reasonably well
    # (for "trange = ['2016-12-07 05:11:00', '2016-12-07 05:21:00']") with the one reported by
    # FuselierJGR2019 (doi:10.1029/2019JA027143, see fig. 4).    
    image_rotated = np.transpose(np.flipud(np.fliplr(shear)))
    # Smoothen the image
    image_smooth = sp.ndimage.filters.gaussian_filter(image_rotated, sigma=[5, 5], mode='nearest')
    im = axs.imshow(np.transpose(np.flipud(np.fliplr(image_smooth))), extent=[-15, 15, -15, 15],
    origin='lower', cmap=plt.cm.viridis)
    divider = make_axes_locatable(axs)

    patch = patches.Circle((0, 0), radius=15, transform=axs.transData, fc='none', ec='k', lw=0.1)
    axs.add_patch(patch)
    im.set_clip_path(patch)

    axs.tick_params(axis="both", direction="in", top=True, labeltop=False, bottom=True,
                    labelbottom=True, left=True, labelleft=True, right=True, labelright=False, labelsize=14)
    axs.set_xlabel(r'Y [GSM, $R_\oplus$]', fontsize=18)
    axs.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=18)

    cax = divider.append_axes("top", size="5%", pad=0.01)
    cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=None, fraction=0.05,
                     pad=0.01)
    cbar.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=True,
                        labelbottom=False, pad=0.01, labelsize=14)
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.set_xlabel(f'Shear Angle ({angle_units})', fontsize=18)

    # Write the timme range on the plot

    clock_angle_degrees = np.round(clock_angle, 2)
    dipole_angle_degrees = np.round(dipole_tilt_angle * 180 / np.pi, 2)

    axs.text(1.0, 0.5, f'Time of observation: {time_observation}', horizontalalignment='left',
              verticalalignment='center', transform=axs.transAxes, rotation=270, color='r')

    axs.text(0.01, 0.99, f'Clock Angle: {clock_angle_degrees}$^\circ$', horizontalalignment='left',
              verticalalignment='top', transform=axs.transAxes, rotation=0, color='r')

    axs.text(0.99, 0.99, f'Dipole tilt: {dipole_angle_degrees}$^\circ$',
    horizontalalignment='right', verticalalignment='top', transform=axs.transAxes, rotation=0,
    color='r')

    if (plot_figure):
        plt.show()

    if (save_figure):
        fig.savefig(f'{figure_file}_{model_type}_{dr}dr.{figure_format}',
                    bbox_inches='tight', pad_inches=0.05, format=figure_format, dpi=300)
        plt.close()
        print(f"Figure saved to {figure_file}_{data_file}_{model_type}_{dr}dr.png")

    return shear

inputs = {
    "b_imf" : None,
    "np_imf" : None,
    "v_imf" : None,
    "dmp" : None,
    "dr" : None,
    "model_type" : "t96",
    "angle_units" : "degrees",
    "use_real_data" : True,
    "time_observation" : "2016-12-07 05:16:00",
    "dt" : 5,
    "save_data" : False,
    "data_file" : None,
    "plot_figure" : True,
    "save_figure" : True,
    "figure_file" : "shear_angle_calculator",
    "figure_format" : "pdf",
    "verbose" : True
}

shear_angle_calculator(**inputs)
