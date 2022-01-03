from time import strftime
from mpl_toolkits.axes_grid1.axes_divider import Divider
import numpy as np
import matplotlib.pyplot as plt
import warnings
import geopack as gp
import datetime
import logging
import traceback
import os
import h5py
from skimage.util import pad
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyspedas as spd
import pytplot as ptt

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
    time_obsevation=None,
    dt=30,
    save_data=False,
    data_file="shear_data",
    plot_figure=False,
    save_figure=False,
    figure_file="shear_angle_calculator",
    figure_format="png",
):
    r"""
    Calculate the shear angle between the IMF and the Magnetosheath magnetic field. The code also
    saves the data to a hdf5 file and plots the figure.

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
        Model type to use. Default is "t96". Other option is "t01".

    angle_units : str
        Units of the angle returned by the code. Default is "radians".

    use_real_data : bool
        If set to True, then the code will use the real data from the CDAWEB website. If set to
        False, then the code will use the default values for the parameters.

    time_obsevation : datetime.datetime
        Time of observation. If not given, then the code assumes the time of observation is the
        current time

    """
    if (b_imf is None and use_real_data is False):
        raise ValueError("Interplanetary magnetic field b_imf is not defined")
    if use_real_data:
        if time_obsevation is None:
            raise ValueError("Time of observation must be provided in order to use real time data\
                              Please provide time_obsevation in the format YYYY-MM-DD HH:MM:SS")
        else:
            time_obsevation = datetime.datetime.strptime(time_obsevation, "%Y-%m-%d %H:%M:%S")
            if dt is None:
                dt = 30
                time_range = [(time_obsevation - datetime.timedelta(minutes=dt)).strftime(
                    "%Y-%m-%d %H:%M:%S"), (time_obsevation + datetime.timedelta(minutes=dt)).strftime(
                    "%Y-%m-%d %H:%M:%S")]

        omni_vars = spd.omni.data(trange=trange, level="hro")

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

    if use_real_data is False:
        print("Using default solar wind parameter values. If you want to use real time data, please\
               use the function with the argument 'use_real_data=True'")
        v_imf = np.array([-500, 0.0, 0.0])
        np_imf = 5.0
        sym_h_imf = -30

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

    # Compute the dipole tilt angle
    if use_real_data:
        time_obsevation = datetime.datetime.strptime(time_obsevation, "%Y-%m-%d %H:%M:%S")
        time_obsevation = time_obsevation.replace(tzinfo=datetime.timezone.utc).timestamp()
    else:
        print("Using current time in UTC to compute the dipole tilt angle")
        time_obsevation = datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')
        time_obsevation = time_obsevation.replace(tzinfo=datetime.timezone.utc).timestamp()

    # Compute the dipole tilt angle
    dipole_tilt_angle = gp.recalc(time_obsevation)

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
                    b_msy[j, k] = A * (- b_imf[0] * (y0 / (2 * ll)) + b_imf[2] * (2 - y0**2/(
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
                        shear[j, k] = get_shear([bx[j, k], by[j, k], bz[j, k]],
                                                [b_msx[j, k], b_msy[j, k], b_msz[j, k]],
                                                angle_units=angle_units)
                    break

    if (save_data):
        # Check if the data folder exists, if not then create it.
        if not os.path.exists("data_folder"):
            os.makedirs("data_folder")

        # Save shear data to an hdf5 file
        with h5py.File(f'data_folder/{data_file}_{model_type}_{dr}dr.hdf5', 'w') as f:
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

    if (plot_figure):
        fig, axs = plt.subplots(1, 1, figsize=(10, 10))

        im1 = axs.imshow(abs(shear), origin='lower', cmap=plt.cm.viridis)
        divider = make_axes_locatable(axs)
        cax = divider.append_axes("right", size="5%", pad=0.01)
        cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)
        cbar.ax.tick_params(axis='x', direction='in', top=True, bottom=False, labeltop=True, labelbottom=False, labelsize=18)
        cbar.ax.set_xlabel(f'Shear ({angle_units})', fontsize=18)
        axs.set_xlabel(r'Y [GSM, $R_\oplus$]', fontsize=18)
        axs.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=18)
        plt.show()

    if save_figure:
        fig.savefig(f'{figure_file}_{data_file}_{model_type}_{dr}dr.{figure_format}',
                    bbox_inches='tight', pad_inches=0.05, format=figure_format, dpi=300)
        plt.close()
        print(f"Figure saved to figures/{figure_file}_{data_file}_{model_type}_{dr}dr.png")
