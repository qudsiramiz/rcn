import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pyspedas as spd
import pytplot as ptt
import geopack.geopack as gp
from tabulate import tabulate
import warnings


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
        raise ValueError("trange must be specified as a list of start and end times in the format"
                         "'YYYY-MM-DD HH:MM:SS'.")

    # Check if trange is either a list or an array of length 2
    if not isinstance(trange, (list, np.ndarray)) or len(trange) != 2:
        raise ValueError(
            "trange must be specified as a list or array of length 2 in the format" +
            "'YYYY-MM-DD HH:MM:SS.")

    # Download the OMNI data (default level of 'hro_1min') for the specified timerange.
    omni_varnames = ['BX_GSE', 'BY_GSM', 'BZ_GSM', 'proton_density', 'Vx', 'Vy', 'Vz', 'SYM_H', 'T',
                     'flow_speed']
    omni_vars = spd.omni.data(trange=trange, varnames=omni_varnames, level=omni_level,
                              time_clip=time_clip)

    omni_time = ptt.get_data(omni_vars[0])[0]
    # print(f'omni_time: {omni_time}')
    omni_bx_gse = ptt.get_data(omni_vars[0])[1]
    omni_by_gsm = ptt.get_data(omni_vars[1])[1]
    omni_bz_gsm = ptt.get_data(omni_vars[2])[1]
    omni_np = ptt.get_data(omni_vars[3])[1]
    omni_vx = ptt.get_data(omni_vars[4])[1]
    omni_vy = ptt.get_data(omni_vars[5])[1]
    omni_vz = ptt.get_data(omni_vars[6])[1]
    omni_sym_h = ptt.get_data(omni_vars[7])[1]
    omni_t_p = ptt.get_data(omni_vars[8])[1]
    omni_vsw = ptt.get_data(omni_vars[9])[1]

    # Convert omni_time to datetime objects from unix time
    omni_time_datetime = [datetime.datetime.utcfromtimestamp(t) for t in omni_time]
    # Get trange in datetime format
    omni_trange_time_object = [pd.to_datetime(trange[0]).tz_localize("UTC"),
                               pd.to_datetime(trange[1]).tz_localize("UTC")]

    # Create the dataframe for OMNI data using omni_time_datetime as the index
    omni_df = pd.DataFrame({
        "time": omni_time,
        "bx_gsm": omni_bx_gse,
        "by_gsm": omni_by_gsm,
        "bz_gsm": omni_bz_gsm,
        "vx": omni_vx,
        "vy": omni_vy,
        "vz": omni_vz,
        "np": omni_np,
        "sym_h": omni_sym_h,
        "t_p": omni_t_p,
        "vsw": omni_vsw
    }, index=omni_time_datetime)


    #Get the mean values of the parameters from OMNI data for the time range betwwen
    #omni_trange_time_object[0] and omni_trange_time_object[1]
    time_imf = np.nanmean(omni_df["time"].loc[omni_trange_time_object[0]:
                                              omni_trange_time_object[1]])
    b_imf_x = np.nanmean(omni_df["bx_gsm"].loc[omni_trange_time_object[0]:
                                               omni_trange_time_object[1]])
    b_imf_y = np.nanmean(omni_df["by_gsm"].loc[omni_trange_time_object[0]:
                                               omni_trange_time_object[1]])
    b_imf_z = np.nanmean(omni_df["bz_gsm"].loc[omni_trange_time_object[0]:
                                               omni_trange_time_object[1]])
    vx_imf = np.nanmean(omni_df["vx"].loc[omni_trange_time_object[0]:
                                          omni_trange_time_object[1]])
    vy_imf = np.nanmean(omni_df["vy"].loc[omni_trange_time_object[0]:
                                          omni_trange_time_object[1]])
    vz_imf = np.nanmean(omni_df["vz"].loc[omni_trange_time_object[0]:
                                          omni_trange_time_object[1]])
    np_imf = np.nanmean(omni_df["np"].loc[omni_trange_time_object[0]:
                                          omni_trange_time_object[1]])
    sym_h_imf = np.nanmean(omni_df["sym_h"].loc[omni_trange_time_object[0]:
                                                omni_trange_time_object[1]])
    tp_imf = np.nanmean(omni_df["t_p"].loc[omni_trange_time_object[0]:
                                           omni_trange_time_object[1]])
    vsw_imf = np.nanmean(omni_df["vsw"].loc[omni_trange_time_object[0]:
                                            omni_trange_time_object[1]])

    if (b_imf_z > 15 or b_imf_z < -18):
        warnings.warn(
            f"The given parameters produced the z-component of IMF field (b_imf_z) {b_imf_z} nT,"
            f"which is out of range in which model is valid (-18 nT < b_imf_z < 15 nT)"
        )

    time_imf_hrf = datetime.datetime.utcfromtimestamp(time_imf)

    v_imf = [vx_imf, vy_imf, vz_imf]
    b_imf = [b_imf_x, b_imf_y, b_imf_z]

    imf_clock_angle = np.arctan2(b_imf[1], b_imf[2]) * 180 / np.pi
    if imf_clock_angle < 0:
        imf_clock_angle += 360
    print("IMF parameters found:")
    if (verbose):
        print(tabulate(
            [["Time of observation (UTC)", f"{time_imf_hrf}"],
             ["IMF Magnetic field [GSM] (nT)", f"[{b_imf[0]:.2f}, {b_imf[1]:.2f}, {b_imf[2]:.2f}]"],
             ["IMF Proton density (1/cm^-3)", f"{np_imf:.2f}"],
             ["IMF Plasma velocity (km/sec)", f"[{v_imf[0]:.2f}, {v_imf[1]:.2f}, {v_imf[2]:.2f}]"],
             ["IMF clock angle (degrees)", f"{imf_clock_angle:.2f}"],
             ["IMF Sym H", f"{sym_h_imf:.2f}"]],
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

    rho = np_imf * m_proton * 1.15  # NOTE to self: Unit is fine, do not worry about it
    # print(f"Proton density is {np_imf} 1/cm^3")

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

    # Convert time_imf to datetime object
    time_imf = datetime.datetime.utcfromtimestamp(time_imf)
    # Make a dictionary of all the solar wind parameters
    sw_dict = {}
    sw_dict['time'] = time_imf
    sw_dict['b_imf'] = np.linalg.norm(b_imf)
    sw_dict['b_imf_x'] = b_imf_x
    sw_dict['b_imf_y'] = b_imf_y
    sw_dict['b_imf_z'] = b_imf_z
    sw_dict['rho'] = rho
    sw_dict['ps'] = ps
    sw_dict['p_dyn'] = p_dyn
    sw_dict['sym_h'] = sym_h_imf
    sw_dict['t_p'] = tp_imf
    sw_dict['vsw'] = np.linalg.norm(v_imf)
    sw_dict['imf_clock_angle'] = imf_clock_angle
    #sw_dict['param'] = param

    #return sw_dict, df_fgm, df_fpi, df_mec, df_fgm_fpi
    return sw_dict
