import datetime
import os

import matplotlib.pyplot as plt
import more_itertools as mit

import numpy as np
import pandas as pd
import pyspedas as spd
import pyspedas.mms.cotrans.mms_cotrans_lmn as mms_cotrans_lmn
import pytplot as ptt
import pytz


def jet_reversal_check(crossing_time=None, dt=90, probe=3, data_rate='fast', level='l2',
                       coord_type='lmn', data_type='dis-moms', time_clip=True, latest_version=True,
                       jet_len=5, figname='mms_jet_reversal_check', date_obs=None,
                       fname='../data/mms_jet_reversal_times.csv',
                       error_file_log_name="../data/mms_jet_reversal_check_error_log.csv",
                       verbose=True
                       ):
    """
    For a given crossing time and a given probe, the function finds out if MMS observed a jet
    during magnetopause crossing. If there was indeed a jet reversal, then the function saves that
    time to a csv file, along with the probe number, position of the spacecraft, and the time of
    the crossing.

    Parameters
    ----------
    crossing_time : datetime.datetime
        The time of the crossing of the magnetopause.
    dt : int
        The time interval to look for a jet reversal.
    probe : int
        The probe number. Can be 1, 2, 3, or 4. Default is 3. (since the magnetopause crossing
        times are for MMS3)
    data_rate : str
        The data rate. Can be 'fast' or 'srvy' or 'brst. Default is 'fast'.
    level : str
        The data level. Can be either 'l1' or 'l2'. Default is 'l2'.
    coord_type : str
        The coordinate type. Can be either 'lmn' or 'gse'. Default is 'lmn'.
    data_type : str
        The data type. Can be either 'dis-moms' or 'des-moms'. Default is 'dis-moms'.
    time_clip : bool
        Whether or not to clip the data to the time range of the crossing time. Default is True.
    latest_version : bool
        Whether or not to use the latest version of the data. Default is True.
    jet_len : int
        The time length for which the jet should be observed. Default is 5.
    figname : str
        The name of the figure to be saved. Default is 'mms_jet_reversal_check'.
    date_obs : str
        The date of the observation. Default is None.
    fname : str
        The name of the csv file to save the data to. Default is 'mms_jet_reversal_times.csv'.
    error_file_log_name : str
        The name of the csv file to save the error log to. Default is
        '../data/mms_jet_reversal_check_error_log.csv'.
    verbose : bool
        Whether or not to print the status of the function. Default is True.

    Returns
    -------
    None
    """

    # Define the absolute permeability of free space in m^2 kg^-1 s^-1
    mu_0 = 4 * np.pi * 1e-7

    # Define the Boltzmann constant in J K^-1
    k_B = 1.38064852e-23

    # Define the conversion factor for converting temperature from ev to K
    ev_to_K = 11604.525

    crossing_time_min = crossing_time - datetime.timedelta(seconds=dt)
    crossing_time_max = crossing_time + datetime.timedelta(seconds=dt)
    trange = [crossing_time_min, crossing_time_max]

    # Get the index corresponding to the crossing time in the data
    # try:
    #     df_crossing_temp = pd.read_csv("../data/mms_jet_reversal_times_list_20221017_beta.csv",
    #                                    index_col=False)
    #     # Change the formatting of the date
    #     for xx in range(len(df_crossing_temp)):
    #         df_crossing_temp['jet_time'][xx] = df_crossing_temp['jet_time'][xx].split(
    #                                          '+')[0].split('.')[0]
    #     df_crossing_temp.set_index("jet_time", inplace=True)
    #     crossing_time_str = crossing_time.strftime("%Y-%m-%d %H:%M:%S")[:]
    # except:
    # df_crossing_temp = pd.read_csv("../data/mms_magnetopause_crossings.csv", index_col=False)
    # df_crossing_temp.set_index("DateStart", inplace=True)
    # crossing_time_str = crossing_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
    df_crossing_temp = pd.read_csv("../data/event_list_MMS_jet_reversals_from_steve.csv",
                                   index_col=False)
    df_crossing_temp.set_index("jet_time", inplace=True)
    crossing_time_str = crossing_time.strftime("%Y-%m-%d %H:%M:%S")
    ind_crossing = np.where(df_crossing_temp.index == crossing_time_str)[0][0]

    # Get the data from the FPI
    mms_fpi_varnames = [f'mms{probe}_dis_numberdensity_{data_rate}',
                        f'mms{probe}_dis_bulkv_gse_{data_rate}',
                        f'mms{probe}_dis_temppara_{data_rate}',
                        f'mms{probe}_dis_tempperp_{data_rate}',
                        f'mms{probe}_dis_energyspectr_omni_{data_rate}',
                        f'mms{probe}_des_energyspectr_omni_{data_rate}']

    _ = spd.mms.fpi(trange=trange, probe=probe, data_rate=data_rate, level=level,
                    datatype=data_type, varnames=mms_fpi_varnames, time_clip=time_clip,
                    latest_version=latest_version)

    mms_fpi_time_unix = ptt.get_data(mms_fpi_varnames[0])[0]
    # Convert the time to a datetime object
    mms_fpi_time_local = pd.to_datetime(mms_fpi_time_unix, unit='s')
    mms_fpi_time = mms_fpi_time_local.tz_localize(pytz.utc)

    mms_fpi_numberdensity = ptt.get_data(mms_fpi_varnames[0])[1]
    _ = ptt.get_data(mms_fpi_varnames[1])[1:4][0]
    mms_fpi_temppara = ptt.get_data(mms_fpi_varnames[2])[1]
    mms_fpi_tempperp = ptt.get_data(mms_fpi_varnames[3])[1]

    # Store both the temperatures in the ptt
    ptt.store_data('Tp', data=[f'mms{probe}_dis_temppara_{data_rate}',
                               f'mms{probe}_dis_tempperp_{data_rate}'])

    # Covert gse to gsm
    _ = spd.cotrans(name_in=f'mms{probe}_dis_bulkv_gse_{data_rate}',
                    name_out=f'mms{probe}_dis_bulkv_gsm_{data_rate}', coord_in='gse',
                    coord_out='gsm')

    mms_fpi_bulkv_gsm = ptt.get_data(f'mms{probe}_dis_bulkv_gsm_{data_rate}')[1:4][0]
    mms_fpi_bulkv_gse = ptt.get_data(f'mms{probe}_dis_bulkv_gse_{data_rate}')[1:4][0]

    if coord_type == 'lmn':
        # Convert gsm to lmn
        _ = mms_cotrans_lmn.mms_cotrans_lmn(name_in=f'mms{probe}_dis_bulkv_gsm_{data_rate}',
                                            name_out=f'mms{probe}_dis_bulkv_lmn_{data_rate}',
                                            gse=False, gsm=True, probe=str(probe), data_rate=data_rate)

        mms_fpi_bulkv_lmn = ptt.get_data(f'mms{probe}_dis_bulkv_lmn_{data_rate}')[1:4][0]

    # Create a dataframe with the FPI data
    if coord_type == 'lmn':
        df_mms_fpi = pd.DataFrame(index=mms_fpi_time, data={'np': mms_fpi_numberdensity,
                                                            'vp_lmn_l': mms_fpi_bulkv_lmn[:, 0],
                                                            'vp_lmn_m': mms_fpi_bulkv_lmn[:, 1],
                                                            'vp_lmn_n': mms_fpi_bulkv_lmn[:, 2],
                                                            'tp_para': mms_fpi_temppara,
                                                            'tp_perp': mms_fpi_tempperp})
    else:
        df_mms_fpi = pd.DataFrame(index=mms_fpi_time, data={'np': mms_fpi_numberdensity,
                                                            'vp_gsm_x': mms_fpi_bulkv_gsm[:, 0],
                                                            'vp_gsm_y': mms_fpi_bulkv_gsm[:, 1],
                                                            'vp_gsm_z': mms_fpi_bulkv_gsm[:, 2],
                                                            'vp_gse_x': mms_fpi_bulkv_gse[:, 0],
                                                            'vp_gse_y': mms_fpi_bulkv_gse[:, 1],
                                                            'vp_gse_z': mms_fpi_bulkv_gse[:, 2],
                                                            'tp_para': mms_fpi_temppara,
                                                            'tp_perp': mms_fpi_tempperp})

    # Make sure that the time indices are in increasing order
    df_mms_fpi = df_mms_fpi.sort_index()

    if verbose:
        print("\n\033[1;32m FPI dataframe created \033[0m \n")
        print(f"The fpi Datafram:\n {df_mms_fpi.head()}")

    # Add rolling median to the dataframe
    df_mms_fpi['np_rolling_median'] = df_mms_fpi['np'].rolling('60s', center=True).median()

    if coord_type == 'lmn':
        df_mms_fpi['vp_lmn_l_rolling_median'] = df_mms_fpi['vp_lmn_l'].rolling(
                                                                        '60s', center=True).median()
        df_mms_fpi['vp_lmn_m_rolling_median'] = df_mms_fpi['vp_lmn_m'].rolling(
                                                                        '60s', center=True).median()
        df_mms_fpi['vp_lmn_n_rolling_median'] = df_mms_fpi['vp_lmn_n'].rolling(
                                                                        '60s', center=True).median()

    else:
        df_mms_fpi['vp_gsm_x_rolling_median'] = df_mms_fpi['vp_gsm_x'].rolling(
                                                                        '60s', center=True).median()
        df_mms_fpi['vp_gsm_y_rolling_median'] = df_mms_fpi['vp_gsm_y'].rolling(
                                                                        '60s', center=True).median()
        df_mms_fpi['vp_gsm_z_rolling_median'] = df_mms_fpi['vp_gsm_z'].rolling(
                                                                        '60s', center=True).median()

    df_mms_fpi['tp_para_rolling_median'] = df_mms_fpi['tp_para'].rolling(
                                                                        '60s', center=True).median()
    df_mms_fpi['tp_perp_rolling_median'] = df_mms_fpi['tp_perp'].rolling(
                                                                        '60s', center=True).median()

    if coord_type == 'lmn':
        df_mms_fpi['vp_diff_x'] = df_mms_fpi['vp_lmn_l'] - df_mms_fpi['vp_lmn_l_rolling_median']
        df_mms_fpi['vp_diff_y'] = df_mms_fpi['vp_lmn_m'] - df_mms_fpi['vp_lmn_m_rolling_median']
        df_mms_fpi['vp_diff_z'] = df_mms_fpi['vp_lmn_n'] - df_mms_fpi['vp_lmn_n_rolling_median']
        # df_mms_fpi['vp_diff_z'] = df_mms_fpi['vp_lmn_l'] - np.nanmean(df_mms_fpi['vp_lmn_l'])
    else:
        df_mms_fpi['vp_diff_x'] = df_mms_fpi['vp_gsm_x'] - df_mms_fpi['vp_gsm_x_rolling_median']
        df_mms_fpi['vp_diff_y'] = df_mms_fpi['vp_gsm_y'] - df_mms_fpi['vp_gsm_y_rolling_median']
        df_mms_fpi['vp_diff_z'] = df_mms_fpi['vp_gsm_z'] - df_mms_fpi['vp_gsm_z_rolling_median']

    # Get the data from the FGM
    # mms_fgm_varnames = [f'mms{probe}_fgm_b_gsm_srvy_l2_bvec']
    # _ = spd.mms.fgm(trange=trange, probe=probe, time_clip=time_clip, latest_version=True,
    #                 varnames=[f"mms{probe}_fgm_b_gsm_srvy_{level}",
    #                           f"mms{probe}_fgm_r_gsm_srvy_{level}"], get_fgm_ephemeris=True)
    _ = spd.mms.fgm(trange=trange, probe=probe, time_clip=time_clip, latest_version=True,
                    get_fgm_ephemeris=True)
    # Get the time corresponding to the FGM data
    mms_fgm_time = ptt.get_data(f"mms{probe}_fgm_b_gsm_srvy_{level}")[0]
    # Convert the time to a datetime object
    mms_fgm_time = pd.to_datetime(mms_fgm_time, unit='s')
    mms_fgm_time = mms_fgm_time.tz_localize(pytz.utc)

    mms_fgm_b_gsm = ptt.get_data(f'mms{probe}_fgm_b_gsm_srvy_{level}')[1:4][0]
    mms_fgm_b_gse = ptt.get_data(f'mms{probe}_fgm_b_gse_srvy_{level}')[1:4][0]
    mms_fgm_r_gsm = ptt.get_data(f'mms{probe}_fgm_r_gsm_srvy_{level}')[1:4][0]

    if coord_type == 'lmn':
        _ = mms_cotrans_lmn.mms_cotrans_lmn(name_in=f'mms{probe}_fgm_b_gsm_srvy_{level}_bvec',
                                            name_out=f'mms{probe}_fgm_b_lmn_srvy_{level}',
                                            gsm=True, probe=str(probe), data_rate=data_rate)

        mms_fgm_b_lmn = ptt.get_data(f'mms{probe}_fgm_b_lmn_srvy_{level}')[1:4][0]

    # Create a dataframe with the FGM data
    if coord_type == 'lmn':
        df_mms_fgm = pd.DataFrame(index=mms_fgm_time, data={'b_lmn_l': mms_fgm_b_lmn[:, 0],
                                                            'b_lmn_m': mms_fgm_b_lmn[:, 1],
                                                            'b_lmn_n': mms_fgm_b_lmn[:, 2]})
    else:
        df_mms_fgm = pd.DataFrame(index=mms_fgm_time, data={'b_gsm_x': mms_fgm_b_gsm[:, 0],
                                                            'b_gsm_y': mms_fgm_b_gsm[:, 1],
                                                            'b_gsm_z': mms_fgm_b_gsm[:, 2],
                                                            'b_gse_x': mms_fgm_b_gse[:, 0],
                                                            'b_gse_y': mms_fgm_b_gse[:, 1],
                                                            'b_gse_z': mms_fgm_b_gse[:, 2]})

    # Make sure all time indices are in increasing order
    df_mms_fgm = df_mms_fgm.sort_index()

    if verbose:
        print("\n\033[1;32m FGM dataframe created \033[0m \n")
        print(f"FGM Dataframe: \n{df_mms_fgm.head()} \n")
    # Merge the two dataframes
    df_mms = pd.merge_asof(df_mms_fpi, df_mms_fgm, left_index=True, right_index=True)

    # Compute the median time difference between the points
    # NOTE: The reason for factor 10**9 is to convert the time difference to seconds from
    # nanoseconds which is the format in which df_mms.index.astype(np.int64) outputs time
    time_cadence_median = np.median(np.diff(df_mms.index.astype(np.int64) / 10**9))

    # Make the time index timezone aware
    try:
        df_mms.index = df_mms.index.tz_localize(pytz.utc)
    except Exception:
        if verbose:
            print("\033[1;31m Timezone conversion failed \033[0m \n")

    # Check if magnetosheath and magnetopshere are present and if they are then get the indices
    # corresponding to the magnetosheath and magnetosphere
    ind_range_msp, ind_range_msh = check_msp_msh_location(df_mms=df_mms,
                                                          time_cadence_median=time_cadence_median,
                                                          verbose=verbose)

    (jet_detection, delta_v_min, delta_v_max, t_jet_center, ind_jet_center,
     ind_jet_center_minus_1_min, ind_jet_center_plus_1_min, vp_lmn_diff_l) = check_jet_location(
        df_mms=df_mms, jet_len=jet_len, v_thresh=70, ind_msh=ind_range_msh,
        time_cadence_median=time_cadence_median, verbose=verbose, ind_crossing=ind_crossing,
        date_obs=date_obs
    )

    if jet_detection:
        t_delta_v_min = df_mms.index[ind_jet_center_minus_1_min:ind_jet_center]
        t_delta_v_max = df_mms.index[ind_jet_center:ind_jet_center_plus_1_min]

        # Convert t_delta_v_min and t_delta_v_max to from datetime object to UNIX time
        t_delta_v_min_unix = [t_delta_v_min[i].timestamp() for i in range(len(t_delta_v_min))]
        t_delta_v_max_unix = [t_delta_v_max[i].timestamp() for i in range(len(t_delta_v_max))]
        t_delta_v_min_unix = np.array(t_delta_v_min_unix)
        t_delta_v_max_unix = np.array(t_delta_v_max_unix)

        # Add delta_v_min and delta_v_max to ptt
        ptt.store_data('delta_v_min', data={'x': t_delta_v_min_unix,
                                            'y': delta_v_min})
        ptt.store_data('delta_v_max', data={'x': t_delta_v_max_unix,
                                            'y': delta_v_max})
        ptt.store_data('delta_v', data=['delta_v_min', 'delta_v_max'])

        # Add vp_lmn_diff_l to ptt
        ptt.store_data('vp_lmn_diff_l', data={'x': mms_fpi_time_unix,
                                              'y': vp_lmn_diff_l})

        # Add delta_v and vp_lmn_diff_l to ptt
        ptt.store_data('delta_v_vp_lmn_diff_l', data=['delta_v_min', 'delta_v_max', 'vp_lmn_diff_l'])

        # If jet was detected, then check if walen relation is satisfied

        (walen_relation_satisfied_v1, walen_relation_satisfied_v2, theta_w_deg, R_w,
         theta_all_deg, R_all, ind_walen_vals, ind_walen_check) = check_walen_relation(
                                                                 df_mms=df_mms,
                                                                 t_jet_center=t_jet_center,
                                                                 coord_type=coord_type,
                                                                 ind_msh=ind_range_msh,
                                                                 n_points_walen=jet_len,
                                                                 jet_detection=jet_detection
                                                                 )

        # Select the times corresponding to walen check
        t_walen_check = df_mms.index[ind_walen_check]

        # Convert t_walen_check to from datetime object to UNIX time
        t_walen_check_unix = [t_walen_check[i].timestamp() for i in range(len(t_walen_check))]

        # Add theta_w_deg and R_w to ptt
        ptt.store_data('theta_w_deg', data={'x': t_walen_check_unix, 'y': theta_w_deg})
        ptt.store_data('R_w', data={'x': t_walen_check_unix, 'y': R_w})
        ptt.store_data('theta_R_w', data=['theta_w_deg', 'R_w'])
    else:
        walen_relation_satisfied_v1 = False
        walen_relation_satisfied_v2 = False
        theta_w_deg = np.nan
        R_w = np.nan
        # theta_all_deg = np.nan
        # R_all = np.nan
        # ind_walen_vals = np.nan
        ind_walen_check = np.nan
        t_walen_check_unix = np.nan

    # Get different parameters for magnetosphere and magnetosheath
    np_msp = df_mms['np'][ind_range_msp] * 1e6  # Convert to m^-3 from cm^-3
    np_msh = df_mms['np'][ind_range_msh] * 1e6  # Convert to m^-3 from cm^-3

    if coord_type == 'lmn':

        vp_lmn_vec_msp = np.array([df_mms['vp_lmn_n'][ind_range_msp],
                                   df_mms['vp_lmn_m'][ind_range_msp],
                                   df_mms['vp_lmn_l'][ind_range_msp]]) * 1e3  # Convert to m/s
        vp_lmn_vec_msp = vp_lmn_vec_msp.T

        vp_lmn_vec_msh = np.array([df_mms['vp_lmn_n'][ind_range_msh],
                                   df_mms['vp_lmn_m'][ind_range_msh],
                                   df_mms['vp_lmn_l'][ind_range_msh]]) * 1e3  # Convert to m/s
        vp_lmn_vec_msh = vp_lmn_vec_msh.T

        b_lmn_vec_msp = np.array([df_mms['b_lmn_n'][ind_range_msp],
                                  df_mms['b_lmn_m'][ind_range_msp],
                                  df_mms['b_lmn_l'][ind_range_msp]]) * 1e-9  # Convert to T from nT
        b_lmn_vec_msp = b_lmn_vec_msp.T
        b_lmn_vec_msh = np.array([df_mms['b_lmn_n'][ind_range_msh],
                                  df_mms['b_lmn_m'][ind_range_msh],
                                  df_mms['b_lmn_l'][ind_range_msh]]) * 1e-9  # Convert to T from nT
        b_lmn_vec_msh = b_lmn_vec_msh.T

        # Get the mean and median values of np, vp, and b for the magnetosphere and magnetosheath
        np_msp_median = np.nanmedian(np_msp) / 1e6  # Convert to cm^-3 from m^-3
        np_msh_median = np.nanmedian(np_msh) / 1e6  # Convert to cm^-3 from m^-3
        np_msp_mean = np.nanmean(np_msp) / 1e6  # Convert to cm^-3 from m^-3
        np_msh_mean = np.nanmean(np_msh) / 1e6  # Convert to cm^-3 from m^-3

        vp_lmn_vec_msp_median = np.nanmedian(vp_lmn_vec_msp, axis=0)  # km/sec
        vp_lmn_vec_msh_median = np.nanmedian(vp_lmn_vec_msh, axis=0)  # km/sec
        vp_lmn_vec_msp_mean = np.nanmean(vp_lmn_vec_msp, axis=0)  # km/sec
        vp_lmn_vec_msh_mean = np.nanmean(vp_lmn_vec_msh, axis=0)  # km/sec

        b_lmn_vec_msp_median = np.nanmedian(b_lmn_vec_msp, axis=0)  # T
        b_lmn_vec_msh_median = np.nanmedian(b_lmn_vec_msh, axis=0)  # T
        b_lmn_vec_msp_mean = np.nanmean(b_lmn_vec_msp, axis=0)  # T
        b_lmn_vec_msh_mean = np.nanmean(b_lmn_vec_msh, axis=0)  # T

        # Get the angle between the two vectors
        angle_b_lmn_vec_msp_msh_median = np.arccos(np.dot(b_lmn_vec_msp_median,
                                                          b_lmn_vec_msh_median) /
                                                   (np.linalg.norm(b_lmn_vec_msp_median) *
                                                    np.linalg.norm(b_lmn_vec_msh_median))
                                                   ) * 180 / np.pi
    else:
        vp_gse_vec_msp = np.array([df_mms['vp_gse_x'][ind_range_msp],
                                   df_mms['vp_gse_y'][ind_range_msp],
                                   df_mms['vp_gse_z'][ind_range_msp]]) * 1e3  # Convert km/s==> m/s
        vp_gse_vec_msp = vp_gse_vec_msp.T
        vp_gse_vec_msh = np.array([df_mms['vp_gse_x'][ind_range_msh],
                                   df_mms['vp_gse_y'][ind_range_msh],
                                   df_mms['vp_gse_z'][ind_range_msh]]) * 1e3  # Convert km/s==> m/s
        vp_gse_vec_msh = vp_gse_vec_msh.T
        b_gse_vec_msp = np.array([df_mms['b_gse_x'][ind_range_msp],
                                  df_mms['b_gse_y'][ind_range_msp],
                                  df_mms['b_gse_z'][ind_range_msp]]) * 1e-9  # Convert to T from nT
        b_gse_vec_msp = b_gse_vec_msp.T
        b_gse_vec_msh = np.array([df_mms['b_gse_x'][ind_range_msh],
                                  df_mms['b_gse_y'][ind_range_msh],
                                  df_mms['b_gse_z'][ind_range_msh]]) * 1e-9  # Convert to T from nT
        b_gse_vec_msh = b_gse_vec_msh.T

    tp_para_msp = df_mms['tp_para'][ind_range_msp] * ev_to_K  # Convert to K from ev
    tp_para_msh = df_mms['tp_para'][ind_range_msh] * ev_to_K  # Convert to K from ev
    tp_perp_msp = df_mms['tp_perp'][ind_range_msp] * ev_to_K  # Convert to K from ev
    tp_perp_msh = df_mms['tp_perp'][ind_range_msh] * ev_to_K  # Convert to K from ev

    # Get the mean and median values of temperature for the magnetosphere and magnetosheath
    tp_para_msp_median = np.nanmedian(tp_para_msp)
    tp_para_msh_median = np.nanmedian(tp_para_msh)
    tp_para_msp_mean = np.nanmean(tp_para_msp)
    tp_para_msh_mean = np.nanmean(tp_para_msh)
    tp_perp_msp_median = np.nanmedian(tp_perp_msp)
    tp_perp_msh_median = np.nanmedian(tp_perp_msh)
    tp_perp_msp_mean = np.nanmean(tp_perp_msp)
    tp_perp_msh_mean = np.nanmean(tp_perp_msh)

    # Compute the magnetosheath beta value
    beta_msh_mean = 2 * mu_0 * np_msh_mean * 1e6 * k_B * (2 * tp_para_msh_mean + tp_perp_msh_mean
                                                          ) / (3 * np.linalg.norm(
                                                            b_lmn_vec_msh_mean) ** 2)
    beta_msp_mean = 2 * mu_0 * np_msp_mean * 1e6 * k_B * (2 * tp_para_msp_mean + tp_perp_msp_mean
                                                          ) / (3 * np.linalg.norm(
                                                            b_lmn_vec_msp_mean) ** 2)

    # Check if within 2 minutes of crossing time the values went above and below the threshold
    # If ind_vals is not empty, then append the crossing time to the csv file
    # if len(ind_vals) > 0:
    if walen_relation_satisfied_v1 or jet_detection or walen_relation_satisfied_v2:
        # for xxx in range(1):
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        # _ = spd.mms.mec
        # mms_sc_pos = ptt.get_data(mms_mec_varnames[0])[1:3][0][0] / r_e
        x = np.nanmean(mms_fgm_r_gsm[:, 0]) / r_e
        y = np.nanmean(mms_fgm_r_gsm[:, 1]) / r_e
        z = np.nanmean(mms_fgm_r_gsm[:, 2]) / r_e
        r_yz = np.sqrt(y**2 + z**2)  # Projection distance in yz plane.

        # Restrict the range of theta_w from 0 to 90 degrees
        theta_w_deg_rest = theta_w_deg
        theta_w_deg_rest[theta_w_deg_rest > 90] = 180 - theta_w_deg_rest[theta_w_deg_rest > 90]
        theta_w_deg_median = np.round(np.nanmedian(theta_w_deg_rest), 3)
        # TODO: Add magnetic beta for sheath
        # List of variables to be saved in the csv file
        var_list = 'Date,Probe,walen1,walen2,jet_detection,x_gsm,y_gsm,z_gsm,r_spc,r_W,theta_w,'\
                   'jet_time,ind_min_msp,ind_max_msp,ind_min_msh,ind_max_msh,ind_jet_center,'\
                   'angle_b_lmn_vec_msp_msh_median,b_lmn_vec_msp_mean_n,b_lmn_vec_msp_mean_m,'\
                   'b_lmn_vec_msp_mean_l,b_lmn_vec_msp_median_n,b_lmn_vec_msp_median_m,'\
                   'b_lmn_vec_msp_median_l,b_lmn_vec_msh_mean_n,b_lmn_vec_msh_mean_m,'\
                   'b_lmn_vec_msh_mean_l,b_lmn_vec_msh_median_n,b_lmn_vec_msh_median_m,'\
                   'b_lmn_vec_msh_median_l,np_msp_median,np_msp_mean,np_msh_median,np_msh_mean,'\
                   'vp_lmn_vec_msp_mean_n,vp_lmn_vec_msp_mean_m,vp_lmn_vec_msp_mean_l,'\
                   'vp_lmn_vec_msp_median_n,vp_lmn_vec_msp_median_m,vp_lmn_vec_msp_median_l,'\
                   'vp_lmn_vec_msh_mean_n,vp_lmn_vec_msh_mean_m,vp_lmn_vec_msh_mean_l,'\
                   'vp_lmn_vec_msh_median_n,vp_lmn_vec_msh_median_m,vp_lmn_vec_msh_median_l,'\
                   'tp_para_msp_median,tp_para_msh_median,tp_para_msp_mean,tp_para_msh_mean,'\
                   'tp_perp_msp_median,tp_perp_msh_median,tp_perp_msp_mean,tp_perp_msh_mean,'\
                   'beta_msh_mean,beta_msp_mean'

        data_dict = {'Date': crossing_time,
                     'Probe': probe,
                     'walen1': walen_relation_satisfied_v1,
                     'walen2': walen_relation_satisfied_v2,
                     'jet_detection': jet_detection,
                     'x_gsm': np.round(x, 3),
                     'y_gsm': np.round(y, 3),
                     'z_gsm': np.round(z, 3),
                     'r_spc': np.round(r_yz, 3),
                     'r_W': np.round(np.nanmedian(R_w), 3),
                     'theta_w': theta_w_deg_median,
                     'jet_time': t_jet_center,
                     'ind_min_msp': ind_range_msp[0],
                     'ind_max_msp': ind_range_msp[-1],
                     'ind_min_msh': ind_range_msh[0],
                     'ind_max_msh': ind_range_msh[-1],
                     'ind_jet_center': ind_jet_center,
                     'angle_b_lmn_vec_msp_msh_median': np.round(angle_b_lmn_vec_msp_msh_median, 3),
                     'b_lmn_vec_msp_mean_n': np.round(b_lmn_vec_msp_mean[0] * 1e9, 3),
                     'b_lmn_vec_msp_mean_m': np.round(b_lmn_vec_msp_mean[1] * 1e9, 3),
                     'b_lmn_vec_msp_mean_l': np.round(b_lmn_vec_msp_mean[2] * 1e9, 3),
                     'b_lmn_vec_msp_median_n': np.round(b_lmn_vec_msp_median[0] * 1e9, 3),
                     'b_lmn_vec_msp_median_m': np.round(b_lmn_vec_msp_median[1] * 1e9, 3),
                     'b_lmn_vec_msp_median_l': np.round(b_lmn_vec_msp_median[2] * 1e9, 3),
                     'b_lmn_vec_msh_mean_n': np.round(b_lmn_vec_msh_median[0] * 1e9, 3),
                     'b_lmn_vec_msh_mean_m': np.round(b_lmn_vec_msh_median[1] * 1e9, 3),
                     'b_lmn_vec_msh_mean_l': np.round(b_lmn_vec_msh_median[2] * 1e9, 3),
                     'b_lmn_vec_msh_median_n': np.round(b_lmn_vec_msh_median[0] * 1e9, 3),
                     'b_lmn_vec_msh_median_m': np.round(b_lmn_vec_msh_median[1] * 1e9, 3),
                     'b_lmn_vec_msh_median_l': np.round(b_lmn_vec_msh_median[2] * 1e9, 3),
                     'np_msp_median': np.round(np_msp_median, 3),
                     'np_msp_mean': np.round(np_msp_mean, 3),
                     'np_msh_median': np.round(np_msh_median, 3),
                     'np_msh_mean': np.round(np_msh_mean, 3),
                     'vp_lmn_vec_msp_mean_n': np.round(vp_lmn_vec_msp_mean[0] / 1e3, 3),
                     'vp_lmn_vec_msp_mean_m': np.round(vp_lmn_vec_msp_mean[1] / 1e3, 3),
                     'vp_lmn_vec_msp_mean_l': np.round(vp_lmn_vec_msp_mean[2] / 1e3, 3),
                     'vp_lmn_vec_msp_median_n': np.round(vp_lmn_vec_msp_median[0] / 1e3, 3),
                     'vp_lmn_vec_msp_median_m': np.round(vp_lmn_vec_msp_median[1] / 1e3, 3),
                     'vp_lmn_vec_msp_median_l': np.round(vp_lmn_vec_msp_median[2] / 1e3, 3),
                     'vp_lmn_vec_msh_mean_n': np.round(vp_lmn_vec_msh_mean[0] / 1e3, 3),
                     'vp_lmn_vec_msh_mean_m': np.round(vp_lmn_vec_msh_mean[1] / 1e3, 3),
                     'vp_lmn_vec_msh_mean_l': np.round(vp_lmn_vec_msh_mean[2] / 1e3, 3),
                     'vp_lmn_vec_msh_median_n': np.round(vp_lmn_vec_msh_median[0] / 1e3, 3),
                     'vp_lmn_vec_msh_median_m': np.round(vp_lmn_vec_msh_median[1] / 1e3, 3),
                     'vp_lmn_vec_msh_median_l': np.round(vp_lmn_vec_msh_median[2] / 1e3, 3),
                     'tp_para_msp_median': np.round(tp_para_msp_median, 3),
                     'tp_para_msh_median': np.round(tp_para_msh_median, 3),
                     'tp_para_msp_mean': np.round(tp_para_msp_mean, 3),
                     'tp_para_msh_mean': np.round(tp_para_msh_mean, 3),
                     'tp_perp_msp_median': np.round(tp_perp_msp_median, 3),
                     'tp_perp_msh_median': np.round(tp_perp_msh_median, 3),
                     'tp_perp_msp_mean': np.round(tp_perp_msp_mean, 3),
                     'tp_perp_msh_mean': np.round(tp_perp_msh_mean, 3),
                     'beta_msh_mean': np.round(beta_msh_mean, 3),
                     'beta_msp_mean': np.round(beta_msp_mean, 3)
                     }

        if x > -5 and x < 12 and r_yz < 12:
            # Check if the file exists, if not then create it
            if not os.path.isfile(fname):
                with open(fname, 'w') as f:
                    f.write(var_list + '\n')
                    if verbose:
                        print(f'File {fname} created')
            # Check if the file already contains the data corresponding to the crossing time
            if os.path.isfile(fname):
                df_added_list = pd.read_csv(fname, sep=',', index_col=False)
                if not np.any(df_added_list['Date'].values == str(crossing_time)):
                    with open(fname, 'a') as f:
                        for key in data_dict.keys():
                            f.write(f'{data_dict[key]},')
                        f.write('\n')
                    if verbose:
                        print(f'Data for all keys written to file {fname}')
                    f.close()

    tplot_fnc(ptt=ptt, probe=probe, data_rate=data_rate, df_mms=df_mms,
              ind_range_msp=ind_range_msp, ind_range_msh=ind_range_msh, t_jet_center=t_jet_center,
              walen_v1=walen_relation_satisfied_v1, walen_v2=walen_relation_satisfied_v2,
              jet_detection=jet_detection, ind_crossing=ind_crossing,
              shear_val=int(angle_b_lmn_vec_msp_msh_median), date_obs=date_obs)

    return None


def check_walen_relation(df_mms=None, t_jet_center=None, dt_walen=30, coord_type="gsm",
                         ind_msh=None, n_points_walen=10, time_cadence_median=0.15,
                         jet_detection=False,):
    """
    Check if the Walen relation is satisfied for a given index values in the dataframe.

    Parameters
    ----------
    df_mms : pandas.DataFrame
        Dataframe containing the MMS data, both the FPI and FPGM
    t_jet_center : datetime
        Time of the jet center. This is the time of the center of the jet.
    dt_walen : int
        Time interval in seconds to check the Walen relation
    coord_type : str
        Coordinate system to use for the Walen relation. Default is "gsm". Other options are "gse"
        and "lmn".
    ind_msh : list
        List of indices for the magnetosheath.
    n_points_walen : int
        Number of points to use for the Walen relation. Default is 10.
    time_cadence_median : float
        Median time cadence of the data. Default is 0.15 seconds.
    jet_detection : bool
        Whether the jet was detected or not. Default is False.

    Returns
    -------
    walen_relation_satisfied_v1 : bool
        Whether the Walen relation, version 1 is satisfied or not.
    walen_relation_satisfied_v2 : bool
        Whether the Walen relation, version 2 is satisfied or not.
    theta_w_deg : numpy.ndarray
        The theta parameter for the Walen relation, restricted to the region for which the Walen
        relation was checked
    R_w : numpy.ndarray
        The R parameter for the Walen relation, restricted to the region for which the Walen
        relation was checked
    theta_all_deg : numpy.ndarray
        The theta parameter for the Walen relation, for whole dataframe
    R_all : numpy.ndarray
        The R parameter for the Walen relation, for whole dataframe
    ind_walen_vals: numpy.ndarray
        Indicess within the jet region for which the Walen relation was satisfied continuously for
        n_points_walen
    ind_walen_check: numpy.ndarray
        Indices within the jet region for which the Walen relation was checked
    """

    # Define the constants
    mu_0 = 4 * np.pi * 1e-7  # Permeability of free space (H/m)
    k_B = 1.38064852e-23  # Boltzmann constant in J/K
    m_p = 1.6726219e-27  # Mass of proton in kg
    ev_to_K = 11604.505  # Conversion factor from eV to K

    n_points_walen = int(n_points_walen / time_cadence_median)

    # Define ind_walen_check as indices corresponding to 30 seconds before and after t_jet_center
    t_jet_center_minus_30_sec = t_jet_center - datetime.timedelta(seconds=dt_walen)
    ind_jet_center_minus_30_sec = np.argmin(np.abs(df_mms.index - t_jet_center_minus_30_sec))
    t_jet_center_plus_30_sec = t_jet_center + datetime.timedelta(seconds=dt_walen)
    ind_jet_center_plus_30_sec = np.argmin(np.abs(df_mms.index - t_jet_center_plus_30_sec))

    ind_walen_check = np.arange(ind_jet_center_minus_30_sec, ind_jet_center_plus_30_sec)

    alpha_jet = np.full(len(ind_walen_check), np.nan)
    alpha_msh = np.full(len(ind_msh), np.nan)
    v_th_jet = np.full((len(ind_walen_check), 3), np.nan)
    v_th_msh = np.full((len(ind_msh), 3), np.nan)

    if coord_type == 'lmn':
        for i, xx in enumerate(ind_walen_check):
            alpha_jet[i] = (mu_0 * df_mms['np'][xx] * 1e6 * k_B) * (
                df_mms['tp_para'][xx] - df_mms['tp_perp'][xx]) * ev_to_K / (1e-18 * (
                    df_mms['b_lmn_l'][xx] ** 2 + df_mms['b_lmn_m'][xx] ** 2 +
                    df_mms['b_lmn_n'][xx] ** 2))

            b_lmn_vec_jet = np.array([df_mms['b_lmn_l'][xx], df_mms['b_lmn_m'][xx],
                                      df_mms['b_lmn_n'][xx]]) * 1e-9
            for j in range(3):
                v_th_jet[i, j] = b_lmn_vec_jet[j] * ((1 - alpha_jet[i]) / (mu_0 * df_mms['np'][xx]
                                                                           * 1e6 * m_p))**0.5
        for i, xx in enumerate(ind_msh):
            alpha_msh[i] = (mu_0 * df_mms['np'][xx] * 1e6 * k_B) * (
                df_mms['tp_para'][xx] - df_mms['tp_perp'][xx]) * ev_to_K / (
                1e-18 *
                (df_mms['b_lmn_l'][xx] ** 2
                 + df_mms['b_lmn_m'][xx] ** 2
                 + df_mms['b_lmn_n'][xx] ** 2
                 )
            )

            b_lmn_vec_msh = np.array([df_mms['b_lmn_l'][xx], df_mms['b_lmn_m'][xx],
                                      df_mms['b_lmn_n'][xx]]) * 1e-9
            for j in range(3):
                v_th_msh[i, j] = b_lmn_vec_msh[j] * ((1 - alpha_msh[i]) / (mu_0 * df_mms['np'][xx]
                                                                           * 1e6 * m_p))**0.5

    alpha_all = mu_0 * df_mms['np'] * 1e6 * k_B * (df_mms['tp_para'] - df_mms['tp_perp']
                                                   ) * ev_to_K / (1e-18 * (df_mms['b_lmn_l'] ** 2
                                                                  + df_mms['b_lmn_m'] ** 2
                                                                  + df_mms['b_lmn_n'] ** 2))

    v_th_all = np.full((len(df_mms['tp_para']), 3), np.nan)
    for i in range(len(df_mms['tp_para'])):
        v_th_all[i, :] = 1e-9 * np.array([df_mms['b_lmn_l'][i], df_mms['b_lmn_m'][i],
                                          df_mms['b_lmn_n'][i]]) * ((1 - alpha_all[i]) / (
                                                                     mu_0 * df_mms['np'][i] * 1e6
                                                                     * m_p))**0.5
    delta_v_th = np.nanmean(v_th_msh, axis=0) - v_th_jet
    delta_v_th_all = np.nanmean(v_th_msh, axis=0) - v_th_all

    vp_lmn_vec_walen = np.array([df_mms['vp_lmn_l'][ind_walen_check],
                                 df_mms['vp_lmn_m'][ind_walen_check],
                                 df_mms['vp_lmn_n'][ind_walen_check]]).T * 1e3  # km/s to m/s
    vp_lmn_vec_all = np.array([df_mms['vp_lmn_l'], df_mms['vp_lmn_m'], df_mms['vp_lmn_n']]).T * 1e3

    if coord_type == 'lmn':
        # Find the proton bulk velocity in the msh region and then compute its average
        vp_lmn_vec_msh = np.array([df_mms['vp_lmn_l'][ind_msh], df_mms['vp_lmn_m'][ind_msh],
                                   df_mms['vp_lmn_n'][ind_msh]]) * 1e3  # Convert to m/s from km/s
        vp_lmn_vec_msh_mean = np.nanmean(vp_lmn_vec_msh, axis=1)  # km/sec

        delta_v_obs = vp_lmn_vec_msh_mean - vp_lmn_vec_walen
        delta_v_obs_all = vp_lmn_vec_msh_mean - vp_lmn_vec_all

    # Compute the angle between the observed and the theoretical velocity jumps
    theta_w = np.full(len(ind_walen_check), np.nan)
    for i in range(len(ind_walen_check)):
        theta_w[i] = np.arccos(np.dot(delta_v_th[i, :], delta_v_obs[i, :]) / (
                               np.linalg.norm(delta_v_th[i, :]) *
                               np.linalg.norm(delta_v_obs[i, :])
                               )
                               )

    theta_all = np.full(len(df_mms['tp_para']), np.nan)
    for i in range(len(df_mms['tp_para'])):
        theta_all[i] = np.arccos(np.dot(delta_v_th_all[i, :], delta_v_obs_all[i, :]) / (
                               np.linalg.norm(delta_v_th_all[i, :]) *
                               np.linalg.norm(delta_v_obs_all[i, :])
                               )
                               )

    # Convert angle to degrees
    theta_w_deg = theta_w * 180 / np.pi
    theta_all_deg = theta_all * 180 / np.pi

    # Compute the ratio of the observed and the theoretical velocity jumps
    R_w = np.linalg.norm(delta_v_th, axis=1) / np.linalg.norm(delta_v_obs, axis=1)
    R_all = np.linalg.norm(delta_v_th_all, axis=1) / np.linalg.norm(delta_v_obs_all, axis=1)

    # print(f"Values of theta_w_deg: {theta_w_deg}")
    # print(f"Values of R_w: {R_w}")
    # Check if Walen relation is satisfied (0.4 < R_w < 3 and 0 < theta_w < 30 or 150 < theta_w <
    # 180) in green color

    bool_array_gt = np.logical_and(np.logical_and(0.4 < R_w, R_w < 3),
                                   np.logical_and(0 <= theta_w_deg, theta_w_deg <= 30))

    bool_array_lt = np.logical_and(np.logical_and(0.4 < R_w, R_w < 3),
                                   np.logical_and(150 <= theta_w_deg, theta_w_deg <= 180))

    # Check the number of True in the bool_array
    num_true_gt = np.sum(bool_array_gt)
    num_true_lt = np.sum(bool_array_lt)

    # If more than 50% of the data points satisfy the Walen relation, then the Walen relation
    # version 1 is assumed to be satisfied
    if ((num_true_gt + num_true_lt) >= len(ind_walen_check)/2 and num_true_gt > 0 and
            num_true_lt > 0):
        walen_relation_satisfied_v1 = True
    else:
        walen_relation_satisfied_v1 = False

    # Check if the Walen relation is satisfied (0.4 < R_w < 3 and 0 < theta_w < 30 or 150 < theta_w
    # < 180) continuously at least once for same length of time as the jet length

    ind_walen_vals_gt = np.flatnonzero(np.convolve(bool_array_gt > 0,
                                       np.ones(n_points_walen, dtype=int),
                                      'valid') >= n_points_walen)
    ind_walen_vals_lt = np.flatnonzero(np.convolve(bool_array_lt > 0,
                                       np.ones(n_points_walen, dtype=int),
                                      'valid') >= n_points_walen)

    # Set ind_walen_vals to union of ind_walen_vals_gt and ind_walen_vals_lt
    ind_walen_vals = np.union1d(ind_walen_vals_gt, ind_walen_vals_lt)

    # If Walen relation is satisfied continuously for same length of time as the jet length, then
    # Walen relation version 2 is assumed to be satisfied
    if len(ind_walen_vals_gt) > 0 and len(ind_walen_vals_lt) > 0:
        walen_relation_satisfied_v2 = True
    else:
        walen_relation_satisfied_v2 = False

    if walen_relation_satisfied_v1:
        print('\033[92m \n Walen relation satisfied \033[0m \n')
        print(f'\033[92m R_w: {np.nanmedian(R_w):.3f} \033[0m')
        print(f'\033[92m theta_w_deg: {np.nanmedian(theta_w_deg):.3f} \033[0m')
        if jet_detection:
            print(f'\033[92m Jet detection: {jet_detection} \033[0m')
        else:
            print(f'\033[91m Jet detection: {jet_detection} \033[0m')
    else:
        print('\033[91m \n Walen relation not satisfied \033[0m \n')
        print(f'\033[91m R_w: {np.nanmedian(R_w):.3f} \033[0m')
        print(f'\033[91m theta_w_deg: {np.nanmedian(theta_w_deg):.3f} \033[0m')
        if jet_detection:
            print(f'\033[92m Jet detection: {jet_detection} \033[0m')
        else:
            print(f'\033[91m Jet detection: {jet_detection} \033[0m')

    if walen_relation_satisfied_v2:
        print('\033[92m Walen relation version 2 is satisfied \033[0m \n')
    else:
        print('\033[91m Walen relation version 2 is not satisfied \033[0m \n')

    return (walen_relation_satisfied_v1, walen_relation_satisfied_v2, theta_w_deg, R_w,
            theta_all_deg, R_all, ind_walen_vals, ind_walen_check)


def check_jet_location(df_mms=None, jet_len=3, time_cadence_median=0.15, v_thresh=70,
                       ind_msh=None, verbose=True, ind_crossing=None, date_obs=None):
    """
    This function checks if a jet is present in the given data frame.

    Parameters
    ----------
    df_mms : pandas.DataFrame
        Data frame containing the MMS data
    jet_len : float
        Length of the time for which jet must have threshold velocity to be considered as a jet,
        both at the positive and negative side of the crossing point
    time_cadence_median : float
        Median time cadence of the data in seconds
    v_thresh : float
        Threshold velocity in km/s. Default is 70 km/s
    ind_msh : numpy.ndarray
        Indices of the MSH data
    verbose : bool
        If True, prints the jet location and the jet length
    ind_crossing : numpy.ndarray
        Indices of the crossing points
    date_obs : str
        Date of observation in the format 'YYYY-MM-DD'

    Returns
    -------
    jet_present : bool
        True if jet is present, False otherwise
    delta_v_min : float
        Difference between the velocity with respect to a reference magnetosheath velocity and the
        time series from 1 minute before/after the jet center to the point just before/after the jet
        starts/ends
    delta_v_max : float
        Difference between the velocity with respect to a reference magnetosheath velocity and the
        time series from 1 minute before/after the jet center to the point just before/after the jet
        starts/ends
    t_jet_center : datetime.datetime
        Time of the jet center in datetime format
    ind_jet_center : int
        Index of the jet center in the dataframe
    ind_jet_center_minus_1_min: int
        Index of the point 1 minute before the jet center in the dataframe
    ind_jet_center_plus_1_min: int
        Index of the point 1 minute after the jet center in the dataframe
    vp_lmn_diff_l: np.ndarray
        Difference between the velocity with respect to a reference magnetosheath velocity and the
        'l' component of the velocity
    """
    # Compute the number of points corresponding to jet_len
    n_points_jet = int(jet_len / time_cadence_median)

    # Define the distance within which maximum and minimum values of jet velocity must lie
    # to be considered as a jet (within 30 seconds)
    delta_jet_min_max_ind = int(60 / time_cadence_median)
    # Get the median value of the velocity corresponding to the magnetosheath in the lmn coordinate
    # system
    vp_lmn_vec_msh_median = np.array([np.nanmedian(df_mms['vp_lmn_l'][ind_msh]),
                                      np.nanmedian(df_mms['vp_lmn_m'][ind_msh]),
                                      np.nanmedian(df_mms['vp_lmn_n'][ind_msh])])

    # Subtract the median value of the velocity corresponding to the magnetosheath from the velocity
    # in the lmn coordinate system
    vp_lmn_diff_l = df_mms['vp_lmn_l'] - vp_lmn_vec_msh_median[0]
    # vp_lmn_diff_m = df_mms['vp_lmn_m'] - vp_lmn_vec_msh_mean[1]
    # vp_lmn_diff_n = df_mms['vp_lmn_n'] - vp_lmn_vec_msh_mean[2]

    # Add the difference in the velocity in the lmn coordinate system to the dataframe
    df_mms['vp_lmn_diff_l'] = vp_lmn_diff_l

    print('vp_lmn_diff_l added to the dataframe')

    # Find the index where vp_lmn_diff_l has the maximum value
    ind_jet_max = np.argmax(vp_lmn_diff_l)
    t_jet_max = df_mms.index[ind_jet_max]

    # Within 30 seconds of ind_jet_max find the index where vp_lmn_diff_l has the minimum value
    ind_jet_min = ind_jet_max - delta_jet_min_max_ind + \
        np.argmin(
            vp_lmn_diff_l[ind_jet_max - delta_jet_min_max_ind:ind_jet_max + delta_jet_min_max_ind])
    t_jet_min = df_mms.index[ind_jet_min]

    if t_jet_max < t_jet_min:
        # Find the time centered between t_jet_max and t_jet_min, which are datetime objects
        t_jet_center = t_jet_max + (t_jet_min - t_jet_max) / 2
        ind_jet_center = np.argmin(np.abs(df_mms.index - t_jet_center))

        # Find the time difference between t_jet_max and t_jet_center
        delta_t_jet_max_center = t_jet_center - t_jet_max
        # Find the number of points corresponding to delta_t_jet_max_center
        delta_n_jet_max_center = int(delta_t_jet_max_center.total_seconds() / time_cadence_median)

        # Find the time difference between t_jet_min and t_jet_center
        delta_t_jet_min_center = t_jet_min - t_jet_center
        # Find the number of points corresponding to delta_t_jet_min_center
        delta_n_jet_min_center = int(delta_t_jet_min_center.total_seconds() / time_cadence_median)

        # Find the median value of vp_lmn_diff_l from 1 minute before t_jet_center until just before
        # the jet starts
        t_jet_center_minus_1_min = t_jet_center - datetime.timedelta(minutes=1)
        ind_jet_center_minus_1_min = np.argmin(np.abs(df_mms.index - t_jet_center_minus_1_min))
        v_max_median = np.nanmedian(vp_lmn_diff_l[ind_jet_center_minus_1_min:
                                                  (ind_jet_center - 2 * delta_n_jet_max_center)])

        # Find the median value of vp_lmn_diff_l between t_jet_center and 1 minute after
        # t_jet_center
        t_jet_center_plus_1_min = t_jet_center + datetime.timedelta(minutes=1)
        ind_jet_center_plus_1_min = np.argmin(np.abs(df_mms.index - t_jet_center_plus_1_min))
        v_min_median = np.nanmedian(vp_lmn_diff_l[(ind_jet_center + 2 * delta_n_jet_min_center):
                                                  ind_jet_center_plus_1_min])

        # Subtract v_max_median from vp_lmn_diff_l from t_jet_cener_minus_1_min to
        # t_jet_center
        delta_v_max = vp_lmn_diff_l[ind_jet_center_minus_1_min:ind_jet_center] - v_max_median

        # Subtract v_min_median from vp_lmn_diff_l from t_jet_center to t_jet_center_plus_1_min
        delta_v_min = vp_lmn_diff_l[ind_jet_center:ind_jet_center_plus_1_min] - v_min_median

        # Check if delta_v_max has a sustained value greater than v_thresh for n_points_jet
        # points
        ind_v_gt_vthresh = np.flatnonzero(np.convolve(delta_v_max >= v_thresh,
                                                      np.ones(n_points_jet, dtype=int),
                                                      'valid') >= n_points_jet)
        # Check if delta_v_min has a sustained value less than -v_thresh for n_points_jet
        # points
        ind_v_lt_vthresh = np.flatnonzero(np.convolve(delta_v_min <= -v_thresh,
                                                      np.ones(n_points_jet, dtype=int),
                                                      'valid') >= n_points_jet)
        # If both conditions are satisfied then a jet is found
        if len(ind_v_gt_vthresh) > 0 and len(ind_v_lt_vthresh) > 0:
            jet_detection = True
        else:
            jet_detection = False
            # ind_jet_center = []
            # delta_v_min = []
            # delta_v_max = []
    else:
        t_jet_center = t_jet_min + (t_jet_max - t_jet_min) / 2
        ind_jet_center = np.argmin(np.abs(df_mms.index - t_jet_center))

        # Find the time difference between t_jet_max and t_jet_center
        delta_t_jet_max_center = t_jet_max - t_jet_center
        # Find the number of points corresponding to delta_t_jet_max_center
        delta_n_jet_max_center = int(delta_t_jet_max_center.total_seconds() / time_cadence_median)

        # Find the time difference between t_jet_min and t_jet_center
        delta_t_jet_min_center = t_jet_center - t_jet_min
        # Find the number of points corresponding to delta_t_jet_min_center
        delta_n_jet_min_center = int(delta_t_jet_min_center.total_seconds() / time_cadence_median)

        # Find the median value of vp_lmn_diff_l from 1 minute before t_jet_center to just before
        # the jet starts
        t_jet_center_minus_1_min = t_jet_center - datetime.timedelta(minutes=1)
        ind_jet_center_minus_1_min = np.argmin(np.abs(df_mms.index - t_jet_center_minus_1_min))
        v_min_median = np.nanmedian(vp_lmn_diff_l[ind_jet_center_minus_1_min:
                                                  (ind_jet_center - 2 * delta_n_jet_min_center)])

        # Find the median value of vp_lmn_diff_l between just following the jet and 1 minute after
        # t_jet_center
        t_jet_center_plus_1_min = t_jet_center + datetime.timedelta(minutes=1)
        ind_jet_center_plus_1_min = np.argmin(np.abs(df_mms.index - t_jet_center_plus_1_min))
        v_max_median = np.nanmedian(vp_lmn_diff_l[(ind_jet_center + 2 * delta_n_jet_max_center):
                                                  ind_jet_center_plus_1_min])

        # Subtract v_min_median from vp_lmn_diff_l from t_jet_cener_minus_1_min to
        # t_jet_center
        delta_v_min = vp_lmn_diff_l[ind_jet_center_minus_1_min:ind_jet_center] - v_min_median

        # Subtract v_max_median from vp_lmn_diff_l from t_jet_center to t_jet_center_plus_1_min
        delta_v_max = vp_lmn_diff_l[ind_jet_center:ind_jet_center_plus_1_min] - v_max_median

        # Check if delta_v_min has a sustained value less than -v_thresh for n_points_jet
        # points
        ind_v_lt_vthresh = np.flatnonzero(np.convolve(delta_v_min <= -v_thresh,
                                                      np.ones(n_points_jet, dtype=int),
                                                      'valid') >= n_points_jet)
        # Check if delta_v_max has a sustained value greater than v_thresh for n_points_jet
        # points
        ind_v_gt_vthresh = np.flatnonzero(np.convolve(delta_v_max >= v_thresh,
                                                      np.ones(n_points_jet, dtype=int),
                                                      'valid') >= n_points_jet)

        # If both conditions are satisfied then a jet is found
        if len(ind_v_lt_vthresh) > 0 and len(ind_v_gt_vthresh) > 0:
            jet_detection = True
        else:
            jet_detection = False
            # ind_jet_center = []
            # delta_v_min = []
            # delta_v_max = []
    if verbose:
        if jet_detection:
            print(f'\n\033[1;32m Jet found at {t_jet_center}\033[0m \n')
        else:
            print(f'\n\033[1;31m No jet found at {t_jet_center}\033[0m \n')

    # Plot vp_lmn_diff_l, delta_v_max, delta_v_min, v_max_median, and v_min_median
    plt.figure(figsize=(10, 5))
    plt.plot(vp_lmn_diff_l, 'r-', alpha=0.3)
    plt.plot(delta_v_max, 'b-', alpha=1, label='$\\Delta v_{\\rm max}$')
    plt.plot(delta_v_min, 'g-', alpha=1, label='$\\Delta v_{\\rm min}$')
    # Draw a vertical line at t_jet_center
    plt.axvline(t_jet_center, color='k', linestyle='--', alpha=0.5)
    # Draw a vertical line at t_jet_center_minus_1_min
    plt.axvline(t_jet_center_minus_1_min, color='m', linestyle='--', alpha=0.5)
    # Draw a vertical line at t_jet_center_plus_1_min
    plt.axvline(t_jet_center_plus_1_min, color='g', linestyle='--', alpha=0.5)
    # Draw a horizontal line at v_thresh
    plt.axhline(v_thresh, color='k', linestyle='--', alpha=0.5)
    # Draw a horizontal line at -v_thresh
    plt.axhline(-v_thresh, color='k', linestyle='--', alpha=0.5)
    plt.title("Different delta as a function of time at"
              f" {t_jet_center.strftime('%Y-%m-%d %H:%M:%S')}")
    plt.xlabel("Time [UTC]")
    plt.ylabel("$\\Delta V$ [km/s]")
    plt.legend(loc=1)
    # On plot write jet detection status
    if jet_detection:
        plt.text(0.05, 0.95, ind_crossing, horizontalalignment='left',
                 verticalalignment='top', transform=plt.gca().transAxes, color='g')
        save_folder = f"../figures/jet_reversal_checks/check_{date_obs}/delta_v/jet/"
    else:
        plt.text(0.05, 0.95, ind_crossing, horizontalalignment='left',
                 verticalalignment='top', transform=plt.gca().transAxes, color='r')
        save_folder = f"../figures/jet_reversal_checks/check_{date_obs}/delta_v/no_jet/"

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    fig_name = f"{save_folder}/delta_v_{t_jet_center.strftime('%Y-%m-%d_%H-%M-%S')}.png"

    plt.savefig(fig_name, dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True,
                facecolor='w', edgecolor='w', orientation='landscape')

    return (jet_detection, delta_v_min, delta_v_max, t_jet_center, ind_jet_center,
            ind_jet_center_minus_1_min, ind_jet_center_plus_1_min, vp_lmn_diff_l)


def check_msp_msh_location(df_mms=None, time_cadence_median=0.15, verbose=True):
    """
    Try to find the location of the msh and msp

    Parameters
    ----------
    df_mms : pandas.DataFrame
        Dataframe containing the mms data
    time_cadence_median : float
        Median time cadence of the data
    verbose : bool
        If True, print information

    Returns
    -------
    ind_range_msh : numpy.ndarray
        Location of the msh
    ind_range_msp : numpy.ndarray
        Location of the msp
    """
    # TODO: Check if threshold value of 5 and 10 ares fine or if we need to decrease/increase it
    n_thresh_msp = 5
    n_thresh_msh = 10

    # Compute the median time difference between the points
    # NOTE: The reason for factor 10**9 is to convert the time difference to seconds from
    # nanoseconds which is the format in which df_mms.index.astype(np.int64) outputs time
    time_cadence_median = np.median(np.diff(df_mms.index.astype(np.int64) / 10**9))

    # TODO: Instead of one difference between indices, use the median value of differences between
    # all the indices
    n_points_msp = int(5 / time_cadence_median)
    n_points_msh = int(5 / time_cadence_median)

    # n_points_walen = int(10 / time_cadence_median)

    # Find an interval of length at least 'n_points_msp_msh' where 'np' is greater than 'n_thresh'
    # on both sides of the minimum value
    np_msp_bool_array = df_mms['np'] < n_thresh_msp
    np_msh_bool_array = df_mms['np'] > n_thresh_msh

    # ind_np_msp_vals = np.flatnonzero(np.convolve(np_msp_bool_array > 0,
    #                                              np.ones(
    #                                                  n_points_msp, dtype=int),
    #                                              'valid') >= n_points_msp)

    result_msp = list(mit.run_length.encode(np_msp_bool_array))
    # Find the length of longest subsequence of True, and the location if that index in result
    max_true_count_msp = -1
    max_true_idx_msp = -1
    for idx, (val, count) in enumerate(result_msp):
        if val and max_true_count_msp < count:
            max_true_count_msp = count
            max_true_idx_msp = idx
    # Find total elements before and after the longest subsequence tuple
    elems_before_idx_msp = sum((idx[1] for idx in result_msp[:max_true_idx_msp]))
    # elems_after_idx_msp = sum((idx[1] for idx in result_msp[max_true_idx_msp + 1:]))

    # Check if the longest subsequence is greater than the threshold
    if max_true_count_msp >= n_points_msp:
        ind_min_msp = int(elems_before_idx_msp + (max_true_count_msp - n_points_msp) / 2)
        ind_max_msp = int(elems_before_idx_msp + (max_true_count_msp + n_points_msp) / 2)
        ind_range_msp = np.arange(ind_min_msp, ind_max_msp)

    # ind_np_msh_vals = np.flatnonzero(np.convolve(np_msh_bool_array > 0,
    #                                  np.ones(n_points_msh, dtype=int),
    #                                  'valid') >= n_points_msh)

    result_msh = list(mit.run_length.encode(np_msh_bool_array))
    # Find the length of longest subsequence of True, and the location if that index in result
    max_true_count_msh = -1
    max_true_idx_msh = -1
    for idx, (val, count) in enumerate(result_msh):
        if val and max_true_count_msh < count:
            max_true_count_msh = count
            max_true_idx_msh = idx
    # Find total elements before and after the longest subsequence tuple
    elems_before_idx_msh = sum((idx[1]
                               for idx in result_msh[:max_true_idx_msh]))
    # elems_after_idx_msh = sum((idx[1] for idx in result_msh[max_true_idx_msh + 1:]))

    # Check if the longest subsequence is greater than the threshold
    if max_true_count_msh >= n_points_msh:
        ind_min_msh = int(elems_before_idx_msh + (max_true_count_msh - n_points_msh) / 2)
        ind_max_msh = int(elems_before_idx_msh + (max_true_count_msh + n_points_msh) / 2)
        ind_range_msh = np.arange(ind_min_msh, ind_max_msh)

    if verbose:
        try:
            print(f"ind_min_msp: {ind_min_msp}")
            print(f"ind_max_msp: {ind_max_msp}")
            print(f"ind_min_msh: {ind_min_msh}")
            print(f"ind_max_msh: {ind_max_msh}")
        except Exception:
            pass

    return ind_range_msp, ind_range_msh


def tplot_fnc(ptt=None, probe=3, data_rate='brst', df_mms=None, ind_range_msp=None,
              ind_range_msh=None, t_jet_center=None, walen_v1=False, walen_v2=False,
              jet_detection=False, ind_crossing=None, shear_val=None, date_obs=None,):
    """
    Plot the data from the MMS spacecraft along with walen test and jet detection results

    Parameters
    ----------
    ptt : str
        The pytplot object
    probe : int
        The probe number. Default is 3
    data_rate : str
        The data rate. Default is 'brst'
    df_mms : pandas.DataFrame
        The dataframe containing the MMS data
    ind_range_msp : numpy.ndarray
        The indices of the MSP
    ind_range_msh : numpy.ndarray
        The indices of the MSH
    t_jet_center : datetime.datetime
        The time of the jet center
    walen_v1 : bool
        The status of the walen test v1
    walen_v2 : bool
        The status of the walen test v2
    jet_detection : bool
        The status of the jet detection
    ind_crossing : int
        The index of the crossing date
    shear_val : float
        The shear value between the magnetosheath and the magnetosphere
    date_obs : datetime.datetime
        The observation date

    Returns
    -------
    None
    """
    # Set the fontstyle to Times New Roman
    font = {'family': 'serif', 'weight': 'normal', 'size': 12}
    plt.rc('font', **font)
    plt.rc('text', usetex=False)
    # Set tick parameters such that ticks are visible on all sides
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True

    # Set the tick length and widths
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.minor.width'] = 1

    # Define the size of the figure
    plt.rcParams['figure.figsize'] = [20, 6]

    # plt.style.use('dark_background')

    # if jet_detection:

    tplot_global_options = {"show_all_axes": True,
                            "black_background": True,
                            "crosshair": True,
                            "vertical_spacing": 0,
                            "wsize": [2500, 1080],
                            # "title": f"Probe {probe} {data_rate} data",
                            }
    for key in tplot_global_options:
        ptt.tplot_options(key, tplot_global_options[key])

    keys_to_plot = [f'mms{probe}_dis_energyspectr_omni_{data_rate}',
                    f'mms{probe}_des_energyspectr_omni_{data_rate}',
                    f'mms{probe}_dis_numberdensity_{data_rate}',
                    'Tp',
                    'mms3_fgm_b_lmn_srvy_l2',
                    f'mms{probe}_dis_bulkv_lmn_{data_rate}',
                    'delta_v_vp_lmn_diff_l',
                    'R_w',
                    'theta_w_deg',
                    ]

    # ptt.timebar(ptt.get_data(mms_fpi_varnames[0])[0].min() + 200, color='red', dash=True, thick=2)
    # ptt.timespan(ptt.get_data(mms_fpi_varnames[0])[0].min() + 160, 100, keyword="seconds")
    ion_energy_spectr_dict_option = {'Colormap': "Spectral_r",
                                     'ylog': True,
                                     'zlog': True,
                                     'ytitle': "$p^+$",
                                     'ysubtitle': "[eV]",
                                     'ztitle': "Counts",
                                     }

    electron_energy_spectr_dict_option = {'Colormap': "Spectral_r",
                                          'ylog': True,
                                          'zlog': True,
                                          'ytitle': "$e^-$",
                                          'ysubtitle': "[eV]",
                                          'ztitle': "Counts",
                                          }

    number_density_dict_option = {'ytitle': '$n_{\\rm p}$',
                                  'ysubtitle': '[cm$^{-3}$]',
                                  'ylog': True,
                                  }

    tp_dict_option = {'ylog': True,
                      'color': ['red', 'blue'],
                      'linestyle': '-',
                      'legend_names': ['$\\parallel$', '$\\perp$'],
                      'ytitle': '$T_{\\rm p}$',
                      'ysubtitle': '[eV]',
                      }

    b_dict_option = {'ytitle': '$B$',
                     'ysubtitle': '[nT]',
                     'color': ['blue', 'green', 'red'],
                     'legend_names': ['L', 'M', 'N'],
                     'linestyle': '-',
                     }

    bulkv_dict_option = {'ytitle': '$v_{\\rm p}$',
                         'ysubtitle': '[km/s]',
                         'color': ['blue', 'green', 'red'],
                         'linestyle': '-',
                         'legend_names': ['L', 'M', 'N'],
                         }

    print(f"Jet detection: {jet_detection}")
    if jet_detection:
        if walen_v1 and walen_v2:
            delta_v_line_color = 'green'
        elif walen_v1 and not walen_v2:
            delta_v_line_color = 'blue'
        elif not walen_v1 and walen_v2:
            delta_v_line_color = 'magenta'
        else:
            delta_v_line_color = 'red'
    else:
        delta_v_line_color = 'black'

    # Define the limits of the delta v plot
    dv_min_min = np.nanmin(ptt.get_data("delta_v_min")[1:])
    dv_min_max = np.nanmax(ptt.get_data("delta_v_min")[1:])
    dv_max_min = np.nanmin(ptt.get_data("delta_v_max")[1:])
    dv_max_max = np.nanmax(ptt.get_data("delta_v_max")[1:])
    dv_vpl_min = np.nanmin(ptt.get_data("vp_lmn_diff_l")[1:])
    dv_vpl_max = np.nanmax(ptt.get_data("vp_lmn_diff_l")[1:])

    dv_min = 1.1 * min(dv_min_min, dv_max_min, dv_vpl_min)
    dv_max = 1.1 * max(dv_min_max, dv_max_max, dv_vpl_max)

    delta_v_dict_option = {'color': ['green', 'blue', 'red'],
                           'linestyle': ['-', '--', '--'],
                           'lw': [2, 1, 5],
                           'yrange': [dv_min, dv_max],
                           'ytitle': '$\\Delta v$',
                           'ysubtitle': 'km/s',
                           'legend_names': ['$\\Delta v_{min}$', '$\\Delta v_{max}$',
                                      '$\\Delta v_{p,L}$'],
                           }

    r_w_dict_option = {'color': delta_v_line_color,
                       'linestyle': '-',
                       'ytitle': '$R_{\\rm w}$',
                       }
    theta_w_deg_dict_option = {'color': delta_v_line_color,
                               'linestyle': '-',
                               'yrange': [0, 180],
                               'ytitle': '$\\theta_{\\rm w}$',
                               'ysubtitle': 'deg',
                               }

    msp_time = df_mms.index[ind_range_msp[int(len(ind_range_msp)/2)]]
    msh_time = df_mms.index[ind_range_msh[int(len(ind_range_msh)/2)]]

    # Convert msp_time to UNIX time
    msp_time_unix = msp_time.timestamp()
    msh_time_unix = msh_time.timestamp()
    t_jet_center_unix = t_jet_center.timestamp()
    ptt.timebar(msp_time_unix, databar=False, color='red', dash=True, thick=2)
    ptt.timebar(msh_time_unix, databar=False, color='blue', dash=True, thick=2)
    ptt.timebar(t_jet_center_unix, databar=False, color='green', dash=True, thick=1)

    ptt.options(f'mms{probe}_dis_energyspectr_omni_{data_rate}',
                opt_dict=ion_energy_spectr_dict_option)
    ptt.options(f'mms{probe}_des_energyspectr_omni_{data_rate}',
                opt_dict=electron_energy_spectr_dict_option)
    ptt.options(f'mms{probe}_dis_numberdensity_{data_rate}',
                opt_dict=number_density_dict_option)

    ptt.options('Tp', opt_dict=tp_dict_option)
    ptt.options('mms3_fgm_b_lmn_srvy_l2', opt_dict=b_dict_option)
    ptt.options(f'mms{probe}_dis_bulkv_lmn_{data_rate}', opt_dict=bulkv_dict_option)
    ptt.options('delta_v_vp_lmn_diff_l', opt_dict=delta_v_dict_option)
    ptt.options('R_w', opt_dict=r_w_dict_option)
    ptt.options('theta_w_deg', opt_dict=theta_w_deg_dict_option)

    if (walen_v1 or walen_v2) & jet_detection:
        folder_name = f"../figures/jet_reversal_checks/check_{date_obs}/{data_rate}/jet_walen"
        # If the folder doesn't exist, create it
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
    elif (walen_v1 or walen_v2) & (not jet_detection):
        folder_name = f"../figures/jet_reversal_checks/check_{date_obs}/{data_rate}/walen"
        # If the folder doesn't exist, create it
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
    elif (not walen_v1) & (
            not walen_v2) & jet_detection:
        folder_name = f"../figures/jet_reversal_checks/check_{date_obs}/{data_rate}/jet"
        # If the folder doesn't exist, create it
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
    else:
        folder_name = f"../figures/jet_reversal_checks/check_{date_obs}/{data_rate}/no_jet_no_walen"
        # If the folder doesn't exist, create it
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

    try:
        r_w_median = np.round(np.nanmedian(ptt.get_data('R_w')[1]), 1)
        theta_w = ptt.get_data('theta_w_deg')[1]
    except Exception:
        r_w_median = np.nan
        theta_w = np.nan
    # If theta_w is greater than 90 degree, then subtract it from 180 degree
    try:
        theta_w[theta_w > 90] = 180 - theta_w[theta_w > 90]
        theta_w_median = int(np.nanmedian(theta_w))
    except Exception:
        theta_w_median = np.nan

    figname = f"{folder_name}/mms{probe}_{t_jet_center.strftime('%Y%m%d_%H%M')}_" + \
              f"{str(ind_crossing).zfill(5)}_{shear_val}s_{r_w_median}rw_{theta_w_median}th"
    # figname = 'test'
    ptt.tplot(keys_to_plot, save_png=figname, display=False)

    return None
