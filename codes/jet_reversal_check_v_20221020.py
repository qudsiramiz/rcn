import datetime
import os
from turtle import color

import matplotlib.pyplot as plt
from matplotlib.transforms import Transform
import more_itertools as mit

# from matplotlib import rc
import numpy as np
import pandas as pd
import pyspedas as spd
import pyspedas.mms.cotrans.mms_cotrans_lmn as mms_cotrans_lmn
import pytplot as ptt
import pytz


def jet_reversal_check(crossing_time=None, dt=90, probe=3, data_rate='fast', level='l2',
                       coord_type='lmn', data_type='dis-moms', time_clip=True, latest_version=True,
                       jet_len=5, figname='mms_jet_reversal_check',
                       fname='../data/mms_jet_reversal_times.csv',
                       error_file_log_name="../data/mms_jet_reversal_check_error_log.csv",
                       verbose=True
                       ):
    """
    For a given crossing time and a given probe, the function finds out if MMS observed a jet during
    magnetopause crossing. If there was indeed a jet reversal, then the function saves that time to
    a csv file, along with the probe number, position of the spacecraft, and the time of the
    crossing.

    Parameters
    ----------
    crossing_time : datetime.datetime
        The time of the crossing of the magnetopause.
    dt : int
        The time interval to look for a jet reversal.
    probe : int
        The probe number. Can be 1, 2, 3, or 4. Default is 3. (since the magnetopause crossing times
        are for MMS3)
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
    fname : str
        The name of the csv file to save the data to. Default is 'mms_jet_reversal_times.csv'.

    Returns
    -------
    df_mms_fpi : pandas.DataFrame
        The dataframe containing the data from the FPI.

    df_mms_fgm : pandas.DataFrame
        The dataframe containing the data from the FGM.

    df_mms : pandas.DataFrame
        The merged dataframe containing the data from the FPI and FGM.
    """

    # Define the mass of proton in kg
    m_p = 1.6726219e-27

    # Define the absolute permeability of free space in m^2 kg^-1 s^-1
    mu_0 = 4 * np.pi * 1e-7

    # Define the Boltzmann constant in J K^-1
    k_B = 1.38064852e-23

    # Define the conversion factor for converting temperature from ev to K
    ev_to_K = 11604.525

    crossing_time_min = crossing_time - datetime.timedelta(seconds=dt)
    crossing_time_max = crossing_time + datetime.timedelta(seconds=dt)
    trange = [crossing_time_min, crossing_time_max]

    # Get the data from the FPI
    mms_fpi_varnames = [f'mms{probe}_dis_numberdensity_{data_rate}',
                        f'mms{probe}_dis_bulkv_gse_{data_rate}',
                        f'mms{probe}_dis_temppara_{data_rate}',
                        f'mms{probe}_dis_tempperp_{data_rate}']

    _ = spd.mms.fpi(trange=trange, probe=probe, data_rate=data_rate, level=level,
                    datatype=data_type, varnames=mms_fpi_varnames, time_clip=time_clip,
                    latest_version=latest_version)

    mms_fpi_time = ptt.get_data(mms_fpi_varnames[0])[0]
    # Convert the time to a datetime object
    mms_fpi_time = pd.to_datetime(mms_fpi_time, unit='s')
    mms_fpi_time = mms_fpi_time.tz_localize(pytz.utc)

    mms_fpi_numberdensity = ptt.get_data(mms_fpi_varnames[0])[1]
    _ = ptt.get_data(mms_fpi_varnames[1])[1:4][0]
    mms_fpi_temppara = ptt.get_data(mms_fpi_varnames[2])[1]
    mms_fpi_tempperp = ptt.get_data(mms_fpi_varnames[3])[1]

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
                                            gsm=True, probe=str(probe), data_rate=data_rate)

        mms_fpi_bulkv_lmn = ptt.get_data(f'mms{probe}_dis_bulkv_lmn_{data_rate}')[1:4][0]

    # Create a dataframe with the FPI data
    if coord_type == 'lmn':
        df_mms_fpi = pd.DataFrame(index=mms_fpi_time, data={'np': mms_fpi_numberdensity,
                                                            'vp_lmn_n': mms_fpi_bulkv_lmn[:, 2],
                                                            'vp_lmn_m': mms_fpi_bulkv_lmn[:, 1],
                                                            'vp_lmn_l': mms_fpi_bulkv_lmn[:, 0],
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
        df_mms_fpi['vp_diff_x'] = df_mms_fpi['vp_lmn_n'] - df_mms_fpi['vp_lmn_n_rolling_median']
        df_mms_fpi['vp_diff_y'] = df_mms_fpi['vp_lmn_m'] - df_mms_fpi['vp_lmn_m_rolling_median']
        df_mms_fpi['vp_diff_z'] = df_mms_fpi['vp_lmn_l'] - df_mms_fpi['vp_lmn_l_rolling_median']
        # df_mms_fpi['vp_diff_z'] = df_mms_fpi['vp_lmn_l'] - np.nanmean(df_mms_fpi['vp_lmn_l'])
    else:
        df_mms_fpi['vp_diff_x'] = df_mms_fpi['vp_gsm_x'] - df_mms_fpi['vp_gsm_x_rolling_median']
        df_mms_fpi['vp_diff_y'] = df_mms_fpi['vp_gsm_y'] - df_mms_fpi['vp_gsm_y_rolling_median']
        df_mms_fpi['vp_diff_z'] = df_mms_fpi['vp_gsm_z'] - df_mms_fpi['vp_gsm_z_rolling_median']

    # If the absolute of maximum or minimum value of 'vp_gsm_z_diff' is greater than the threshold,
    # then check if we observed a jet
    v_thresh = 70  # Defined based on values in literature (Trattner et al. 2017)
    jet_detection = False
    # Find out the median value of difference between the indices of df_mms_fpi in seconds
    median_diff = np.median(np.diff(df_mms_fpi.index.astype(np.int64) / 10**9))
    print(f"Median difference between the indices of df_mms_fpi in seconds: {median_diff}")
    n_points = int(jet_len / median_diff)

    if np.abs(df_mms_fpi['vp_diff_z']).max() > v_thresh:

        # Find the index where the maximum value is
        ind_max_z_diff = df_mms_fpi.index[
            df_mms_fpi['vp_diff_z'] == df_mms_fpi['vp_diff_z'].max()]

        # Set a time window of +/- 60 seconds around the maximum value
        time_check_range = [ind_max_z_diff[0] - pd.Timedelta(minutes=1),
                            ind_max_z_diff[0] + pd.Timedelta(minutes=1)]

        # Define a velocity which might refer to a jet
        vp_jet = df_mms_fpi['vp_diff_z'].loc[time_check_range[0]:time_check_range[1]]

        # If there is a jet, then it must have sustained velocity above the threshold for at least
        # jet_len seconds. If not, then it is not a jet.
        # Based on the data rate, find out the number of data points in the time window which must
        # be checked for a jet
        # NOTE: 12 seconds comes from Trenchi et al. (2008) based on their study using Double Star
        # observation. Though they used it to confirm the Walen relation, I believe that similar
        # principal can be used for jet detection as well.
        if n_points == 0 or n_points > len(vp_jet):
            n_points = 2
            if verbose:
                print(f"Since time interval was greater than jet_len, setting it to {n_points}\n")
        # Find out the indices of all such points
        ind_right_vals = np.flatnonzero(np.convolve(vp_jet > v_thresh,
                                        np.ones(n_points, dtype=int), 'same') >= n_points)
        ind_left_vals = np.flatnonzero(np.convolve(vp_jet < -v_thresh,
                                       np.ones(n_points, dtype=int), 'same') >= n_points)

        if (len(ind_right_vals) and len(ind_left_vals)) > 0:
            jet_detection = True
            # Set the jet location to union of the positive and negative indices
            ind_jet = np.union1d(ind_right_vals, ind_left_vals)
            # Find time corresponding to the center of the jet
            jet_time = vp_jet.index[ind_jet[0]] + (
                       vp_jet.index[ind_jet[-1]] - vp_jet.index[ind_jet[0]]
                       ) / 2
            if verbose:
                print(f"\n\033[1;32m Jet detected at {jet_time} \033[0m \n")
        else:
            ind_jet = np.array([])
            if verbose:
                print("\033[1;31m No jet detected \033[0m \n")
    else:
        ind_jet = np.array([])
        if verbose:
            print("\033[1;31m No jet detected \033[0m \n")

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
        _ = mms_cotrans_lmn.mms_cotrans_lmn(name_in=f'mms{probe}_fgm_b_gsm_srvy_{level}',
                                            name_out=f'mms{probe}_fgm_b_lmn_srvy_{level}',
                                            gsm=True, probe=str(probe), data_rate=data_rate)

        mms_fgm_b_lmn = ptt.get_data(f'mms{probe}_fgm_b_lmn_srvy_{level}')[1:4][0]

    # Create a dataframe with the FGM data
    if coord_type == 'lmn':
        df_mms_fgm = pd.DataFrame(index=mms_fgm_time, data={'b_lmn_n': mms_fgm_b_lmn[:, 2],
                                                            'b_lmn_m': mms_fgm_b_lmn[:, 1],
                                                            'b_lmn_l': mms_fgm_b_lmn[:, 0]})
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
    # Make the time index timezone aware
    try:
        df_mms.index = df_mms.index.tz_localize(pytz.utc)
    except Exception:
        if verbose:
            print("\033[1;31m Timezone conversion failed \033[0m \n")

    # Standardise the proton density in the dataframe
    # df_mms['np_std'] = (df_mms['np'] - df_mms['np'].mean()) / df_mms['np'].std()

    # Find out the index for which 'np_std' is minimum
    # ind_min_np = np.where(df_mms['np_std'] == df_mms['np_std'].min())[0][0]
    if len(ind_jet) > 0:
        # Get the time check range centered around jet location
        time_check_center = vp_jet.index[ind_jet[0]] + (
                            vp_jet.index[ind_jet[-1]] - vp_jet.index[ind_jet[0]]) / 2
        time_check_range = [time_check_center - pd.Timedelta(60, unit='s'),
                            time_check_center + pd.Timedelta(60, unit='s')]

        # Find out where np is minimum in a time range around the jet location
        np_check = df_mms['np'].loc[time_check_range[0]:time_check_range[-1]]
        # Get the time corresponding to the minimum np
        np_min_time = np_check[np_check == np_check.min()].index[0]
        # Get the index of the minimum np time in that range
        ind_min_np = np.where(df_mms.index == np_min_time)[0][0]
        if verbose:
            print(f"ind_min_np: {ind_min_np}")
    else:
        ind_min_np = np.where(df_mms['np'] == df_mms['np'].min())[0][0]

    if verbose:
        print(f"Magnetopause location: {df_mms.index[ind_min_np]} at location index {ind_min_np}\n")

    # Find the index value closest to 'ind_min_np' where 'np' is greater than 'n_thresh' on both
    # sides of the minimum value
    # TODO: Check if threshold value of 3 is fine or if we need to decrease/increase it
    # n_thresh = 3
    n_thresh_msp = 5
    n_thresh_msh = 10

    # TODO: Instead of one difference between indices, use the median value of differences between
    # all the indices
    n_points_msp = int(5 / median_diff)
    n_points_msh = int(10 / median_diff)

    n_points_walen = n_points * 3

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

    # Set the index values to the full range where we have decided magnetosphere and magnetosheath
    # are.
    ind_msp = ind_range_msp
    ind_walen_check_min = np.where(df_mms.index == vp_jet.index[ind_jet[0]])[0][0]
    ind_walen_check_max = np.where(df_mms.index == vp_jet.index[ind_jet[-1]])[0][0]
    ind_walen_check = np.arange(ind_walen_check_min, ind_walen_check_max)
    ind_msh = ind_range_msh

    print(f"ind jet is {ind_jet}\n min value of walen is {ind_walen_check_min}\n max value of walen is {ind_walen_check_max}\n")
    # Get different parameters for magnetosphere and magnetosheath
    np_msp = df_mms['np'][ind_msp] * 1e6  # Convert to m^-3 from cm^-3
    np_msh = df_mms['np'][ind_msh] * 1e6  # Convert to m^-3 from cm^-3

    if coord_type == 'lmn':

        vp_lmn_vec_msp = np.array([df_mms['vp_lmn_n'][ind_msp], df_mms['vp_lmn_m'][ind_msp],
                                   df_mms['vp_lmn_l'][ind_msp]]) * 1e3  # Convert to m/s from km/s
        vp_lmn_vec_msp = vp_lmn_vec_msp.T

        vp_lmn_vec_msh = np.array([df_mms['vp_lmn_n'][ind_msh], df_mms['vp_lmn_m'][ind_msh],
                                   df_mms['vp_lmn_l'][ind_msh]]) * 1e3  # Convert to m/s from km/s
        vp_lmn_vec_msh = vp_lmn_vec_msh.T

        b_lmn_vec_msp = np.array([df_mms['b_lmn_n'][ind_msp], df_mms['b_lmn_m'][ind_msp],
                                  df_mms['b_lmn_l'][ind_msp]]) * 1e-9  # Convert to T from nT
        b_lmn_vec_msp = b_lmn_vec_msp.T
        b_lmn_vec_msh = np.array([df_mms['b_lmn_n'][ind_msh], df_mms['b_lmn_m'][ind_msh],
                                  df_mms['b_lmn_l'][ind_msh]]) * 1e-9  # Convert to T from nT
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
        vp_gse_vec_msp = np.array([df_mms['vp_gse_x'][ind_msp], df_mms['vp_gse_y'][ind_msp],
                                   df_mms['vp_gse_z'][ind_msp]]) * 1e3  # Convert to m/s from km/s
        vp_gse_vec_msp = vp_gse_vec_msp.T
        vp_gse_vec_msh = np.array([df_mms['vp_gse_x'][ind_msh], df_mms['vp_gse_y'][ind_msh],
                                   df_mms['vp_gse_z'][ind_msh]]) * 1e3  # Convert to m/s from km/s
        vp_gse_vec_msh = vp_gse_vec_msh.T
        b_gse_vec_msp = np.array([df_mms['b_gse_x'][ind_msp], df_mms['b_gse_y'][ind_msp],
                                  df_mms['b_gse_z'][ind_msp]]) * 1e-9  # Convert to T from nT
        b_gse_vec_msp = b_gse_vec_msp.T
        b_gse_vec_msh = np.array([df_mms['b_gse_x'][ind_msh], df_mms['b_gse_y'][ind_msh],
                                  df_mms['b_gse_z'][ind_msh]]) * 1e-9  # Convert to T from nT
        b_gse_vec_msh = b_gse_vec_msh.T

    tp_para_msp = df_mms['tp_para'][ind_msp] * ev_to_K  # Convert to K from ev
    tp_para_msh = df_mms['tp_para'][ind_msh] * ev_to_K  # Convert to K from ev
    tp_perp_msp = df_mms['tp_perp'][ind_msp] * ev_to_K  # Convert to K from ev
    tp_perp_msh = df_mms['tp_perp'][ind_msh] * ev_to_K  # Convert to K from ev

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

    method = 2
    if method == 1:
        alpha_jet = np.full(len(ind_walen_check), np.nan)
        alpha_msh = np.full(len(ind_msh), np.nan)
        v_th_jet = np.full((len(ind_walen_check), 3), np.nan)
        v_th_msh = np.full((len(ind_msh), 3), np.nan)

        if coord_type == 'lmn':
            for i, xx in enumerate(ind_walen_check):
                alpha_jet[i] = (mu_0 * df_mms['np'][xx] * 1e6 * k_B) * (
                                df_mms['tp_para'][xx] - df_mms['tp_perp'][xx]) * ev_to_K / (1e-18 * (
                                df_mms['b_lmn_n'][xx] ** 2 + df_mms['b_lmn_m'][xx] ** 2 +
                                df_mms['b_lmn_l'][xx] ** 2))

                b_lmn_vec_jet = np.array([df_mms['b_lmn_n'][xx], df_mms['b_lmn_m'][xx],
                                            df_mms['b_lmn_l'][xx]]) * 1e-9
                for j in range(3):
                    v_th_jet[i, j] = b_lmn_vec_jet[j] * ((1 - alpha_jet[i]) / (mu_0 * df_mms['np'][xx]
                                                          * 1e6 * m_p))**0.5
            for i in range(len(ind_msh)):
                alpha_msh[i] = (mu_0 * np_msh[i] * k_B) * (tp_para_msh[i] - tp_perp_msh[i]) / (
                    np.linalg.norm(b_lmn_vec_msh[i, :])**2)
                for j in range(3):
                    v_th_msh[i, j] = b_lmn_vec_msh[i, j] * ((1 - alpha_msh[i]) / (
                                                             mu_0 * np_msh[i] * m_p))**0.5
        else:
            for i in range(len(ind_jet)):
                alpha_jet[i] = (mu_0 * np_msp[i] * k_B) * (tp_para_msp[i] - tp_perp_msp[i]) / (
                    np.linalg.norm(b_gse_vec_msp[i, :])**2)
                alpha_msh[i] = (mu_0 * np_msh[i] * k_B) * (tp_para_msh[i] - tp_perp_msh[i]) / (
                    np.linalg.norm(b_gse_vec_msh[i, :])**2)
                for j in range(3):
                    v_th_jet[i, j] = b_gse_vec_msp[i, j] * ((1 - alpha_jet[i]) / (
                                            mu_0 * np_msp[i] * m_p))**0.5
                    v_th_msh[i, j] = b_gse_vec_msh[i, j] * (1 - alpha_msh[i]) / (
                        mu_0 * np_msh[i] * m_p * (1 - alpha_msh[i])
                    )**0.5

        alpha_all = mu_0 * df_mms['np'] * 1e6 * k_B * (df_mms['tp_para'] - df_mms['tp_perp']
                    ) * ev_to_K / (1e-18 * (df_mms['b_lmn_n'] ** 2 + df_mms['b_lmn_m'] ** 2
                    + df_mms['b_lmn_l'] ** 2))

        v_th_all = np.full((len(df_mms['tp_para']), 3), np.nan)
        for i in range(len(df_mms['tp_para'])):
            v_th_all[i, :] = 1e-9 * np.array([df_mms['b_lmn_n'][i], df_mms['b_lmn_m'][i],
                                              df_mms['b_lmn_l'][i]]) * ((1 - alpha_all[i]) / (
                                                                         mu_0 * df_mms['np'][i] * 1e6
                                                                         * m_p))**0.5
        delta_v_th = np.nanmean(v_th_msh, axis=0) - v_th_jet
        delta_v_th_all = np.nanmean(v_th_msh, axis=0) - v_th_all

        vp_lmn_vec_walen = np.array([df_mms['vp_lmn_n'][ind_walen_check],
                                     df_mms['vp_lmn_m'][ind_walen_check],
                                     df_mms['vp_lmn_l'][ind_walen_check]]).T * 1e3  # km/s to m/s
        vp_lmn_vec_all = np.array([df_mms['vp_lmn_n'], df_mms['vp_lmn_m'], df_mms['vp_lmn_l']]).T * 1e3

        if coord_type == 'lmn':
            delta_v_obs = vp_lmn_vec_msh_mean - vp_lmn_vec_walen
            delta_v_obs_all = vp_lmn_vec_msh_mean - vp_lmn_vec_all
            print(f"Shape of vp_lmn_vec_msh_mean: {vp_lmn_vec_msh_mean.shape}")
        else:
            delta_v_obs = vp_gse_vec_msh - vp_gse_vec_msp
        # delta_v_obs_mag = np.linalg.norm(delta_v_obs, axis=1)

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

        bool_array = np.logical_and(np.logical_and(0.4 < R_w, R_w < 3),
                                    np.logical_and(np.logical_and(0 < theta_w_deg, theta_w_deg < 30),
                                    np.logical_and(150 < theta_w_deg, theta_w_deg < 180)))

        # Check the number of True in the bool_array
        num_true = np.sum(bool_array)
        if num_true >= n_points_walen/2:
            walen_relation_satisfied = True
        else:
            walen_relation_satisfied = False

        ind_walen_vals = np.flatnonzero(np.convolve(bool_array > 0,
                                        np.ones(n_points_walen, dtype=int),
                                        'valid') >= n_points_walen/3)

        if len(ind_walen_vals) > 0:
            walen_relation_satisfied_v2 = True
        else:
            walen_relation_satisfied_v2 = False

        if walen_relation_satisfied:
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

    else:

        walen_relation_satisfied, walen_relation_satisfied_v2, theta_w_deg, R_w, theta_all_deg, R_all, ind_walen_vals = check_walen_relation(df_mms=df_mms, ind_walen_check=ind_walen_check, coord_type=coord_type, ind_msh=ind_msh, n_points_walen=n_points_walen, jet_detection=jet_detection)

    # Check if within 2 minutes of crossing time the values went above and below the threshold
    # If ind_vals is not empty, then append the crossing time to the csv file
    # if len(ind_vals) > 0:
    if walen_relation_satisfied or jet_detection or walen_relation_satisfied_v2:
        # for xxx in range(1):
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        # _ = spd.mms.mec
        # mms_sc_pos = ptt.get_data(mms_mec_varnames[0])[1:3][0][0] / r_e
        x = np.nanmean(mms_fgm_r_gsm[:, 0]) / r_e
        y = np.nanmean(mms_fgm_r_gsm[:, 1]) / r_e
        z = np.nanmean(mms_fgm_r_gsm[:, 2]) / r_e
        r_yz = np.sqrt(y**2 + z**2)  # Projection distance in yz plane.

        # TODO: Add magnetic beta for sheath
        # List of variables to be saved in the csv file
        var_list = 'Date,Probe,walen1,walen2,jet_detection,x_gsm,y_gsm,z_gsm,r_spc,r_W,theta_w,'\
                   'jet_time,ind_min_msp,ind_max_msp,ind_min_msh,ind_max_msh,'\
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
                     'walen1': walen_relation_satisfied,
                     'walen2': walen_relation_satisfied_v2,
                     'jet_detection': jet_detection,
                     'x_gsm': np.round(x, 3),
                     'y_gsm': np.round(y, 3),
                     'z_gsm': np.round(z, 3),
                     'r_spc': np.round(r_yz, 3),
                     'r_W': np.round(np.nanmedian(R_w), 3),
                     'theta_w': np.round(np.nanmedian(theta_w_deg), 3),
                     'jet_time': jet_time,
                     'ind_min_msp': ind_min_msp,
                     'ind_max_msp': ind_max_msp,
                     'ind_min_msh': ind_min_msh,
                     'ind_max_msh': ind_max_msh,
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

    # Get the index corresponding to the crossing time in the data
    try:
        df_crossing_temp = pd.read_csv("../data/mms_jet_reversal_times_list_20221017_beta.csv",
                                       index_col=False)
        # Change the formatting of the date
        for xx in range(len(df_crossing_temp)):
            df_crossing_temp['jet_time'][xx] = df_crossing_temp['jet_time'][xx].split('+')[0].split('.')[0]
        df_crossing_temp.set_index("jet_time", inplace=True)
        crossing_time_str = crossing_time.strftime("%Y-%m-%d %H:%M:%S")[:]
    except:
        df_crossing_temp = pd.read_csv("../data/mms_magnetopause_crossings.csv", index_col=False)
        df_crossing_temp.set_index("DateStart", inplace=True)
        crossing_time_str = crossing_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]

    ind_crossing = np.where(df_crossing_temp.index == crossing_time_str)[0][0]

    # Set the fontstyle to Times New Roman
    font = {'family': 'serif', 'weight': 'normal', 'size': 16 }
    plt.rc('font', **font)
    plt.rc('text', usetex=True)
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

    # plt.style.use('dark_background')

    plt.close("all")

    if coord_type == 'lmn':
        # Make a figure with the B, np, and vp in L, M, N coordinates
        lw = 0.75
        # Set the tickmark location inside the plot for all subplots
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        # TODO: Consider adding time series of R_w and Theta_w to the plot
        fig, axs = plt.subplots(5, 1, figsize=(10, 12), sharex=True)
        fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0, hspace=0.0)
        axs[0].plot(df_mms.index, df_mms['b_lmn_l'].rolling(1).mean(), label=r'$B_L$', color='b',
                    lw=lw)
        axs[0].plot(df_mms.index, df_mms['b_lmn_m'].rolling(1).mean(), label=r'$B_M$', color='g',
                    lw=lw)
        axs[0].plot(df_mms.index, df_mms['b_lmn_n'].rolling(1).mean(), label=r'$B_N$', color='r',
                    lw=lw)
        axs[0].set_ylabel(r'$B$ (nT)')
        # axs[0].legend(loc='upper right')
        # Write a text box with B_l and B_m and B_n in different colors rotated by 0 degrees
        axs[0].text(1.01, 0.80, r'$B_{\rm L}$', color='b', transform=axs[0].transAxes, rotation=0, 
                    ha='left', va='center')
        axs[0].text(1.01, 0.50, r'$B_{\rm M}$', color='g', transform=axs[0].transAxes, rotation=0, 
                    ha='left', va='center')
        axs[0].text(1.01, 0.2, r'$B_{\rm N}$', color='r', transform=axs[0].transAxes, rotation=0, 
                    ha='left', va='center')

        # axs[0].grid(True)

        # MMS Density plots
        axs[1].plot(df_mms.index, df_mms['np'].rolling(1).mean(), color='b', lw=lw)
        axs[1].set_ylabel(r'$n_{\rm p}$ (cm$^{-3}$)', color='b')
        axs[1].spines['left'].set_color('blue')
        axs[1].tick_params(which="both", axis='y', colors='b')
        # axs[1].grid(True)
        axs[1].set_yscale('log')

        # MMS temperature plots
        axs_t = axs[1].twinx()
        axs_t.plot(df_mms.index, df_mms['tp_para'].rolling(1).mean(), color='r', lw=lw)
        axs_t.plot(df_mms.index, df_mms['tp_perp'].rolling(1).mean(), color='g', lw=lw)
        # axs_t.plot(df_mms.index, (2 * df_mms['tp_para'] + df_mms['tp_perp']) / 3, label='T', 
        #            color='k', lw=lw)
        axs_t.set_ylabel(r'$T_{\rm p}$ (K)', color='r')
        # Write a text box with T_para and T_perp in different colors rotated by 0 degrees
        axs_t.text(0.99, 0.55, r'$T_{\rm p \parallel}$', color='r', transform=axs_t.transAxes,
                   rotation=0, ha='right', va='bottom')
        axs_t.text(0.99, 0.45, r'$T_{\rm p \perp}$', color='g', transform=axs_t.transAxes,
                   rotation=0, ha='right', va='top')
        axs_t.spines['right'].set_color('red')
        axs_t.spines['left'].set_visible(False)
        axs_t.tick_params(which='both', axis='y', colors='r')
        # axs_t.legend(loc='upper right')

        # axs_t.grid(True)
        axs_t.set_yscale('log')

        axs[1].axvspan(df_mms.index[ind_msp[0]], df_mms.index[ind_msp[-1]], color="r",
                       alpha=0.2, label="Magnetosphere")
        axs[1].axvspan(df_mms.index[ind_msh[0]], df_mms.index[ind_msh[-1]], color="b",
                       alpha=0.2, label="Magnetosheath")
        # axs[1].legend(loc=2)
        axs[1].text(0.01, 0.55, 'Magnetosphere', color='r', transform=axs[1].transAxes, ha='left',
                    va='bottom')
        axs[1].text(0.01, 0.45, 'Magnetosheath', color='b', transform=axs[1].transAxes, ha='left',
                    va='top')

        axs[2].plot(df_mms.index, df_mms['vp_lmn_l'].rolling(1).mean(), label=r'$V_L$', color='b',
                    lw=lw)
        axs[2].plot(df_mms.index, df_mms['vp_lmn_m'].rolling(1).mean(), label=r'$V_M$', color='g',
                    lw=lw)
        axs[2].plot(df_mms.index, df_mms['vp_lmn_n'].rolling(1).mean(), label=r'$V_N$', color='r',
                    lw=lw)
        axs[2].set_ylabel(r'$V_{\rm p}$ (km/s)')
        # Wrtie a text box with V_L, V_M and V_N in different colors rotated by 0 degrees
        axs[2].text(1.01, 0.80, r'$V_{\rm L}$', color='b', transform=axs[2].transAxes, rotation=0,
                    ha='left', va='center')
        axs[2].text(1.01, 0.50, r'$V_{\rm M}$', color='g', transform=axs[2].transAxes, rotation=0,
                    ha='left', va='center')
        axs[2].text(1.01, 0.2, r'$V_{\rm N}$', color='r', transform=axs[2].transAxes, rotation=0,
                    ha='left', va='center')
        # axs[2].legend(loc='upper right')
        # axs[2].grid(True)

        if walen_relation_satisfied:
            axs[3].plot(df_mms.index, df_mms.vp_diff_z.rolling(1).mean(), 'g-', lw=lw)
        elif walen_relation_satisfied_v2:
            axs[3].plot(df_mms.index, df_mms.vp_diff_z.rolling(1).mean(), 'b-', lw=lw)
        else:
            axs[3].plot(df_mms.index, df_mms.vp_diff_z.rolling(1).mean(), 'r-', lw=lw)
        if jet_detection:
            # Make a box around the jet
            axs[3].axvspan(vp_jet.index[ind_jet[0]], vp_jet.index[ind_jet[-1]], color='c',
                           alpha=0.2)
            axs[4].axvspan(vp_jet.index[ind_jet[0]], vp_jet.index[ind_jet[-1]], color='c',
                           alpha=0.2, label='Jet Location')
            # Plot a vertical line at the jet time center
            jet_center = vp_jet.index[ind_jet[0]] + (vp_jet.index[ind_jet[-1]] -
                                                     vp_jet.index[ind_jet[0]]) / 2
            axs[3].axvline(jet_center, color='k', lw=0.5, label='Jet Center')
            axs[3].legend(loc=1)
            axs[3].text(0.01, 0.8, 'Jet Location', color='c', transform=axs[3].transAxes,
                        ha='left', va='center')
        # plt.plot(vp_jet, 'r.', ms=2, lw=1, label="jet")
        # Draw a dashed line at +/- v_thres
        axs[3].axhline(y=v_thresh, color='k', linestyle='--', lw=lw)
        axs[3].axhline(y=-v_thresh, color='k', linestyle='--', lw=lw)
        axs[3].set_ylabel(r'$\Delta V_{\rm L}$ (km/s)')
        axs[3].set_ylim(-300, 300)

        df_mms["theta_all_deg"] = theta_all_deg
        df_mms['r_all'] = R_all

        # Get a 10 seconds rolling mean of theta_all_deg
        theta_all_deg_rolling = df_mms["theta_all_deg"].rolling(1).mean()
        # Get a 10 seconds rolling mean of r_all
        r_all_rolling = df_mms['r_all'].rolling(1).mean()

        # Make the theta_all_deg and R_all plot
        axs[4].plot(df_mms.index, theta_all_deg_rolling, 'b-', lw=lw, label=r'$\theta$', alpha=0.75)
        # axs[4].plot(df_mms.index[ind_walen_check], theta_w_deg, 'k-', lw=1.5*lw, alpha=1)
        axs[4].set_ylabel(r'$\theta_{\rm W}$ (deg)', color='b')
        axs[4].tick_params(which='both', axis='y', labelcolor='b')
        # Draw a horizontal line at theta = 0, 30, 150 and 180
        axs[4].axhline(y=30, color='b', linestyle='--', lw=lw)
        axs[4].axhline(y=150, color='b', linestyle='--', lw=lw)
        axs[4].set_ylim(0, 180)

        axs_r = axs[4].twinx()
        axs_r.plot(df_mms.index, r_all_rolling, 'r-', lw=lw, label=r'$R_{\rm W}$', alpha=0.75)
        # axs_r.plot(df_mms.index[ind_walen_check], R_w, '-', color='k', alpha=0)
        axs_r.set_ylabel(r'$R_{\rm W}$', color='r')
        # Draw a horizontal line at R = 0.4 and 3
        axs_r.axhline(y=0.4, color='r', linestyle='--', lw=lw)
        axs_r.axhline(y=3, color='r', linestyle='--', lw=lw)
        # Set the color of y-axis spine to red
        axs_r.spines['right'].set_color('red')
        # Set the color of the tick labels to red
        axs_r.tick_params(which='both', axis='y', colors='red')
        axs_r.set_yscale('log')


        axs[4].set_xlabel('Time (UTC)')
        # axs[4].set_xlim(df_mms.index[0], df_mms.index[-1])

        # Set the x-axis limits to 2 minutes before and after the jet
        # axs[2].set_xlim(jet_center - pd.Timedelta(minutes=2), jet_center +
        # pd.Timedelta(minutes=2))
        axs[4].set_xlim(df_mms.index.min(), df_mms.index.max())

        # Make a shaded region for all three plots where the magnetopause and magnetosheath
        # are
        axs[0].axvspan(df_mms.index[ind_msp[0]], df_mms.index[ind_msp[-1]], color="r",
                       alpha=0.2, label="Magnetosphere")
        axs[0].axvspan(df_mms.index[ind_msh[0]], df_mms.index[ind_msh[-1]], color="b",
                       alpha=0.2, label="Magnetosheath")

        axs[1].axvspan(df_mms.index[ind_msp[0]], df_mms.index[ind_msp[-1]], color="r",
                       alpha=0.2, label="Magnetosphere")
        axs[1].axvspan(df_mms.index[ind_msh[0]], df_mms.index[ind_msh[-1]], color="b",
                       alpha=0.2, label="Magnetosheath")
        axs[2].axvspan(df_mms.index[ind_msp[0]], df_mms.index[ind_msp[-1]], color="r",
                       alpha=0.2, label="Magnetosphere")
        axs[2].axvspan(df_mms.index[ind_msh[0]], df_mms.index[ind_msh[-1]], color="b",
                       alpha=0.2, label="Magnetosheath")
        axs[4].axvspan(df_mms.index[ind_msp[0]], df_mms.index[ind_msp[-1]], color="r",
                       alpha=0.2, label="Magnetosphere")
        axs[4].axvspan(df_mms.index[ind_msh[0]], df_mms.index[ind_msh[-1]], color="b",
                       alpha=0.2, label="Magnetosheath")

        temp3 = crossing_time.strftime('%Y-%m-%d %H:%M:%S')
        axs[0].set_title(f"MMS {probe} Jet Reversal Check at {temp3}")

        # Add text to the plot
        if len(ind_walen_vals) > 0:
            theta_w_deg_new = np.where(theta_w_deg[ind_walen_vals] > 90,
                                       180 - theta_w_deg[ind_walen_vals],
                                       theta_w_deg[ind_walen_vals])
            R_w_jet = R_w[ind_walen_vals]
        else:
            # For all the points in theta_w_deg, if it is greater than 90 then subtract 180 from it
            theta_w_deg_new = np.where(theta_w_deg > 90, 180 - theta_w_deg, theta_w_deg)
            R_w_jet = R_w

        text = f"$R_w$ = {np.nanmedian(R_w_jet):.2f}\n" +\
               f"$\\Theta_w$ = {np.nanmedian(theta_w_deg_new):.2f}"
        axs[4].text(0.02, 1, text, transform=axs[4].transAxes, ha='left', va='top')

        # Add the index of crossing time to the top of the first axis
        axs[0].text(-0.05, 1.07, f"{ind_crossing}", transform=axs[0].transAxes, ha='left',
                    va='top', color='r')

        axs[0].text(
            1.0, 1.01, f"$\\theta_{{B_{{msh}},B_{{msp}}}}=${angle_b_lmn_vec_msp_msh_median:.0f}",
            transform=axs[0].transAxes, ha='center',  va='bottom', color='r')
        if (walen_relation_satisfied or walen_relation_satisfied_v2) & jet_detection:
            folder_name = "../figures/jet_reversal_checks/check_20221020/jet_walen"
            # If the folder doesn't exist, create it
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)
        elif (walen_relation_satisfied or walen_relation_satisfied_v2) & (not jet_detection):
            folder_name = "../figures/jet_reversal_checks/check_20221020/walen"
            # If the folder doesn't exist, create it
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)
        elif (not walen_relation_satisfied) & (
                not walen_relation_satisfied_v2) & jet_detection:
            folder_name = "../figures/jet_reversal_checks/check_20221020/jet"
            # If the folder doesn't exist, create it
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)
        else:
            folder_name = "../figures/jet_reversal_checks/check_20221020/no_jet_no_walen"
            # If the folder doesn't exist, create it
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)

        ttt = str(crossing_time.strftime('%Y%m%d_%H%M'))
        fig_name = f"{folder_name}/mms{probe}_jet_reversal_{ttt}.png"
        plt.savefig(f"{fig_name}", dpi=150, bbox_inches='tight', pad_inches=0.1)
        # plt.savefig(f"{fig_name.replace('.png', '.pdf')}", dpi=300, bbox_inches='tight',
        #             pad_inches=0.1)
        print(f"{fig_name}")
        plt.close("all")

    else:
        ptt.tplot([f'mms{probe}_fgm_b_gsm_srvy_l2_bvec',
                   f'mms{probe}_dis_numberdensity_{data_rate}',
                   f'mms{probe}_dis_bulkv_gsm_{data_rate}'],
                  combine_axes=True, save_png='b_n_v_fig', display=False)
        plt.close("all")

    return df_mms_fpi, df_mms_fgm, df_mms


def check_walen_relation(df_mms=None, ind_walen_check=None, coord_type="gsm", ind_msh=None,
                         n_points_walen=None, jet_detection=False,):

    
    # Define the constants
    mu_0 = 4 * np.pi * 1e-7  # Permeability of free space (H/m)
    k_B = 1.38064852e-23  # Boltzmann constant in J/K
    m_p = 1.6726219e-27  # Mass of proton in kg
    ev_to_K = 11604.505  # Conversion factor from eV to K

    alpha_jet = np.full(len(ind_walen_check), np.nan)
    alpha_msh = np.full(len(ind_msh), np.nan)
    v_th_jet = np.full((len(ind_walen_check), 3), np.nan)
    v_th_msh = np.full((len(ind_msh), 3), np.nan)

    if coord_type == 'lmn':
        for i, xx in enumerate(ind_walen_check):
            alpha_jet[i] = (mu_0 * df_mms['np'][xx] * 1e6 * k_B) * (
                            df_mms['tp_para'][xx] - df_mms['tp_perp'][xx]) * ev_to_K / (1e-18 * (
                            df_mms['b_lmn_n'][xx] ** 2 + df_mms['b_lmn_m'][xx] ** 2 +
                            df_mms['b_lmn_l'][xx] ** 2))

            b_lmn_vec_jet = np.array([df_mms['b_lmn_n'][xx], df_mms['b_lmn_m'][xx],
                                        df_mms['b_lmn_l'][xx]]) * 1e-9
            for j in range(3):
                v_th_jet[i, j] = b_lmn_vec_jet[j] * ((1 - alpha_jet[i]) / (mu_0 * df_mms['np'][xx]
                                                      * 1e6 * m_p))**0.5
        for i,xx in enumerate(ind_msh):
            alpha_msh[i] = (mu_0 * df_mms['np'][xx] * 1e6 * k_B) * (
                            df_mms['tp_para'][xx] - df_mms['tp_perp'][xx]) * ev_to_K / (1e-18 * (
                            df_mms['b_lmn_n'][xx] ** 2 + df_mms['b_lmn_m'][xx] ** 2 +
                            df_mms['b_lmn_l'][xx] ** 2))

            b_lmn_vec_msh = np.array([df_mms['b_lmn_n'][xx], df_mms['b_lmn_m'][xx],
                                        df_mms['b_lmn_l'][xx]]) * 1e-9
            for j in range(3):
                v_th_msh[i, j] = b_lmn_vec_msh[j] * ((1 - alpha_msh[i]) / (mu_0 * df_mms['np'][xx]
                                                      * 1e6 * m_p))**0.5

    alpha_all = mu_0 * df_mms['np'] * 1e6 * k_B * (df_mms['tp_para'] - df_mms['tp_perp']
                ) * ev_to_K / (1e-18 * (df_mms['b_lmn_n'] ** 2 + df_mms['b_lmn_m'] ** 2
                + df_mms['b_lmn_l'] ** 2))

    v_th_all = np.full((len(df_mms['tp_para']), 3), np.nan)
    for i in range(len(df_mms['tp_para'])):
        v_th_all[i, :] = 1e-9 * np.array([df_mms['b_lmn_n'][i], df_mms['b_lmn_m'][i],
                                          df_mms['b_lmn_l'][i]]) * ((1 - alpha_all[i]) / (
                                                                     mu_0 * df_mms['np'][i] * 1e6
                                                                     * m_p))**0.5
    delta_v_th = np.nanmean(v_th_msh, axis=0) - v_th_jet
    delta_v_th_all = np.nanmean(v_th_msh, axis=0) - v_th_all

    vp_lmn_vec_walen = np.array([df_mms['vp_lmn_n'][ind_walen_check],
                                 df_mms['vp_lmn_m'][ind_walen_check],
                                 df_mms['vp_lmn_l'][ind_walen_check]]).T * 1e3  # km/s to m/s
    vp_lmn_vec_all = np.array([df_mms['vp_lmn_n'], df_mms['vp_lmn_m'], df_mms['vp_lmn_l']]).T * 1e3

    if coord_type == 'lmn':
        # Find the proton bulk velocity in the msh region and then compute its average
        vp_lmn_vec_msh = np.array([df_mms['vp_lmn_n'][ind_msh], df_mms['vp_lmn_m'][ind_msh],
                                   df_mms['vp_lmn_l'][ind_msh]]) * 1e3  # Convert to m/s from km/s
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

    bool_array = np.logical_and(np.logical_and(0.4 < R_w, R_w < 3),
                                np.logical_and(np.logical_and(0 < theta_w_deg, theta_w_deg < 30),
                                np.logical_and(150 < theta_w_deg, theta_w_deg < 180)))

    # Check the number of True in the bool_array
    num_true = np.sum(bool_array)
    if num_true >= n_points_walen/2:
        walen_relation_satisfied = True
    else:
        walen_relation_satisfied = False

    ind_walen_vals = np.flatnonzero(np.convolve(bool_array > 0,
                                    np.ones(n_points_walen, dtype=int),
                                    'valid') >= n_points_walen/3)

    if len(ind_walen_vals) > 0:
        walen_relation_satisfied_v2 = True
    else:
        walen_relation_satisfied_v2 = False

    if walen_relation_satisfied:
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

    return walen_relation_satisfied, walen_relation_satisfied_v2, theta_w_deg, R_w, theta_all_deg, R_all, ind_walen_vals