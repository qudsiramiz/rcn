import datetime
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyspedas as spd
import pytplot as ptt
import pytz


def jet_reversal_check(crossing_time=None, dt=90, probe=3, data_rate='fast', level='l2',
                       data_type='dis-moms', time_clip=True, latest_version=True, jet_len=5,
                       figname='mms_jet_reversal_check',
                       fname='../data/mms_jet_reversal_times.csv', verbose=True
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

    # Get the data from the FGM
    # mms_fgm_varnames = [f'mms{probe}_fgm_b_gsm_srvy_l2_bvec']
    # _ = spd.mms.fgm(trange=trange, probe=probe, time_clip=time_clip, latest_version=True,
    #                 varnames=[f"mms{probe}_fgm_b_gsm_srvy_{level}",
    #                           f"mms{probe}_fgm_r_gsm_srvy_{level}"], get_fgm_ephemeris=True)
    _ = spd.mms.fgm(trange=trange, probe=probe, time_clip=time_clip, latest_version=True,
                    get_fgm_ephemeris=True)
    # Get the time corresponding to the FGM data
    mms_fgm_time = ptt.get_data(f"mms{probe}_fgm_b_gsm_srvy_{level}")[0]
    mms_fgm_time = pd.to_datetime(mms_fgm_time, unit='s')
    mms_fgm_time = mms_fgm_time.tz_localize(pytz.utc)

    mms_fgm_b_gsm = ptt.get_data(f'mms{probe}_fgm_b_gsm_srvy_{level}')[1:4][0]
    mms_fgm_b_gse = ptt.get_data(f'mms{probe}_fgm_b_gse_srvy_{level}')[1:4][0]
    mms_fgm_r_gsm = ptt.get_data(f'mms{probe}_fgm_r_gsm_srvy_{level}')[1:4][0]

    # Create a dataframe with the FPI data
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
        print(f"\n\033[1;32m FPI dataframe created \033[0m \n")
        print(f"The fpi Datafram:\n {df_mms_fpi.head()}")

    # Add rolling mean to the dataframe
    df_mms_fpi['np_rolling_mean'] = df_mms_fpi['np'].rolling('60s', center=True).mean()
    df_mms_fpi['vp_gsm_x_rolling_mean'] = df_mms_fpi['vp_gsm_x'].rolling('60s', center=True).mean()
    df_mms_fpi['vp_gsm_y_rolling_mean'] = df_mms_fpi['vp_gsm_y'].rolling('60s', center=True).mean()
    df_mms_fpi['vp_gsm_z_rolling_mean'] = df_mms_fpi['vp_gsm_z'].rolling('60s', center=True).mean()
    df_mms_fpi['tp_para_rolling_mean'] = df_mms_fpi['tp_para'].rolling('60s', center=True).mean()
    df_mms_fpi['tp_perp_rolling_mean'] = df_mms_fpi['tp_perp'].rolling('60s', center=True).mean()

    # Find the difference wrt the rolling mean
    df_mms_fpi['np_diff'] = df_mms_fpi['np'] - df_mms_fpi['np_rolling_mean']
    df_mms_fpi['vp_gsm_x_diff'] = df_mms_fpi['vp_gsm_x'] - df_mms_fpi['vp_gsm_x_rolling_mean']
    df_mms_fpi['vp_gsm_y_diff'] = df_mms_fpi['vp_gsm_y'] - df_mms_fpi['vp_gsm_y_rolling_mean']
    df_mms_fpi['vp_gsm_z_diff'] = df_mms_fpi['vp_gsm_z'] - df_mms_fpi['vp_gsm_z_rolling_mean']

    # If the absolute of maximum or minimum value of 'vp_gsm_z_diff' is greater than the threshold,
    # then check if we observed a jet
    v_thresh = 70  # Defined based on values in literature (Trattner et al. 2017)
    jet_detection = False
    n_points = int(jet_len / (df_mms_fpi.index[1] - df_mms_fpi.index[0]).total_seconds())
    if n_points == 0:
            n_points = 3
            if verbose:
                print(f"Since time interval was greater than jet_len, setting it to {n_points}\n")
    if np.abs(df_mms_fpi['vp_gsm_z_diff']).max() > v_thresh:

        # Find the index where the maximum value is
        ind_max_z_diff = df_mms_fpi.index[
            df_mms_fpi['vp_gsm_z_diff'] == df_mms_fpi['vp_gsm_z_diff'].max()]

        # Set a time window of +/- 60 seconds around the maximum value
        time_check_range = [ind_max_z_diff[0] - pd.Timedelta(minutes=1),
                            ind_max_z_diff[0] + pd.Timedelta(minutes=1)]

        # Define a velocity which might refer to a jet
        vp_jet = df_mms_fpi['vp_gsm_z_diff'].loc[time_check_range[0]:time_check_range[1]]

        # If there is a jet, then it must have sustained velocity above the threshold for at least
        # jet_len seconds. If not, then it is not a jet.
        # Based on the data rate, find out the number of data points in the time window which must
        # be checked for a jet
        # NOTE: 12 seconds comes from Trenchi et al. (2008) based on their study using Double Star
        # observation. Though they used it to confirm the Walen relation, I believe that similar
        # principal can be used for jet detection as well.
        n_points = int(jet_len / (df_mms_fpi.index[1] - df_mms_fpi.index[0]).total_seconds())
        if n_points == 0:
            n_points = 3
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
        else:
            ind_jet = np.array([])
            if verbose:
                print("\033[1;31m No jet detected \033[0m \n")

    # Create a dataframe with the FGM data
    df_mms_fgm = pd.DataFrame(index=mms_fgm_time, data={'b_gsm_x': mms_fgm_b_gsm[:, 0],
                                                        'b_gsm_y': mms_fgm_b_gsm[:, 1],
                                                        'b_gsm_z': mms_fgm_b_gsm[:, 2],
                                                        'b_gse_x': mms_fgm_b_gse[:, 0],
                                                        'b_gse_y': mms_fgm_b_gse[:, 1],
                                                        'b_gse_z': mms_fgm_b_gse[:, 2]})

    # Make sure all time indices are in increasing order
    df_mms_fgm = df_mms_fgm.sort_index()

    if verbose:
        print(f"\n\033[1;32m FGM dataframe created \033[0m \n")
        print(f"FGM Dataframe: \n{df_mms_fgm.head()} \n")
    # Merge the two dataframes
    df_mms = pd.merge_asof(df_mms_fpi, df_mms_fgm, left_index=True, right_index=True)
    # Make the time index timezone aware
    df_mms.index = df_mms.index.tz_localize(pytz.utc)

    # Standardise the proton density in the dataframe
    df_mms['np_std'] = (df_mms['np'] - df_mms['np'].mean()) / df_mms['np'].std()

    # Find out the index for which 'np_std' is minimum
    ind_min_np_std = np.where(df_mms['np_std'] == df_mms['np_std'].min())[0][0]

    if verbose:
        print(f"Magnetopause location: {df_mms.index[ind_min_np_std]} at location index {ind_min_np_std}\n")

    # Find the index value closest to 'ind_min_np_std' where 'np_std' is greater than 'n_thresh' on
    # both sides of the minimum value
    # TODO: Check if 1 threshold is fine or if we need to decrease/increase it
    # FIXME: When the location of minimum value of 'np_std' is at the very start or end of the
    # dataframe, then ind_min_np_std_gt_05_left or ind_min_np_std_gt_05_right tends to be empty
    # since the where command returns an empty tuple. Think of a way to fix this.
    n_thresh = 1
    # try:
    #     ind_min_np_std_gt_05_right = np.where(
    #         df_mms['np_std'][ind_min_np_std:] > n_thresh)[0][0] + ind_min_np_std
    # except Exception:
    #     ind_min_np_std_gt_05_right = np.nan
    # try:
    #     ind_min_np_std_gt_05_left = np.where(df_mms['np_std'][:ind_min_np_std] > n_thresh)[0][-1]
    # except Exception:
    #     ind_min_np_std_gt_05_left = np.nan

    ind_min_np_std_gt_05_right = np.where(
             df_mms['np_std'][ind_min_np_std:] > n_thresh)[0][0] + ind_min_np_std

    ind_min_np_std_gt_05_left = np.where(df_mms['np_std'][:ind_min_np_std] > n_thresh)[0][-1]
    # Find the distance of the left and right hanfd indices in terms of number of indices
    diff_right = ind_min_np_std_gt_05_right - ind_min_np_std
    diff_left = ind_min_np_std - ind_min_np_std_gt_05_left

    n_points_walen = n_points * 3
    # Set the index value of magnetosheath to whichever one is closer to the minimum value.
    # 'ind_msp' is the index value of magnetosphere and 'ind_msh' is the index value of
    # magnetosheath.
    if diff_right < diff_left:
        ind_msp = ind_min_np_std
        ind_msh = ind_min_np_std_gt_05_right
        ind_range_msp = np.arange(max(0, ind_msp - n_points_walen), ind_msp)
        ind_range_msh = np.arange(ind_msh, min(len(df_mms), ind_msh + n_points_walen))
        # Compare the length of the two ranges. If one is larger than the other, then starting at
        # their starting point, set their length to be the same.
        if len(ind_range_msp) > len(ind_range_msh):
            ind_range_msp = ind_range_msp[(len(ind_range_msp) - len(ind_range_msh)):]
        elif len(ind_range_msp) < len(ind_range_msh):
            ind_range_msh = ind_range_msh[:len(ind_range_msp)]
    else:
        ind_msp = ind_min_np_std
        ind_msh = ind_min_np_std_gt_05_left
        ind_range_msp = np.arange(ind_msp, min(len(df_mms), ind_msp + n_points_walen))
        ind_range_msh = np.arange(max(0, ind_msh - n_points_walen), ind_msh)
        # Compare the length of the two ranges. If one is larger than the other, then starting at
        # their starting point, set their length to be the same.
        if len(ind_range_msp) > len(ind_range_msh):
            ind_range_msh = ind_range_msh[(len(ind_range_msp) - len(ind_range_msh)):]
        elif len(ind_range_msp) < len(ind_range_msh):
            ind_range_msp = ind_range_msp[:len(ind_range_msh)]

    # Set the index values to the full range where we have decided magnetosphere and magnetosheath
    # are.
    ind_msp = ind_range_msp
    ind_msh = ind_range_msh

    # Get different parameters for magnetosphere and magnetosheath
    np_msp = df_mms['np'][ind_msp] * 1e6  # Convert to m^-3 from cm^-3
    np_msh = df_mms['np'][ind_msh] * 1e6  # Convert to m^-3 from cm^-3
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
    tp_para_msp = df_mms['tp_para'][ind_msp] * 1160  # Convert to K from ev
    tp_para_msh = df_mms['tp_para'][ind_msh] * 1160  # Convert to K from ev
    tp_perp_msp = df_mms['tp_perp'][ind_msp] * 1160  # Convert to K from ev
    tp_perp_msh = df_mms['tp_perp'][ind_msh] * 1160  # Convert to K from ev

    # Define the mass of proton in kg
    m_p = 1.6726219e-27

    # Define the absolute permeability of free space in m^2 kg^-1 s^-1
    mu_0 = 4 * np.pi * 1e-7

    # Define the Boltzmann constant in J K^-1
    k_B = 1.38064852e-23

    alpha_msp = np.full(len(ind_msp), np.nan)
    alpha_msh = np.full(len(ind_msp), np.nan)
    v_th_msp = np.full((len(ind_msp), 3), np.nan)
    v_th_msh = np.full((len(ind_msp), 3), np.nan)

    for i in range(len(ind_msp)):
        alpha_msp[i] = (mu_0 * np_msp[i] * k_B) * (tp_para_msp[i] - tp_perp_msp[i]) / (
            np.linalg.norm(b_gse_vec_msp[i,:])**2)
        alpha_msh[i] = (mu_0 * np_msh[i] * k_B) * (tp_para_msh[i] - tp_perp_msh[i]) / (
            np.linalg.norm(b_gse_vec_msh[i,:])**2)
        for j in range(3):
            v_th_msp[i, j] = b_gse_vec_msp[i,j] * (1 - alpha_msp[i]) / (
                mu_0 * np_msp[i] * m_p * (1 - alpha_msp[i])
            )**0.5
            v_th_msh[i, j] = b_gse_vec_msh[i, j] * (1 - alpha_msh[i]) / (
                mu_0 * np_msh[i] * m_p * (1 - alpha_msh[i])
            )**0.5

    delta_v_th = v_th_msh - v_th_msp
    # delta_v_th_mag = np.linalg.norm(delta_v_th, axis=1)

    # Check on which side the density is smaller and assign it to be magnetopause
    # Check to see if b_gse_z_msp has same sign as vp_gse_z_msp
    for i in range(len(ind_msp)):
        if b_gse_vec_msp[i, 2] * vp_gse_vec_msp[i, 2] > 0:
            delta_v_th[i] = delta_v_th[i]
        else:
            delta_v_th[i] = - delta_v_th[i]

    delta_v_obs = vp_gse_vec_msh - vp_gse_vec_msp
    # delta_v_obs_mag = np.linalg.norm(delta_v_obs, axis=1)

    # Compute the angle between the observed and the theoretical velocity jumps
    theta_w = np.full(len(ind_msp), np.nan)
    for i in range(len(ind_msp)):
        theta_w[i] = np.arccos(np.dot(delta_v_th[i, :], delta_v_obs[i, :]) / (
                        np.linalg.norm(delta_v_th[i, :]) * np.linalg.norm(delta_v_obs[i, :])))

    # Convert angle to degrees
    theta_w_deg = theta_w * 180 / np.pi

    # Compute the ratio of the observed and the theoretical velocity jumps
    R_w = np.linalg.norm(delta_v_th, axis=1) / np.linalg.norm(delta_v_obs, axis=1)

    #print(f"Values of theta_w_deg: {theta_w_deg}")
    #print(f"Values of R_w: {R_w}")
    # Check if Walen relation is satisfied (0.4 < R_w < 3 and 0 < theta_w < 30 or 150 < theta_w <
    # 180) in green color

    bool_array = np.logical_and(np.logical_and(0.4 < R_w, R_w < 3),
                                np.logical_or(np.logical_and(0 < theta_w_deg, theta_w_deg < 30),
                                np.logical_and(150 < theta_w_deg, theta_w_deg < 180)))

    # Check the number of True in the bool_array
    num_true = np.sum(bool_array)
    if num_true >= n_points_walen/2:
        walen_relation_satisfied = True
    else:
        walen_relation_satisfied = False

    ind_walen_vals = np.flatnonzero(np.convolve(bool_array > 0,
                                    np.ones(n_points_walen, dtype=int),
                                    'valid') >= n_points_walen/10)
    
    if len(ind_walen_vals) > 0:
        walen_relation_satisfied_v2 = True
    else:
        walen_relation_satisfied_v2 = False

    if walen_relation_satisfied:
        print('\033[92m \n Walen relation satisfied \033[0m \n')
        print(f'\033[92m R_w: {np.nanmean(R_w):.3f} \033[0m')
        print(f'\033[92m theta_w_deg: {np.nanmean(theta_w_deg):.3f} \033[0m')
        if jet_detection:
            print(f'\033[92m Jet detection: {jet_detection} \033[0m')
        else:
            print(f'\033[91m Jet detection: {jet_detection} \033[0m')
    else:
        print('\033[91m \n Walen relation not satisfied \033[0m \n')
        print(f'\033[91m R_w: {np.nanmean(R_w):.3f} \033[0m')
        print(f'\033[91m theta_w_deg: {np.nanmean(theta_w_deg):.3f} \033[0m')
        if jet_detection:
            print(f'\033[92m Jet detection: {jet_detection} \033[0m')
        else:
            print(f'\033[91m Jet detection: {jet_detection} \033[0m')

    if walen_relation_satisfied_v2:
        print('\033[92m Walen relation version 2 is satisfied \033[0m \n')
    else:
        print('\033[91m Walen relation version 2 is not satisfied \033[0m \n')

    # Check if within 2 minutes of crossing time the values went above and below the threshold
    # If ind_vals is not empty, then append the crossing time to the csv file
    # if len(ind_vals) > 0:
    if walen_relation_satisfied or jet_detection or walen_relation_satisfied_v2:
        # for xxx in range(1):
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        #_ = spd.mms.mec
        #mms_sc_pos = ptt.get_data(mms_mec_varnames[0])[1:3][0][0] / r_e
        x = np.nanmean(mms_fgm_r_gsm[:, 0]) / r_e
        y = np.nanmean(mms_fgm_r_gsm[:, 1]) / r_e
        z = np.nanmean(mms_fgm_r_gsm[:, 2]) / r_e
        r_yz = np.sqrt(y**2 + z**2)  # Projection distance in yz plane.

        if x > -5 and x < 12 and r_yz < 12:
            # Check if the file exists, if nto then create it
            if not os.path.isfile(fname):
                with open(fname, 'w') as f:
                    f.write('Date,Probe,walen1,walen2,jet_detection,x_gsm,y_gsm,z_gsm,r_spc\n')
            # Append the crossing time to the csv file if it does not exist already
            df_temp = pd.read_csv(fname)
            old_crossing_times = np.array([xx[:19] for xx in df_temp['Date'].values])
            ttt = datetime.datetime.strftime(crossing_time, '%Y-%m-%d %H:%M:%S.%f')[:19]
            if ttt not in old_crossing_times:
                with open(fname, 'a') as f:
                    f.write(f'{crossing_time},{probe},{walen_relation_satisfied},{walen_relation_satisfied_v2},{jet_detection},{x:.3f},{y:.3f},{z:.3f},{r_yz:.3f}\n')
                f.close()

            # Set the fontstyle to Times New Roman
            font = {'family': 'serif', 'weight': 'normal', 'size': 10}
            plt.rc('font', **font)
            plt.rc('text', usetex=False)

            plt.close("all")
            # Add time to the figname
            if (walen_relation_satisfied or walen_relation_satisfied_v2) & jet_detection:
                folder_name = "../figures/jet_reversal_checks/jet_walen/b_n_v"
            elif (walen_relation_satisfied or walen_relation_satisfied_v2) & (not jet_detection):
                folder_name = "../figures/jet_reversal_checks/walen/b_n_v"
            elif (not walen_relation_satisfied) & (not walen_relation_satisfied_v2) & jet_detection:
                folder_name = "../figures/jet_reversal_checks/jet/b_n_v"
            bnv_figname = f"{folder_name}/b_n_v_{str(crossing_time.strftime('%Y%m%d_%H%M%S'))}"
            print(f"{bnv_figname}.png")
            # ptt.xlim(datetime.datetime.strftime(df_mms.index[0], "%Y-%m-%d %H:%M:%S.%f"),
            #          datetime.datetime.strftime(df_mms.index[-1], "%Y-%m-%d %H:%M:%S.%f"))
            ptt.tplot([f'mms{probe}_fgm_b_gsm_srvy_l2_bvec',
                       f'mms{probe}_dis_numberdensity_{data_rate}',
                       f'mms{probe}_dis_bulkv_gsm_{data_rate}'],
                      combine_axes=True, save_png=bnv_figname, display=False)
            plt.close("all")

            plt.figure(figsize=(6, 3))
            if (walen_relation_satisfied or walen_relation_satisfied_v2):
                plt.plot(df_mms.index, df_mms.vp_gsm_z_diff, 'g-', lw=1)
            else:
                plt.plot(df_mms.index, df_mms.vp_gsm_z_diff, 'b-', lw=1)
            if jet_detection:
                # Make a box around the jet
                plt.axvspan(vp_jet.index[ind_jet[0]], vp_jet.index[ind_jet[-1]], color='r',
                            alpha=0.2, label='Jet Location')
                plt.legend(loc=1)
            # plt.plot(vp_jet, 'r.', ms=2, lw=1, label="jet")
            # Draw a dashed line at +/- v_thres
            plt.axhline(y=v_thresh, color='k', linestyle='--', lw=1)
            plt.axhline(y=-v_thresh, color='k', linestyle='--', lw=1)
            plt.ylim(-200, 200)
            plt.xlim(df_mms.index[0], df_mms.index[-1])
            plt.xlabel("Time (UTC)")
            plt.ylabel("$v_p - <v_p>$ \n $(km/s, GSM, Z)$")
            temp3 = crossing_time.strftime('%Y-%m-%d %H:%M:%S')
            plt.title(f"MMS {probe} Jet Reversal Check at {temp3}")

            # Add text to the plot
            # plt.text(0.02, 0.98, f"$R_w$ = {R_w:.2f}\n $\\Theta_w$ = {theta_w_deg:.2f} \n $W_v$ =
            # {walen_relation_satisfied} \n $j_v$ = {jet_detection}",
            #          transform=plt.gca().transAxes, ha='left', va='top')
            if (walen_relation_satisfied or walen_relation_satisfied_v2) & jet_detection:
                folder_name = "../figures/jet_reversal_checks/jet_walen"
            elif (walen_relation_satisfied or walen_relation_satisfied_v2) & (not jet_detection):
                folder_name = "../figures/jet_reversal_checks/walen"
            elif (not walen_relation_satisfied) & (not walen_relation_satisfied_v2) & jet_detection:
                folder_name = "../figures/jet_reversal_checks/jet"
            ttt = str(crossing_time.strftime('%Y%m%d_%H%M%S'))
            fig_name = f"{folder_name}/mms{probe}_jet_reversal_check_{ttt}.png"
            plt.savefig(f"{fig_name}", dpi=150, bbox_inches='tight', pad_inches=0.1)
            print(f"{fig_name}")
            plt.close("all")

    return df_mms_fpi, df_mms_fgm, df_mms
