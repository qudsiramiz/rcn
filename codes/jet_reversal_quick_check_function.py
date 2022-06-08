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
                       fname='../data/mms_jet_reversal_times.csv'
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

    # Get the data from the FGM
    # mms_fgm_varnames = [f'mms{probe}_fgm_b_gsm_srvy_l2_bvec']
    _ = spd.mms.fgm(trange=trange, probe=probe, time_clip=time_clip, latest_version=True,
                    varnames=[f"mms{probe}_fgm_b_gsm_srvy_{level}",
                              f"mms{probe}_fgm_r_gsm_srvy_{level}"], get_fgm_ephemeris=True)

    # Get the time corresponding to the FGM data
    mms_fgm_time = ptt.get_data(f"mms{probe}_fgm_b_gsm_srvy_{level}")[0]
    mms_fgm_time = pd.to_datetime(mms_fgm_time, unit='s')
    mms_fgm_time = mms_fgm_time.tz_localize(pytz.utc)

    mms_fgm_b_gsm = ptt.get_data(f'mms{probe}_fgm_b_gsm_srvy_{level}')[1:4][0]
    mms_fgm_r_gsm = ptt.get_data(f'mms{probe}_fgm_r_gsm_srvy_{level}')[1:4][0]

    # Create a dataframe with the FPI data
    df_mms_fpi = pd.DataFrame(index=mms_fpi_time, data={'np': mms_fpi_numberdensity,
                                                        'vp_gsm_x': mms_fpi_bulkv_gsm[:, 0],
                                                        'vp_gsm_y': mms_fpi_bulkv_gsm[:, 1],
                                                        'vp_gsm_z': mms_fpi_bulkv_gsm[:, 2],
                                                        'tp_para': mms_fpi_temppara,
                                                        'tp_perp': mms_fpi_tempperp})

    # Make sure that the time indices are in increasing order
    df_mms_fpi = df_mms_fpi.sort_index()

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
        # 12 seconds. If not, then it is not a jet.
        # Based on the data rate, find out the number of data points in the time window which must
        # be checked for a jet
        # NOTE: 12 seconds comes from Trenchi et al. (2008) based on their study using Double Star
        # observation. Though they used it to confirm the Walen relation, I believe that similar
        # principal can be used for jet detection as well.
        n_points = int(jet_len / (df_mms_fpi.index[1] - df_mms_fpi.index[0]).total_seconds())

        # Find out the indices of all such points
        ind_pos_vals = np.flatnonzero(np.convolve(vp_jet > v_thresh,
                                      np.ones(n_points, dtype=int), 'valid') >= n_points)
        ind_neg_vals = np.flatnonzero(np.convolve(vp_jet < -v_thresh,
                                      np.ones(n_points, dtype=int), 'valid') >= n_points)

        if (len(ind_pos_vals) and len(ind_neg_vals)) > 0:
            jet_detection = True
            # Set the jet location to union of the positive and negative indices
            ind_jet = np.union1d(ind_pos_vals, ind_neg_vals)
        else:
            ind_jet = np.array([])

    # Create a dataframe with the FGM data
    df_mms_fgm = pd.DataFrame(index=mms_fgm_time, data={'b_gsm_x': mms_fgm_b_gsm[:, 0],
                                                        'b_gsm_y': mms_fgm_b_gsm[:, 1],
                                                        'b_gsm_z': mms_fgm_b_gsm[:, 2]})

    # Make sure all time indices are in increasing order
    df_mms_fgm = df_mms_fgm.sort_index()

    # Merge the two dataframes
    df_mms = pd.merge_asof(df_mms_fpi, df_mms_fgm, left_index=True, right_index=True)
    # Make the time index timezone aware
    df_mms.index = df_mms.index.tz_localize(pytz.utc)

    # Separate the data into two part, before and after the crossing time
    df_mms_before = df_mms[df_mms.index < crossing_time]
    df_mms_after = df_mms[df_mms.index >= crossing_time]

    # The following equations are based on the following paper:
    # Trenchi2008 (https://ui.adsabs.harvard.edu/abs/2008JGRA..113.7S10T)
    # Compute the average values for density, velocity and magnetic fields
    # before the crossing time
    np_before = df_mms_before['np'].mean() * 1e6  # in m^-3
    vp_gsm_x_before = df_mms_before['vp_gsm_x'].mean() * 1e3   # in m/s
    vp_gsm_y_before = df_mms_before['vp_gsm_y'].mean() * 1e3   # in m/s
    vp_gsm_z_before = df_mms_before['vp_gsm_z'].mean() * 1e3   # in m/s
    vp_gsm_vec_before = np.array([vp_gsm_x_before, vp_gsm_y_before, vp_gsm_z_before])
    tp_para_before = df_mms_before['tp_para'].mean() * 11600  # in K
    tp_perp_before = df_mms_before['tp_perp'].mean() * 11600  # in K
    b_gsm_x_before = df_mms_before['b_gsm_x'].mean() * 1e-9  # in T
    b_gsm_y_before = df_mms_before['b_gsm_y'].mean() * 1e-9  # in T
    b_gsm_z_before = df_mms_before['b_gsm_z'].mean() * 1e-9  # in T
    b_gsm_vec_before = np.array([b_gsm_x_before, b_gsm_y_before, b_gsm_z_before])

    # Compute the average values for density, velocity, magnetic fields and parallel and
    # perpendicular temperature
    # after the crossing time
    np_after = df_mms_after['np'].mean() * 1e6  # in m^-3
    vp_gsm_x_after = df_mms_after['vp_gsm_x'].mean() * 1e3  # in m/s
    vp_gsm_y_after = df_mms_after['vp_gsm_y'].mean() * 1e3  # in m/s
    vp_gsm_z_after = df_mms_after['vp_gsm_z'].mean() * 1e3  # in m/s
    vp_gsm_vec_after = np.array([vp_gsm_x_after, vp_gsm_y_after, vp_gsm_z_after])
    tp_para_after = df_mms_after['tp_para'].mean() * 11600  # in K
    tp_perp_after = df_mms_after['tp_perp'].mean() * 11600  # in K
    b_gsm_x_after = df_mms_after['b_gsm_x'].mean() * 1e-9  # in T
    b_gsm_y_after = df_mms_after['b_gsm_y'].mean() * 1e-9  # in T
    b_gsm_z_after = df_mms_after['b_gsm_z'].mean() * 1e-9  # in T
    b_gsm_vec_after = np.array([b_gsm_x_after, b_gsm_y_after, b_gsm_z_after])

    # Define the mass of proton in kg
    m_p = 1.6726219e-27

    # Define the absolute permeability of free space in m^2 kg^-1 s^-1
    mu_0 = 4 * np.pi * 1e-7

    # Define the Boltzmann constant in J K^-1
    k_B = 1.38064852e-23

    alpha_before = (mu_0 * np_before * k_B) * (tp_para_before - tp_perp_before) / (
        np.linalg.norm(b_gsm_vec_before)**2)

    alpha_after = (mu_0 * np_after * k_B) * (tp_para_after - tp_perp_after) / (
        np.linalg.norm(b_gsm_vec_after)**2)

    v_th_before = b_gsm_vec_before * (1 - alpha_before) / (
        mu_0 * np_before * m_p * (1 - alpha_before)
    )**0.5

    v_th_after = b_gsm_vec_after * (1 - alpha_after) / (
        mu_0 * np_after * m_p * (1 - alpha_after)
    )**0.5

    delta_v_th = v_th_after - v_th_before
    delta_v_th_mag = np.linalg.norm(delta_v_th)

    delta_v_obs = vp_gsm_vec_after - vp_gsm_vec_before
    delta_v_obs_mag = np.linalg.norm(delta_v_obs)

    # Compute the angle between the observed and the theoretical velocity jumps
    theta_w = np.arccos(np.dot(delta_v_th, delta_v_obs) / (
                        np.linalg.norm(delta_v_th) * np.linalg.norm(delta_v_obs)))

    # Convert angle to degrees
    theta_w_deg = theta_w * 180 / np.pi

    # Compute the ratio of the observed and the theoretical velocity jumps
    R_w = np.linalg.norm(delta_v_th) / np.linalg.norm(delta_v_obs)

    # Check if Walen relation is satisfied (0.4 < R_w < 3 and 0 < theta_w < 30 or 150 < theta_w <
    # 180) in green color

    if 0.4 < R_w < 3 and (0 < theta_w_deg < 30 or 150 < theta_w_deg < 180):
        print('\033[92m Walen relation satisfied \033[0m \n')
        print(f'\033[92m R_w: {R_w} \033[0m')
        print(f'\033[92m theta_w: {theta_w_deg} \033[0m')
        print(f'\033[92m delta_v_th: {delta_v_th_mag / 1e3} km/s \033[0m')
        print(f'\033[92m delta_v_obs: {delta_v_obs_mag / 1e3} km/s \033[0m')
        print(f'\033[92m Jet detection: {jet_detection} \033[0m')
        print('\n')
        walen_relation_satisfied = True
    else:
        print('\033[91m Walen relation not satisfied \033[0m \n')
        print(f'\033[91m R_w: {R_w} \033[0m')
        print(f'\033[91m theta_w: {theta_w_deg} \033[0m')
        print(f'\033[91m delta_v_th: {delta_v_th_mag / 1e3} km/s \033[0m')
        print(f'\033[91m delta_v_obs: {delta_v_obs_mag / 1e3} km/s \033[0m')
        print(f'\033[91m Jet detection: {jet_detection} \033[0m')
        print('\n')
        walen_relation_satisfied = False

    # Check if within 2 minutes of crossing time the values went above and below the threshold

    # If ind_vals is not empty, then append the crossing time to the csv file
    # if len(ind_vals) > 0:
    if walen_relation_satisfied or jet_detection:
        # for xxx in range(1):
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        # mms_sc_pos = ptt.get_data(mms_mec_varnames[0])[1:3][0][0] / r_e
        x = np.nanmean(mms_fgm_r_gsm[0]) / r_e
        y = np.nanmean(mms_fgm_r_gsm[1]) / r_e
        z = np.nanmean(mms_fgm_r_gsm[2]) / r_e
        r_yz = np.sqrt(y**2 + z**2)  # Projection distance in yz plane.

        if x > -5 and x < 12 and r_yz < 12:
            # Check if the file exists, if nto then create it
            if not os.path.isfile(fname):
                with open(fname, 'w') as f:
                    f.write('Date, Probe,walen,jet_detection,R_w,theta_w,x_gsm,y_gsm,z_gsm,\
                            r_spc\n')
            # Append the crossing time to the csv file if it does not exist already
            df_temp = pd.read_csv(fname)
            old_crossing_times = np.array([xx[:19] for xx in df_temp['Date'].values])
            ttt = datetime.datetime.strftime(crossing_time, '%Y-%m-%d %H:%M:%S.%f')[:19]
            if ttt not in old_crossing_times:
                with open(fname, 'a') as f:
                    f.write(f'{crossing_time},{probe},{walen_relation_satisfied},{jet_detection}\
                              ,{R_w:.3f},{theta_w_deg:.3f},{x:.3f},{y:.3f},{z:.3f},{r_yz:.3f}\
                              \n')
                f.close()

            # Set the fontstyle to Times New Roman
            font = {'family': 'serif', 'weight': 'normal', 'size': 10}
            plt.rc('font', **font)
            plt.rc('text', usetex=False)

            plt.close("all")
            # Add time to the figname
            figname = f"../figures/jet_reversal_checks/{figname}_{str(crossing_time.strftime('%Y%m%d_%H%M%S'))}"
            print(figname)
            ptt.tplot([f'mms{probe}_fgm_b_gsm_srvy_l2_bvec',
                       f'mms{probe}_dis_numberdensity_{data_rate}',
                       f'mms{probe}_dis_bulkv_gsm_{data_rate}'],
                      combine_axes=True, save_png=figname, display=False)

            plt.close("all")

            plt.figure(figsize=(6, 3))
            plt.plot(df_mms.index, df_mms.vp_gsm_z_diff, 'b-', lw=1)
            if jet_detection:
                plt.plot(vp_jet[ind_jet], 'rd', ms=2, lw=1, label="jet location")
            # plt.plot(vp_jet, 'r.', ms=2, lw=1, label="jet")
            # Draw a dashed line at +/- v_thres
            plt.axhline(y=v_thresh, color='k', linestyle='--', lw=1)
            plt.axhline(y=-v_thresh, color='k', linestyle='--', lw=1)
            plt.ylim(-200, 200)
            plt.xlabel("Time (UTC)")
            plt.ylabel("$v_p - <v_p>$ \n $(km/s, GSM, Z)$")
            plt.title(f"MMS {probe} Jet Reversal Check at {crossing_time.strftime('%Y-%m-%d %H:%M:%S')}")
            plt.legend(loc=1)

            # Add text to the plot
            plt.text(0.02, 0.98, f"$R_w$ = {R_w:.2f}\n $\Theta_w$ = {theta_w_deg:.2f} \n $W_v$ = \
                {walen_relation_satisfied}, \n $j_v$ = {jet_detection}",
                     transform=plt.gca().transAxes, ha='left', va='top')
            fig_name = f"../figures/jet_reversal_checks/jet/mms{probe}_jet_reversal_check_{str(crossing_time.strftime('%Y%m%d_%H%M%S'))}.png"
            plt.savefig(f"{fig_name}", dpi=150, bbox_inches='tight', pad_inches=0.1)
            print(f"{fig_name}")

    return df_mms_fpi, df_mms_fgm, df_mms, df_mms_before, df_mms_after
