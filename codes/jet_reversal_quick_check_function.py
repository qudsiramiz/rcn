import datetime
import os
import pytz
import tabulate
import importlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyspedas as spd
import pytplot as ptt

import matplotlib.gridspec as gridspec

from matplotlib.pyplot import MaxNLocator


def jet_reversal_check(crossing_time=None, dt=90, probe=3, data_rate='fast', level='l2',
                       data_type='dis-moms', time_clip=True, latest_version=True,
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
    mms_fgm_varnames = [f'mms{probe}_fgm_b_gsm_srvy_l2_bvec']
    _ = spd.mms.fgm(trange=trange, probe=probe, time_clip=time_clip, latest_version=True,
                               varnames=[f"mms{probe}_fgm_b_gsm_srvy_{level}",
                                        f"mms{probe}_fgm_r_gsm_srvy_{level}"],
                                        get_fgm_ephemeris=True)

    # Get the time corresponding to the FGM data
    mms_fgm_time = ptt.get_data(f"mms{probe}_fgm_b_gsm_srvy_{level}")[0]
    mms_fgm_time = pd.to_datetime(mms_fgm_time, unit='s')
    mms_fgm_time = mms_fgm_time.tz_localize(pytz.utc)

    mms_fgm_b_gsm = ptt.get_data(f'mms{probe}_fgm_b_gsm_srvy_{level}')[1:4][0]
    mms_fgm_r_gsm = ptt.get_data(f'mms{probe}_fgm_r_gsm_srvy_{level}')[1:4][0]
    #print(f'\033[92m The length of the b_gsm is {len(mms_fgm_b_gsm)} \033[0m')

    # Create a dataframe with the FPI data
    df_mms_fpi = pd.DataFrame(index=mms_fpi_time, data={'np': mms_fpi_numberdensity,
                                                        'vp_gsm_x': mms_fpi_bulkv_gsm[:,0],
                                                        'vp_gsm_y': mms_fpi_bulkv_gsm[:,1],
                                                        'vp_gsm_z': mms_fpi_bulkv_gsm[:,2],
                                                        'tp_para': mms_fpi_temppara,
                                                        'tp_perp': mms_fpi_tempperp})

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

    # Create a dataframe with the FGM data
    df_mms_fgm = pd.DataFrame(index=mms_fgm_time, data={'b_gsm_x': mms_fgm_b_gsm[:,0],
                                                        'b_gsm_y': mms_fgm_b_gsm[:,1],
                                                        'b_gsm_z': mms_fgm_b_gsm[:,2]})

    # Compute the difference in velocity between the two points separated by 2 minutes
    periods = int(dt / (df_mms_fpi.index[1] - df_mms_fpi.index[0]).total_seconds())
    df_mms_fpi['vp_diff'] = abs(df_mms_fpi['vp_gsm_z_diff'] -
                                df_mms_fpi['vp_gsm_z_diff'].shift(periods=periods))

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
    vp_gsm_x_before = df_mms_before['vp_gsm_x'].mean() * 1e3  # in m/s
    vp_gsm_y_before = df_mms_before['vp_gsm_y'].mean() * 1e3  # in m/s
    vp_gsm_z_before = df_mms_before['vp_gsm_z'].mean() * 1e3  # in m/s
    vp_gsm_vec_before = np.array([vp_gsm_x_before, vp_gsm_y_before, vp_gsm_z_before])
    tp_para_before = df_mms_before['tp_para'].mean() * 11600 # in K
    tp_perp_before = df_mms_before['tp_perp'].mean() * 11600 # in K
    b_gsm_x_before = df_mms_before['b_gsm_x'].mean() * 1e-9 # in T
    b_gsm_y_before = df_mms_before['b_gsm_y'].mean() * 1e-9 # in T
    b_gsm_z_before = df_mms_before['b_gsm_z'].mean() * 1e-9 # in T
    b_gsm_vec_before = np.array([b_gsm_x_before, b_gsm_y_before, b_gsm_z_before])

    # Compute the average values for density, velocity, magnetic fields and parallel and
    # perpendicular temperature
    # after the crossing time
    np_after = df_mms_after['np'].mean() * 1e6  # in m^-3
    vp_gsm_x_after = df_mms_after['vp_gsm_x'].mean() * 1e3  # in m/s
    vp_gsm_y_after = df_mms_after['vp_gsm_y'].mean() * 1e3  # in m/s
    vp_gsm_z_after = df_mms_after['vp_gsm_z'].mean() * 1e3  # in m/s
    vp_gsm_vec_after = np.array([vp_gsm_x_after, vp_gsm_y_after, vp_gsm_z_after])
    tp_para_after = df_mms_after['tp_para'].mean() * 11600 # in K
    tp_perp_after = df_mms_after['tp_perp'].mean() * 11600 # in K
    b_gsm_x_after = df_mms_after['b_gsm_x'].mean() * 1e-9 # in T
    b_gsm_y_after = df_mms_after['b_gsm_y'].mean() * 1e-9 # in T
    b_gsm_z_after = df_mms_after['b_gsm_z'].mean() * 1e-9 # in T
    b_gsm_vec_after = np.array([b_gsm_x_after, b_gsm_y_after, b_gsm_z_after])

    # # Print all the after and before values
    # print('Before crossing time:')
    # print('np:', np_before)
    # print('vp_gsm_x:', vp_gsm_x_before)
    # print('vp_gsm_y:', vp_gsm_y_before)
    # print('vp_gsm_z:', vp_gsm_z_before)
    # print('tp_para:', tp_para_before)
    # print('tp_perp:', tp_perp_before)
    # print('b_gsm_x:', b_gsm_x_before)
    # print('b_gsm_y:', b_gsm_y_before)
    # print('b_gsm_z:', b_gsm_z_before)
    # print('After crossing time:')
    # print('np:', np_after)
    # print('vp_gsm_x:', vp_gsm_x_after)
    # print('vp_gsm_y:', vp_gsm_y_after)
    # print('vp_gsm_z:', vp_gsm_z_after)
    # print('tp_para:', tp_para_after)
    # print('tp_perp:', tp_perp_after)
    # print('b_gsm_x:', b_gsm_x_after)
    # print('b_gsm_y:', b_gsm_y_after)
    # print('b_gsm_z:', b_gsm_z_after)

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

    v_th_before = b_gsm_vec_before * ( 1 - alpha_before) / (
                  mu_0 * np_before * m_p * (1 - alpha_before)
    )**0.5

    v_th_after = b_gsm_vec_after * ( 1 - alpha_after) / (
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

    if 0.4 < R_w < 3 and (0 < theta_w < 30 or 150 < theta_w < 180):
        print('\033[92m Walen relation satisfied \033[0m \n\n\n')
        print('R_w = ', R_w)
        print('theta_w = ', theta_w_deg)
        print('delta_v_th = ', delta_v_th_mag)
        print('delta_v_obs = ', delta_v_obs_mag)
        print('\n')
        walen_relation_satisfied = True
    else:
        print('\033[91m Walen relation not satisfied \033[0m \n\n\n')
        print('R_w = ', R_w)
        print('theta_w = ', theta_w_deg)
        print('delta_v_th = ', delta_v_th_mag)
        print('delta_v_obs = ', delta_v_obs_mag)
        print('\n')
        walen_relation_satisfied = False

    # Check if within 2 minutes of crossing time the values went above and below the threshold
    
    # Compute the angle between the observed and the
    # Check if at least three consecutive points have a difference greater than 70
    # km/s
    v_thresh = 70
    N = 3
    ind_vals = np.flatnonzero(np.convolve(df_mms_fpi['vp_diff'] > v_thresh,
                              np.ones(N, dtype=int), 'valid')>=N)

    # If ind_vals is not empty, then append the crossing time to the csv file
    if len(ind_vals) > 0:

        # Check if at the time of crossing, the probe is actually on the dayside of the
        # magnetopause.
        # NOTE: This is a very rough check, and should be improved. Also, the way I describe the
        # "dayside" of the magnetopause is not based on a scientific understanding of the
        # magnetopause.

        # Get the position of mms spacecraft in gsm coordinates
        #mms_mec_varnames = [f'mms{probe}_mec_r_gsm']
        #mms_mec_vars = spd.mms.mec(trange=trange, varnames=mms_mec_varnames, probe=probe,
        #                       data_rate='srvy', level='l2', time_clip=time_clip,
        #                       latest_version=latest_version)
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
                    f.write('Date, Probe, walen, R_w, theta_w, x_gsm, y_gsm, z_gsm, r_spc\n')
            # Append the crossing time to the csv file if it does not exist already
            df_temp = pd.read_csv(fname)
            old_crossing_times = np.array([xx[:19] for xx in df_temp['Date'].values])
            ttt = datetime.datetime.strftime(crossing_time, '%Y-%m-%d %H:%M:%S.%f')[:19]
            if ttt not in old_crossing_times:
                with open(fname, 'a') as f:
                    f.write(f'{crossing_time}, {probe}, {walen_relation_satisfied}, {R_w},\
                              {theta_w}, {x}, {y}, {z}, {r_yz}\n')
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
        plt.figure(figsize=(10,6))
        plt.plot(df_mms.index, df_mms.vp_gsm_z_diff, 'b-', lw=1, label="diff")
        plt.ylim(-100, 100)
        plt.xlabel("Time (UTC)")
        plt.ylabel("Vp (km/s)")
        plt.title(f"MMS {probe} Jet Reversal Check at {crossing_time.strftime('%Y-%m-%d %H:%M:%S')}")
        plt.legend()
        fig_name= f"../figures/jet_reversal_checks/jet/mms{probe}_jet_reversal_check_{str(crossing_time.strftime('%Y%m%d_%H%M%S'))}.png"
        plt.savefig(f"{fig_name}", dpi=150, bbox_inches='tight', pad_inches=0.1)
        print(f"{fig_name}")

    return df_mms_fpi, df_mms_fgm, df_mms, df_mms_before, df_mms_after
