import datetime
import os
import time
import pytz

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyspedas as spd
import pytplot as ptt

start_time = time.time()


def jet_reversal_check(crossing_time=None, dt=90, probe=3, data_rate='fast', level='l2',
                       data_type='dis-moms', time_clip=True, latest_version=True,
                       fname='mms_jet_reversal_times.csv'):
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
    None.
    """
    crossing_time_min = crossing_time - datetime.timedelta(seconds=dt)
    crossing_time_max = crossing_time + datetime.timedelta(seconds=dt)
    trange = [crossing_time_min, crossing_time_max]

    mms_fpi_varnames = [f'mms{probe}_dis_numberdensity_{data_rate}',
                        f'mms{probe}_dis_bulkv_gse_{data_rate}']
    mms_fpi_vars = spd.mms.fpi(trange=trange, probe=probe, data_rate=data_rate, level=level,
                               datatype=data_type, varnames=mms_fpi_varnames, time_clip=time_clip,
                               latest_version=latest_version)
    # Convert from GSE to GSM data

    mms_fpi_time = ptt.get_data(mms_fpi_varnames[0])[0]
    mms_fpi_numberdensity = ptt.get_data(mms_fpi_varnames[0])[1]
    mms_fpi_bulkv_gse = ptt.get_data(mms_fpi_varnames[1])[1:4][0]
    # Covert gse to gsm
    _ = spd.cotrans(name_in=f'mms{probe}_dis_bulkv_gse_{data_rate}', name_out='v_gsm',
                            coord_in='gse', coord_out='gsm')
    
    mms_fpi_bulkv_gsm = ptt.get_data('v_gsm')[1:4][0]

    # Convert the time to a datetime object
    mms_fpi_time = pd.to_datetime(mms_fpi_time, unit='s')

    # Create a dataframe with the data
    df_mms_fpi = pd.DataFrame(index=mms_fpi_time, data={'np': mms_fpi_numberdensity,
                                                        'vp_gsm_x': mms_fpi_bulkv_gsm[:,0],
                                                        'vp_gsm_y': mms_fpi_bulkv_gsm[:,1],
                                                        'vp_gsm_z': mms_fpi_bulkv_gsm[:,2]})

    # 
    # Compute the difference in velocity between the two points separated by 2 minutes
    periods = int(dt / (df_mms_fpi.index[1] - df_mms_fpi.index[0]).total_seconds())
    df_mms_fpi['vp_diff'] = abs(df_mms_fpi['vp_gsm_z'] -
                                df_mms_fpi['vp_gsm_z'].shift(periods=periods))

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
        mms_mec_varnames = [f'mms{probe}_mec_r_gsm']
        mms_mec_vars = spd.mms.mec(trange=trange, varnames=mms_mec_varnames, probe=probe,
                               data_rate='srvy', level='l2', time_clip=time_clip,
                               latest_version=latest_version)
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        mms_sc_pos = ptt.get_data(mms_mec_varnames[0])[1:3][0][0] / r_e
        x = mms_sc_pos[0]
        y = mms_sc_pos[1]
        z = mms_sc_pos[2]
        r_yz = np.sqrt(y**2 + z**2)  # Projection distance in yz plane.

        if x > -5 and x < 12 and r_yz < 12:
            # Check if the file exists, if nto then create it
            if not os.path.isfile(fname):
                with open(fname, 'w') as f:
                    f.write('Date, Probe, x_gsm, y_gsm, z_gsm, r_spc\n')
            # Append the crossing time to the csv file
            with open(fname, 'a') as f:
                f.write(f'{crossing_time},{probe}, {x},{y},{z},{r_yz}\n')
                f.close()

    return ind_vals

# Read the list of dates from the csv file
df = pd.read_csv("../data/mms_magnetopause_crossings.csv")

# Convert the dates to datetime objects
#df["DateStart"] = pd.to_datetime(df["DateStart"])

# Set the index to the date column
df.set_index("DateStart", inplace=True)

# Set the timezone to UTC
#df.index = df.index.tz_localize(pytz.utc)

jet_reversal_inputs = {
    "dt" : 120,
    "probe" : 3,
    "data_rate" : 'fast',
    "level" : 'l2',
    "data_type" : 'dis-moms',
    "time_clip" : True,
    "latest_version" : True,
    "fname" : '../data/mms_jet_reversal_times.csv'
}

for i, crossing_time in enumerate(df.index[1000:]):
    # Convert the crossing time to a datetime object
    crossing_time = datetime.datetime.strptime(crossing_time, "%Y-%m-%d %H:%M:%S.%f")
    # Set the timezone to UTC
    crossing_time = crossing_time.replace(tzinfo=pytz.utc)
    ind_vals = jet_reversal_check(crossing_time=crossing_time, **jet_reversal_inputs)
    # Print i in the terminal in green color
    if len(ind_vals) > 0:
        print(f'\033[92m{i}\033[0m\n\n')


#plt.figure()
## Convert t from unix seconds to datetime
#t = pd.to_datetime(t, unit='s')
#plt.plot(df.np)
#plt.show()

#plt.figure()
#plt.plot(df.vp_gsm_x, 'b-', label='Vx GSM')
#plt.plot(df.vp_gsm_y, 'g-', label='Vy GSM')
#plt.plot(df.vp_gsm_z, 'r-', label='Vz GSM')
#plt.show()

