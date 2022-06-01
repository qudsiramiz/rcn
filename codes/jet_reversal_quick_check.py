import datetime
import os
import time
import pytz

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyspedas as spd
import pytplot as ptt

import matplotlib.gridspec as gridspec

from matplotlib.pyplot import MaxNLocator

start_time = time.time()


def jet_reversal_check(crossing_time=None, dt=90, probe=3, data_rate='fast', level='l2',
                       data_type='dis-moms', time_clip=True, latest_version=True,
                       figname='mms_jet_reversal_check'):
    """
    """
    crossing_time_min = crossing_time - datetime.timedelta(seconds=dt)
    crossing_time_max = crossing_time + datetime.timedelta(seconds=dt)
    trange = [crossing_time_min, crossing_time_max]
    print(trange)
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

    _ = spd.mms.fgm(trange=trange, probe=probe, time_clip=time_clip, latest_version=True)

    # Create a dataframe with the data
    df_mms_fpi = pd.DataFrame(index=mms_fpi_time, data={'np': mms_fpi_numberdensity,
                                                        'vp_gsm_x': mms_fpi_bulkv_gsm[:,0],
                                                        'vp_gsm_y': mms_fpi_bulkv_gsm[:,1],
                                                        'vp_gsm_z': mms_fpi_bulkv_gsm[:,2]})

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
              'v_gsm'],
              combine_axes=True, save_png=figname, display=False)
  
    return df_mms_fpi

# Read the list of dates from the csv file
df = pd.read_csv("../data/mms_jet_reversal_times.csv")

# Convert the dates to datetime objects
#df["Date"] = pd.to_datetime(df["Date"])

# Set the index to the date column
df.set_index("Date", inplace=True)

# Set the timezone to UTC
#df.index = df.index.tz_localize(pytz.utc)

r_e = 6378.137  # Earth radius in km

dt = 120
probe = 3
data_rate = 'fast'
level = 'l2'
data_type = 'dis-moms'
time_clip = True
latest_version = True

for crossing_time in df.index[0:]:
    # Convert the crossing time to a datetime object
    crossing_time = datetime.datetime.strptime(crossing_time.split('+')[0], "%Y-%m-%d %H:%M:%S.%f")
    # Set the timezone to UTC
    crossing_time = crossing_time.replace(tzinfo=pytz.utc)
    df = jet_reversal_check(crossing_time=crossing_time, dt=dt, probe=probe,
                            data_rate=data_rate, level=level, data_type=data_type,
                            time_clip=time_clip, latest_version=latest_version,
                            figname='mms_jet_reversal_check')
