import datetime
import pytz
import importlib

import numpy as np
import pandas as pd

import jet_reversal_quick_check_function_beta as jrcfb
importlib.reload(jrcfb)

# Read the list of dates from the csv file
# df = pd.read_csv("../data/mms_jet_reversal_check_error_log_20220922.csv")

# Convert the dates to datetime objects
# df["DateStart"] = pd.to_datetime(df["DateStart"])

# Set the index to the date column
# df.set_index("DateStart", inplace=True)

# Set the timezone to UTC
# df.index = df.index.tz_localize(pytz.utc)

r_e = 6378.137  # Earth radius in km

dt = 120
probe = 3
data_rate = 'brst'
level = 'l2'
data_type = 'dis-moms'
time_clip = True
latest_version = True
jet_len = 3

# trange_list = [datetime.datetime(2016, 12, 28, 5, 38)]
# Read the data from csv files
df_crossings = pd.read_csv("../data/mms_magnetopause_crossings.csv")
# Set the index to the date column
df_crossings.set_index("DateStart", inplace=True)

from contextlib import contextmanager,redirect_stderr,redirect_stdout
from os import devnull

@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

# Example usage
import sys

#with suppress_stdout_stderr():
for xxxx in range(1):
    indx_number = 179
    indx_max = indx_number + 1
    # for crossing_time in df_crossings.index[indx_number:indx_max]:
    for xx, crossing_time in enumerate(df_crossings.index[indx_number:indx_max], start=indx_number):
        # Convert the crossing time to a datetime object
        # TODO: Something weird is happening with the timestamp. Check it later: crossing_time =
        # '2017-01-02 02:58:13.0+00:00'
        # crossing_time = '2017-01-02 02:58:13.0+00:00'
        crossing_time = datetime.datetime.strptime(crossing_time.split('+')[0], "%Y-%m-%d %H:%M:%S.%f")
        print(crossing_time)
        # Set the timezone to UTC
        crossing_time = crossing_time.replace(tzinfo=pytz.utc)
        print(crossing_time)
        # Try with 'brst' data rate, if that fails then try with 'fast'
        inputs = {'crossing_time': crossing_time,
              'dt': 120,
              'probe': 3,
              'jet_len': 3,
              'level': 'l2',
              'coord_type': 'lmn',
              'data_type': 'dis-moms',
              'time_clip': True,
              'latest_version': True,
              'figname': 'mms_jet_reversal_check_lmn_mean',
              'fname': '../data/mms_jet_reversal_times_list_20220922_beta.csv',
              #'fname': '../data/test.csv',
              'error_file_log_name': "../data/mms_jet_reversal_check_error_log_20220922_beta.csv",
              "verbose": True
        }
        inputs["data_rate"] = 'brst'
        df_fpi, df_fgm, df_mms, rw, tw = jrcfb.jet_reversal_check(**inputs)
    #try:
    #    try:
    #        inputs["data_rate"] = 'brst'
    #        df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)
    #    except:
    #        inputs["data_rate"] = 'fast'
    #        df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)
    #except Exception as e:
    #    print(f"\033[91;31m\n{e} for date {crossing_time}\n\033[0m")
    #    pass

    # try:
    #     inputs["data_rate"] = 'brst'
    #     df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)
    # except:
    #     inputs["data_rate"] = 'fast'
    #     df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)

# np_std = (df_mms.np - df_mms.np.mean())/df_mms.np.std()
# vp_std = (df_mms.vp_gsm_z - df_mms.vp_gsm_z.mean())/df_mms.vp_gsm_z.std()
# plt.figure()
# plt.plot(np_std, 'r.', ms=1)
# plt.axvline(df_mms.index[ind_pos], color='k', linestyle='--', label='Positive')
# plt.axvline(df_mms.index[ind_neg], color='r', linestyle='--', label='Negative')
# plt.axvline(df_mms.index[ind_msp], color='c', linestyle='--', label='MSP')
# plt.axvline(df_mms.index[ind_msh], color='g', linestyle='-', lw=5, alpha=0.5, label='MSH')
# plt.legend()
# #plt.plot(vp_std, 'b.', ms=1)
# 
# plt.yscale('linear')
# plt.show()
