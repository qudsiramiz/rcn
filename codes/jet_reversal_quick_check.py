import datetime
import pytz
import importlib

import numpy as np
import pandas as pd

import jet_reversal_quick_check_function as jrcf
importlib.reload(jrcf)

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

trange_list = [
    '2016-12-24 15:10:00.00',
    '2016-12-07 05:15:00.00',
    '2015-10-16 10:33:30.00',
    '2015-10-16 13:07:02.00',
    '2015-10-22 06:05:22.00',
    '2015-11-01 15:08:06.00',
    '2015-11-12 07:19:21.00',
    '2015-12-06 23:38:31.00',
    '2015-12-08 11:20:44.00',
    '2015-12-09 01:06:11.00',
    '2015-12-14 01:17:40.00',
    '2016-01-07 09:36:15.00',
    '2016-01-10 09:13:37.00',
    '2016-10-22 12:58:41.00',
    '2016-11-02 14:46:18.00',
    '2016-11-06 08:40:58.00',
    '2016-11-12 17:48:47.00',
    '2016-11-13 09:10:41.00',
    '2016-11-18 12:08:11.00',
    '2016-11-23 07:50:30.00',
    '2016-11-28 15:47:00.00',
    '2016-12-11 04:41:50.00',
    '2016-12-19 14:15:02.00',
    '2017-01-02 02:58:13.00',
    '2017-01-11 04:22:43.00',
    '2017-01-20 12:32:07.00',
    '2017-01-22 10:15:46.00',
    '2017-01-22 10:47:33.00',
    '2017-01-27 12:05:43.00',
    '2015-09-25 12:09:00.00',
    '2015-09-11 15:23:00.00',
    '2016-01-17 06:29:40.00',
    '2015-09-19 07:41:38.00',
    '2015-09-18 11:52:02.00'
]

trange_list = np.sort(trange_list)

# trange_list = [datetime.datetime(2016, 12, 28, 5, 38)]
# Read the data from csv files
df_crossings = pd.read_csv("../data/mms_jet_reversal_check_error_log_20220922.csv")
df_o = pd.read_csv("../data/mms_magnetopause_crossings.csv")
# Set the index to the date column
df_crossings.set_index("DateStart", inplace=True)
df_o.set_index("DateStart", inplace=True)

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
    indx_number = 4
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
                  'dt': dt,
                  'probe': probe,
                  'jet_len': jet_len,
                  'level': level,
                  'coord_type': 'lmn',
                  'data_type': data_type,
                  'time_clip': time_clip,
                  'latest_version': latest_version,
                  'figname': 'mms_jet_reversal_check_lmn',
                  #'fname': '../data/mms_jet_reversal_times_list_20220920.csv',
                  'fname': '../data/tst2.csv',
                  'error_file_log_name': "../data/mms_jet_reversal_check_error_log_20220922.csv",
                  "verbose": True
            }
        inputs["data_rate"] = 'brst'
        df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)
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
