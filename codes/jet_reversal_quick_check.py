import datetime
import pytz
import importlib

import numpy as np
import pandas as pd

import jet_reversal_quick_check_function as jrcf
importlib.reload(jrcf)

# Read the list of dates from the csv file
df = pd.read_csv("../data/mms_magnetopause_crossings.csv")

# Convert the dates to datetime objects
#df["DateStart"] = pd.to_datetime(df["DateStart"])

# Set the index to the date column
df.set_index("DateStart", inplace=True)

# Set the timezone to UTC
#df.index = df.index.tz_localize(pytz.utc)

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

# Read the data from csv files
df_crossings = pd.read_csv("../data/mms_magnetopause_crossings.csv")
# Set the index to the date column
df_crossings.set_index("DateStart", inplace=True)

#for crossing_time in trange_list[28:29]:
for crossing_time in df_crossings.index[100:]:
    # Convert the crossing time to a datetime object
    # TODO: Something weird is happening with the timestamp. Check it later: crossing_time =
    # '2017-01-02 02:58:13.0+00:00'
    #crossing_time = '2017-01-02 02:58:13.0+00:00'
    crossing_time = datetime.datetime.strptime(crossing_time.split('+')[0], "%Y-%m-%d %H:%M:%S.%f")
    # Set the timezone to UTC
    crossing_time = crossing_time.replace(tzinfo=pytz.utc)
    # Try with 'brst' data rate, if that fails then try with 'fast'
    inputs = {'crossing_time': crossing_time,
              'dt': dt,
              'probe': probe,
              'jet_len': jet_len,
              'data_rate': data_rate,
              'level': level,
              'data_type': data_type,
              'time_clip': time_clip,
              'latest_version': latest_version,
              'figname': 'mms_jet_reversal_check',
              'fname': '../data/mms_jet_reversal_times.csv'
                }
    try:
        data_rate = 'brst'
        df_fpi, df_fgm, df_mms, df_mms_before, df_mms_after = jrcf.jet_reversal_check(**inputs)
    except:
        data_rate = 'fast'
        df_fpi, df_fgm, df_mms, df_mms_before, df_mms_after = jrcf.jet_reversal_check(**inputs)
