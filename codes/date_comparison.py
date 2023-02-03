import pandas as pd
import numpy as np

df_r = pd.read_csv("../data/rx_d/reconnection_line_data_mms3_20221109.csv")
df_s = pd.read_csv("../data/rx_d/reconnection_line_data_mms3_20221116.csv")

dates_r = df_r["date_from"].unique()
dates_s = df_s["date_from"].unique()

# Find x_gsm for each date dates_r
x_gsm_r = []
for date_r in dates_r:
    df_r_i = df_r[df_r["date_from"] == date_r]
    x_gsm_r.append(df_r_i["spc_pos_x"].unique()[0])

# Find x_gsm for each date dates_s
x_gsm_s = []
for date_s in dates_s:
    df_s_i = df_s[df_s["date_from"] == date_s]
    x_gsm_s.append(df_s_i["spc_pos_x"].unique()[0])

date_datetime_r = pd.to_datetime(dates_r)
date_datetime_s = pd.to_datetime(dates_s)

# Find the closest date in df_s to each date in df_r
ind = 0
for x_gsm_r_i,date_datetime in zip(x_gsm_r, date_datetime_r):
    # Find the minimum time difference between date_datetime and date_datetime_s
    arg_s = np.argmin(np.abs(date_datetime_s - date_datetime))

    # Find the difference between spacecraft location
    r_x_diff = np.abs(x_gsm_r_i - x_gsm_s[arg_s])
    min_time = np.nanmin(np.abs(date_datetime_s - date_datetime))

    # If min_time is less than 5 minutes, then increase the value of ind by 1
    if min_time < pd.Timedelta("10 minutes"):
        ind += 1
        print(f"r_x_diff = {r_x_diff}")

print(f"ind = {ind}")