import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import importlib
import datetime

import rx_model_funcs as rmf
importlib.reload(rmf)


trange_list = ["2015-09-02 16:45:59"]
dt = 0.5 # minutes

trange_date = datetime.datetime.strptime(trange_list[0], '%Y-%m-%d %H:%M:%S')
trange_date_min = trange_date - datetime.timedelta(minutes=dt)
trange_date_max = trange_date + datetime.timedelta(minutes=dt)
trange = [trange_date_min.strftime('%Y-%m-%d %H:%M:%S'),
          trange_date_max.strftime('%Y-%m-%d %H:%M:%S')]


sw_dict, df_fgm, df_fpi, df_mec, df_fgm_fpi = rmf.get_sw_params(trange=trange, mms_probe_num=3, verbose=True)

# Plot the df_mec data
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.plot(df_mec.index, df_mec.X, 'rd', label="$R_{x}$")
ax.plot(df_mec.index, df_mec.Y, 'gd', label="$R_{y}$")
ax.plot(df_mec.index, df_mec.Z, 'bd', label="$R_{z}$")
ax.set_xlabel("Time [UTC]")
ax.set_ylabel("Location [Re] (GSM)")
ax.set_xlim(df_mec.index[0], df_mec.index[-1])
ax.legend()
ax.set_title(f"MMS{3} MEC data for {trange[0]} to {trange[1]}")

plt.savefig(f"../figures/mec_data.png", dpi=300, bbox_inches="tight")
