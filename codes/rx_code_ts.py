# This the python version of IDL code named "RX_model_batch.pro"
import datetime
import importlib
import os
import time
from contextlib import contextmanager, redirect_stderr, redirect_stdout

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import rx_model_funcs as rmf

importlib.reload(rmf)

# Set the fontstyle to Times New Roman
font = {"family": "serif", "weight": "normal", "size": 10}
plt.rc("font", **font)
plt.rc("text", usetex=True)


@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


start = time.time()

today_date = datetime.datetime.today().strftime("%Y-%m-%d")


# Define a center time (centered on the jet location)
# center_time_list = ["2015-09-07 14:12:33", "2015-09-19 09:37:00", "2016-12-26 11:10:00",
#                     "2016-11-15 13:28:00", "2017-02-04 01:04:00", "2017-01-24 03:56:00",
#                     "2017-01-21 00:09:00", "2017-01-09 02:58:00", "2017-01-05 01:21:00"]

center_time_list = ["2015-09-02 16:47:33"]
ind_min = 0
ind_max = 1

for center_time_str in center_time_list[ind_min:ind_max]:
    center_time = pd.to_datetime(center_time_str).tz_localize("UTC")

    # Get a start/end time 1 hour before/after the center time and set the time zone to UTC
    start_time = (center_time - pd.Timedelta("2 minute")).tz_convert("UTC")
    end_time = (center_time + pd.Timedelta("2 minute")).tz_convert("UTC")

    # Define an array of time ranges to loop over, between the start and end times with 1 minute
    # intervals
    ind_t_min = 200
    ind_t_max = 241
    trange_list = pd.date_range(start_time, end_time, freq="1s").tz_convert("UTC")
    for ind_range, trange in enumerate(trange_list[ind_t_min:ind_t_max:1], start=ind_t_min):
        trange_min = trange - pd.Timedelta("0.5 minute")
        trange_max = trange + pd.Timedelta("0.5 minute")
        trange = [trange_min.strftime('%Y-%m-%d %H:%M:%S'), trange_max.strftime('%Y-%m-%d %H:%M:%S')]
        # with suppress_stdout_stderr():
        #     sw_dict, df_mms, df_mec, df_fgm, df_fpi = rmf.get_sw_params(trange=trange,
        #                                                                 mms_probe_num=3,
        #                                                                 verbose=False)

        # print(f"\033[1;32m{sw_dict['mms_b_gsm'][:]} \033[0m")
        # print(f"\033[1;32m{sw_dict['mms_v_gsm'][:]} \033[0m")
        # print(f"\033[1;32m{df_fgm.index.min()} \033[0m")
        # print(f"\033[1;32m{df_fgm.index.max()} \033[0m")
        # print(f"\033[1;32m{trange} \033[0m")

        # Print the ind_range  which is being processed in green color every 10th time
        if ind_range % 5 == 0:
            print(f"\n Processing for index \033[92m{ind_range}\033[0m of "
                  f"\033[91m {len(trange_list)}\033[0m\n")

        # Convert trange to string to format '%Y-%m-%d %H:%M:%S'
        # trange = trange.strftime('%Y-%m-%d %H:%M:%S')
        # trange = [trange.split("+")[0].split(".")[0]]
        # trange = ["2015-9-9 14:11:14"]
        # print(trange)
        with suppress_stdout_stderr():
        # for foo in range(1):
        #     for bar in range(1):
            try:
                mms_probe_num = str(3)
                min_max_val = 20
                dr = 0.25
                y_min = - min_max_val
                y_max = min_max_val
                z_min = - min_max_val
                z_max = min_max_val
                model_type = "t96"

                model_inputs = {
                    "trange": trange,
                    "dt": 2,
                    "probe": None,
                    "omni_level": "hro",
                    "mms_probe_num": 3,
                    "model_type": model_type,
                    "m_p": 1,
                    "dr": dr,
                    "min_max_val": min_max_val,
                    "y_min": y_min,
                    "y_max": y_max,
                    "z_min": z_min,
                    "z_max": z_max,
                    "save_data": False,
                    "nprocesses": None,
                }
                (bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh, sw_params, x_shu, y_shu,
                 z_shu, b_msx, b_msy, b_msz) = rmf.rx_model(**model_inputs)

                # mask = shear > 175
                # shear[mask] = 0
                # Normalize each quantity to the range [0, 1]
                shear_norm = (shear - np.nanmin(shear)) / (np.nanmax(shear) - np.nanmin(shear))
                rx_en_norm = (rx_en - np.nanmin(rx_en)) / (np.nanmax(rx_en) - np.nanmin(rx_en))
                va_cs_norm = (va_cs - np.nanmin(va_cs)) / (np.nanmax(va_cs) - np.nanmin(va_cs))
                bisec_msp_norm = (bisec_msp - np.nanmin(bisec_msp)) / (np.nanmax(bisec_msp) -
                                                                       np.nanmin(bisec_msp))

                # shear_norm = (shear - np.nanmin(shear)) / (np.std(shear))
                # rx_en_norm = (rx_en - np.nanmin(rx_en)) / (np.std(rx_en))
                # va_cs_norm = (va_cs - np.nanmin(va_cs)) / (np.std(va_cs))
                # bisec_msp_norm = (bisec_msp - np.nanmin(bisec_msp)) / (np.std(bisec_msp))
                # bisec_msp_norm = (bisec_msh - np.nanmin(bisec_msh)) / (np.nanmax(bisec_msh) -
                # np.nanmin(bisec_msh))
                figure_inputs = {
                    "image": [shear_norm, rx_en_norm, va_cs_norm, bisec_msp_norm],
                    "convolution_order": [1, 1, 1, 1],
                    "t_range": trange,
                    "dt": 2,
                    "b_imf": np.round(sw_params["b_imf"], 2),
                    "b_msh": np.round(sw_params["mms_b_gsm"], 2),
                    "v_msh": np.round(sw_params["mms_v_gsm"], 2),
                    "xrange": [y_min, y_max],
                    "yrange": [z_min, z_max],
                    "mms_probe_num": 3,
                    "mms_sc_pos": np.round(sw_params["mms_sc_pos"], 2),
                    "dr": dr,
                    "dipole_tilt_angle": sw_params["ps"],
                    "p_dyn": np.round(sw_params["p_dyn"], 2),
                    "imf_clock_angle": sw_params["imf_clock_angle"],
                    "sigma": [2, 2, 2, 2],
                    "mode": "nearest",
                    "alpha": 1,
                    "vmin": [0, 0, 0, 0],
                    "vmax": [1, 1, 1, 1],
                    "cmap_list": ["viridis", "cividis", "plasma", "magma"],
                    "draw_patch": [True, True, True, True],
                    "draw_ridge": [True, True, True, True],
                    "save_fig": True,
                    "fig_name": "crossing_all_ridge_plots",
                    # "fig_format": "png",
                    "c_label": ["Shear", "Reconnection Energy", "Exhaust Velocity",
                                "Bisection Field"],
                    # "c_unit": [r"${}^\circ$", "nPa", "km/s", "nT"],
                    "c_unit": ["", "nPa", "km/s", "nT"],
                    "wspace": 0.0,
                    "hspace": 0.17,
                    "fig_size": (8.775, 10),
                    "box_style": dict(boxstyle="round", color="k", alpha=0.8),
                    # "box_style": dict(boxstyle="round", color=None, alpha=0.8),
                    "title_y_pos": 1.09,
                    "interpolation": "None",
                    "tsy_model": model_type,
                    "dark_mode": True,
                    "rc_file_name": f"rx_line_data_time_series_mms{3}_"
                                    f"{center_time.strftime('%Y%m%d_%H%M%S')}.csv",
                    "rc_folder": "../data/rx_d/time_series/",
                    "save_rc_file": True,
                    "fig_version": f"time_series_{center_time.strftime('%Y%m%d_%H%M%S')}",
                }

                y_vals, x_intr_vals_list, y_intr_vals_list = rmf.ridge_finder_multiple(
                                                            **figure_inputs, fig_format="png")
                print(f"\033[92m \n Everything saved for Figure number {ind_range} \033[0m \n")
            except Exception as e:
                print(f"\033[91m \n Figure not plotted for time range {trange} \n because of"
                      f"following exception: {e} \n \033[0m")
        # except Exception as e:
        #     # Print the error in green
        #     print("\033[92m", f"Figure not plotted for {trange} and index value of {ind_range}\n",
        #           "\033[0m")
        #     continue

print(f"Took {round(time.time() - start, 3)} seconds")
