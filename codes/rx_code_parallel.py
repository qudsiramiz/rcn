# This the python version of IDL code named "RX_model_batch.pro"
import datetime
import importlib
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import joblib as jl

#from rx_model_funcs import rx_model, ridge_finder_multiple
import rx_model_funcs as rmf

importlib.reload(rmf)

# Set the fontstyle to Times New Roman
font = {"family": "serif", "weight": "normal", "size": 10}
plt.rc("font", **font)
plt.rc("text", usetex=True)


def rx_model_parallel_run(mms_probe_num, trange):

    mms_probe_num = str(mms_probe_num)
    min_max_val = 20
    dr = 0.25
    y_min = - min_max_val
    y_max = min_max_val
    z_min = - min_max_val
    z_max = min_max_val
    model_type = "t96"

    model_inputs = {
        "trange": trange,
        "probe": None,
        "omni_level": "hro",
        "mms_probe_num": mms_probe_num,
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
    _, _, _, shear, rx_en, va_cs, bisec_msp, _, sw_params, _, _, _, _, _, _ = rmf.rx_model(
                                                                                     **model_inputs)

    shear_norm = (shear - np.nanmin(shear)) / (np.nanmax(shear) - np.nanmin(shear))
    rx_en_norm = (rx_en - np.nanmin(rx_en)) / (np.nanmax(rx_en) - np.nanmin(rx_en))
    va_cs_norm = (va_cs - np.nanmin(va_cs)) / (np.nanmax(va_cs) - np.nanmin(va_cs))
    bisec_msp_norm = (bisec_msp - np.nanmin(bisec_msp)) / (np.nanmax(bisec_msp) -
                                                           np.nanmin(bisec_msp))

    return shear_norm, rx_en_norm, va_cs_norm, bisec_msp_norm, sw_params, model_inputs


def make_ridge_figures(mms_probe_num, trange, ind_range, df_jet_reversal):
    trange = trange.strftime('%Y-%m-%d %H:%M:%S')
    trange = [trange.split("+")[0].split(".")[0]]
    (shear_norm, rx_en_norm, va_cs_norm, bisec_msp_norm, sw_params,
                                        model_inputs) = rx_model_parallel_run(mms_probe_num, trange)
    figure_inputs = {
        "image": [shear_norm, rx_en_norm, va_cs_norm, bisec_msp_norm],
        "convolution_order": [1, 1, 1, 1],
        "t_range": trange,
        "b_imf": np.round(sw_params["b_imf"], 2),
        "b_msh": np.round(sw_params["mms_b_gsm"], 2),
        "xrange": [model_inputs["y_min"], model_inputs["y_max"]],
        "yrange": [model_inputs["z_min"], model_inputs["z_max"]],
        "mms_probe_num": mms_probe_num,
        "mms_sc_pos": np.round(np.nanmean(sw_params["mms_sc_pos"], axis=0), 2),
        "dr": model_inputs["dr"],
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
        "c_label": ["Shear", "Reconnection Energy", "Exhaust Velocity", "Bisection Field"],
        #"c_unit": [r"${}^\circ$", "nPa", "km/s", "nT"],
        "c_unit": ["", "nPa", "km/s", "nT"],
        "wspace": 0.0,
        "hspace": 0.17,
        "fig_size": (8.775, 10),
        "box_style": dict(boxstyle="round", color="k", alpha=0.8),
        # "box_style": dict(boxstyle="round", color=None, alpha=0.8),
        "title_y_pos": 1.09,
        "interpolation": "None",
        "tsy_model": model_inputs["model_type"],
        "dark_mode": True,
        "rc_file_name": f"reconnection_line_data_mms{mms_probe_num}_20220927.csv",
        "rc_folder": "../data/rx_d/",
        "save_rc_file": True,
        "walen1": df_jet_reversal["walen1"][ind_range],
        "walen2": df_jet_reversal["walen2"][ind_range],
        "jet_detection": df_jet_reversal["jet_detection"][ind_range],
        "fig_version": "v11",
        "r_W": df_jet_reversal["r_W"][ind_range],
        "theta_W": df_jet_reversal["theta_w"][ind_range],
        #"jet_time": df_jet_reversal["jet_time"][ind_range],
        "np_median_msp": df_jet_reversal["np_msp_median"][ind_range],
        "np_median_msh": df_jet_reversal["np_msh_median"][ind_range],
        "df_jet_reversal": df_jet_reversal.iloc[ind_range],
    }

    y_vals, x_intr_vals_list, y_intr_vals_list = rmf.ridge_finder_multiple(**figure_inputs,
                                                                                   fig_format="png")
    print(f"\033[92m \n Everything saved for Figure number {ind_range} \033[0m \n")
        #except Exception as e:
        #    print(f"\033[91m \n Figure not plotted for time range {trange} \n because of following exception: {e} \n \033[0m")
        #except Exception as e:
        #    # Print the error in green
        #    print("\033[92m", f"Figure not plotted for {trange} and index value of {ind_range}\n",
        #          "\033[0m")
        #    continue

start = time.time()

today_date = datetime.datetime.today().strftime("%Y-%m-%d")

# Read the file with jet reversal times
df_jet_reversal = pd.read_csv("../data/mms_jet_reversal_times_list_20220922_beta.csv")
# Set the index to Date in UTC
df_jet_reversal.set_index("jet_time", inplace=True)
# Sort the dataframe by the index
df_jet_reversal.sort_index(inplace=True)
# Convert the index to datetime
df_jet_reversal.index = pd.to_datetime(df_jet_reversal.index)

trange_list = df_jet_reversal.index.tolist()

mms_probe_num_list = [3]
#mms_probe_num = 3
ind_min = 0
ind_max = 5
#ind_range = range(ind_min, ind_max)

use_parallel = True

if use_parallel:
    # Parallel run
    jl.Parallel(n_jobs=4, backend="multiprocessing")(
        jl.delayed(make_ridge_figures)(mms_probe_num, trange_list[ind_range], ind_range, df_jet_reversal)
        for ind_range in range(ind_min, ind_max)
        for mms_probe_num in mms_probe_num_list
    )
else:
    # Serial run
    for ind_range in range(ind_min, ind_max):
        for mms_probe_num in mms_probe_num_list:
            make_ridge_figures(mms_probe_num, trange_list[ind_range], ind_range, df_jet_reversal)

print(f"Took {round(time.time() - start, 3)} seconds")