# This the python version of IDL code named "RX_model_batch.pro"
import datetime
import importlib
import os
import time
from contextlib import contextmanager, redirect_stderr, redirect_stdout

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import rx_model_funcs_sim as rmf

importlib.reload(rmf)

# Set the fontstyle to Times New Roman
font = {"family": "serif", "weight": "normal", "size": 10}
plt.rc("font", **font)
plt.rc("text", usetex=True)

start = time.time()

today_date = datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

plt.close("all")

@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


b_imf_list = np.array([np.array([-2, 2, 0]), np.array([0, 0, 5]), np.array([0, 0, -5])])

m_proton = 1.672e-27  # Mass of proton in SI unit

np_imf = 5
v_imf = np.array([-500, 0, 0])
sym_h_imf = - 5
tp_imf = 1.5E6

rho = np_imf * m_proton * 1.15

p_dyn = 1.6726e-6 * 1.15 * np_imf * (v_imf[0]**2 + v_imf[1]**2 + v_imf[2]**2)

for b_imf in b_imf_list[1:]:
    param = [p_dyn, sym_h_imf, b_imf[1], b_imf[2], 0, 0, 0, 0, 0, 0]

    imf_clock_angle = np.arctan2(b_imf[1], b_imf[2]) * 180 / np.pi
    if imf_clock_angle < 0:
        imf_clock_angle += 360

    # Make a dictionary of all the solar wind parameters
    sw_params = {}
    sw_params['time'] = today_date
    sw_params['b_imf'] = b_imf
    sw_params['rho'] = rho
    sw_params['ps'] = 0.0
    sw_params['p_dyn'] = p_dyn
    sw_params['sym_h'] = sym_h_imf
    sw_params['t_p'] = tp_imf
    sw_params['imf_clock_angle'] = imf_clock_angle
    sw_params['param'] = param
    sw_params['mms_time'] = np.nan
    sw_params['mms_sc_pos'] = [np.nan, np.nan, np.nan]
    sw_params['mms_b_gsm'] = [np.nan, np.nan, np.nan]
    sw_params['mms_v_gsm'] = [np.nan, np.nan, np.nan]

    # with suppress_stdout_stderr():
    for foo in range(1):
        for bar in range(1):
        # try:
            mms_probe_num = "mms3"
            min_max_val = 20
            dr = 0.25
            y_min = - min_max_val
            y_max = min_max_val
            z_min = - min_max_val
            z_max = min_max_val
            model_type = "t96"

            model_inputs = {
                "trange": today_date,
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
                "sw_params": sw_params,
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
                "t_range": today_date,
                "b_imf": np.round(sw_params["b_imf"], 2),
                "b_msh": None,
                "xrange": [y_min, y_max],
                "yrange": [z_min, z_max],
                "mms_probe_num": mms_probe_num,
                "mms_sc_pos": sw_params["mms_sc_pos"],
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
                "c_label": ["Shear", "Reconnection Energy", "Exhaust Velocity",
                            "Bisection Field"],
                "c_unit": ["", "nPa", "km/s", "nT"],
                "wspace": 0.0,
                "hspace": 0.17,
                "fig_size": (12, 3),
                "box_style": dict(boxstyle="round", color="k", alpha=0.8),
                "title_y_pos": 1.15,
                "interpolation": "None",
                "tsy_model": model_type,
                "dark_mode": False,
                "rc_file_name": f"reconnection_line_data_mms{mms_probe_num}_20230330.csv",
                "rc_folder": "../data/rx_d/",
                "save_rc_file": True,
                "fig_version": "v_test_20230330",
            }

            y_vals, x_intr_vals_list, y_intr_vals_list = rmf.ridge_finder_multiple(
                                                        **figure_inputs, fig_format="pdf")
            print(f"\033[92m \n Everything saved for Figure number \033[0m \n")

        # except Exception as e:
        #     print(f"\033[91m \n Figure not plotted for time range {trange} \n because of"
        #             f"following exception: {e} \n \033[0m")
            # except Exception as e:
            #     # Print the error in green
            #     print("\033[92m", f"Figure not plotted for {trange} and index value of {ind_range}\n",
            #           "\033[0m")
            #     continue

print(f"Took {round(time.time() - start, 3)} seconds")
