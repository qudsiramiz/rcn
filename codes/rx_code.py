# This the python version of IDL code named 'RX_model_batch.pro'
import datetime
import importlib
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#from rx_model_funcs import rx_model, ridge_finder_multiple
import rx_model_funcs as rmf

importlib.reload(rmf)

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

start = time.time()

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

# trange_list = [
# ['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
# ['2016-12-07 05:11:00', '2016-12-07 05:21:00'],
# #['2015-09-08 11:05:00', '2015-09-08 11:15:00'],
# #['2015-09-19 07:43:30'],
# ['2015-10-16 10:33:30'],
# ['2015-10-16 13:07:02'],
# ['2015-10-22 06:05:22'],
# ['2015-11-01 15:08:06'],
# ['2015-11-12 07:19:21'],
# ['2015-12-06 23:38:31'],
# ['2015-12-08 11:20:44'],
# ['2015-12-09 01:06:11'],
# ['2015-12-14 01:17:40'],
# ['2016-01-07 09:36:15'],
# ['2016-01-10 09:13:37'],
# ['2016-10-22 12:58:41'],
# ['2016-11-02 14:46:18'],
# ['2016-11-06 08:40:58'],
# ['2016-11-12 17:48:47'],
# ['2016-11-13 09:10:41'],
# ['2016-11-18 12:08:11'],
# #['2016-11-23 07:49:33'],
# #['2016-11-23 07:49:52'],
# ['2016-11-23 07:50:30'],
# ['2016-11-28 15:47:00'],
# ['2016-12-11 04:41:50'],
# ['2016-12-19 14:15:02'],
# ['2017-01-02 02:58:13'],
# ['2017-01-11 04:22:43'],
# ['2017-01-20 12:32:07'],
# ['2017-01-22 10:15:46'],
# #['2017-01-22 10:15:58'],
# ['2017-01-22 10:47:33'],
# ['2017-01-27 12:05:43'],
# ['2015-09-25 12:05:00', '2015-09-25 12:13:00'],
# ['2015-09-11 15:23:00'],
# ['2016-01-17 06:29:40'],
# ['2015-09-19 07:41:38'],
# ['2015-09-18 11:54:02', '2015-09-18 11:55:59']
# ]

# trange_list = [['2015-09-11 15:18:00', '2015-09-11 15:28:00']]

# trange_list = [['2017-01-22 10:47:33']]
# trange_list = [['2016-12-24 15:08:00', '2016-12-24 15:12:00']]
# trange_list = [['2017-01-22 10:47:33'],
#                #['2015-12-06 23:38:31'],
#                ['2016-12-24 15:08:00', '2016-12-24 15:12:00']
#                ]
#  Sort the trange_list by the start time
# trange_list.sort(key=lambda x: x[0])
trange_ind_list = np.array([2273, 70, 75, 104, 2263, 2259, 2257, 2142, 2052, 2053])

df_jet_reversal = pd.read_csv("../data/mms3_jet_reversal_times_list.csv")
# Set the index to Date in UTC
df_jet_reversal.set_index('Date', inplace=True)
# Sort the dataframe by the index
df_jet_reversal.sort_index(inplace=True)
# Convert the index to datetime
df_jet_reversal.index = pd.to_datetime(df_jet_reversal.index)

trange_list = df_jet_reversal.index.tolist()
#trange_list_new = trange_list[trange_ind_list]
mms_probe_num_list = [1, 2, 3, 4]
ind_min = trange_ind_list[9]
ind_max = ind_min + 1
for mms_probe_num in mms_probe_num_list[2:3]:
    for ind_range, trange in enumerate(trange_list[ind_min:ind_max], start=ind_min):
        # Convert trange to string to format '%Y-%m-%d %H:%M:%S'
        trange = trange.strftime('%Y-%m-%d %H:%M:%S')
        trange = [trange.split("+")[0].split(".")[0]]
        trange = ["2015-12-06 23:38:31"]
        print(trange)
        try:
            for something in range(1):
                mms_probe_num = str(mms_probe_num)
                min_max_val = 20
                dr = 0.25
                y_min = - min_max_val
                y_max = min_max_val
                z_min = - min_max_val
                z_max = min_max_val
                model_type = 't96'

                model_inputs = {
                    'trange': trange,
                    "probe": None,
                    "omni_level": 'hro',
                    "mms_probe_num": mms_probe_num,
                    "model_type": model_type,
                    "m_p": 0.5,
                    "dr": dr,
                    "min_max_val": min_max_val,
                    "y_min": y_min,
                    "y_max": y_max,
                    "z_min": z_min,
                    "z_max": z_max,
                    "save_data": False,
                }
                (bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh, sw_params, x_shu, y_shu, z_shu,
                 b_msx, b_msy, b_msz) = rmf.rx_model(**model_inputs)

                # mask = shear > 175
                # shear[mask] = 0
                # Normalize each quantity to the range [0, 1]
                shear_norm = (shear - np.nanmin(shear)) / (np.nanmax(shear) - np.nanmin(shear))
                rx_en_norm = (rx_en - np.nanmin(rx_en)) / (np.nanmax(rx_en) - np.nanmin(rx_en))
                va_cs_norm = (va_cs - np.nanmin(va_cs)) / (np.nanmax(va_cs) - np.nanmin(va_cs))
                bisec_msp_norm = (bisec_msp - np.nanmin(bisec_msp)) / (np.nanmax(bisec_msp) - np.nanmin(bisec_msp))
                
                figure_inputs = {
                    "image": [shear_norm, rx_en_norm, va_cs_norm, bisec_msp_norm],
                    "convolution_order": [1, 1, 1, 1],
                    "t_range": trange,
                    "b_imf": np.round(sw_params['b_imf'], 2),
                    "b_msh": np.round(sw_params['mms_b_gsm'], 2),
                    "xrange": [y_min, y_max],
                    "yrange": [z_min, z_max],
                    "mms_probe_num": mms_probe_num,
                    "mms_sc_pos": np.round(np.nanmean(sw_params['mms_sc_pos'], axis=0), 2),
                    "dr": dr,
                    "dipole_tilt_angle": sw_params['ps'],
                    "imf_clock_angle": sw_params['imf_clock_angle'],
                    "sigma": [2, 2, 2, 2],
                    "mode": "nearest",
                    "alpha": 1,
                    "vmin": [0, 0, 0, 0],
                    "vmax": [1, 1, 1, 1],
                    "cmap_list": ["viridis", "cividis", "plasma", "magma"],
                    "draw_patch": [True, True, True, True],
                    "draw_ridge": [True, True, True, True],
                    "save_fig": True,
                    "fig_name": 'crossing_all_ridge_plots',
                    # "fig_format": 'png',
                    "c_label": ['Shear', 'Reconnection Energy', 'Exhaust Velocity', 'Bisection Field'],
                    #"c_unit": [r'${}^\circ$', 'nPa', 'km/s', 'nT'],
                    "c_unit": ['', 'nPa', 'km/s', 'nT'],
                    "wspace": 0.0,
                    "hspace": 0.17,
                    "fig_size": (8.775, 10),
                    "box_style": dict(boxstyle='round', color='k', alpha=0.8),
                    # "box_style": dict(boxstyle='round', color=None, alpha=0.8),
                    "title_y_pos": 1.09,
                    "interpolation": 'None',
                    "tsy_model": model_type,
                    "dark_mode": False,
                    "rc_file_name": f"reconnection_line_data_mms{mms_probe_num}_20220612.csv",
                    "rc_folder": "../data/rx_d/",
                    "save_rc_file": False,
                    "walen1": df_jet_reversal['walen1'][ind_range],
                    "walen2": df_jet_reversal["walen2"][ind_range],
                    "jet_detection": df_jet_reversal['jet_detection'][ind_range],
                    "fig_version": 'v07',
                }

                y_vals, x_intr_vals_list, y_intr_vals_list = rmf.ridge_finder_multiple(**figure_inputs,
                                                                                       fig_format='pdf')
            print(f"\033[92m \n Everything saved for Figure number {ind_range} \033[0m \n")
        except Exception as e:
            print(f"\033[91m \n Figure not plotted for time range {trange} \n because of following exception: {e} \n \033[0m")
        #except Exception as e:
        #    # Print the error in green
        #    print('\033[92m', f'Figure not plotted for {trange} and index value of {ind_range}\n',
        #          '\033[0m')
        #    continue

print(f'Took {round(time.time() - start, 3)} seconds')
