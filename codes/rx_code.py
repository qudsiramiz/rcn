# This the python version of IDL code named 'RX_model_batch.pro'
import datetime
import time

import matplotlib.pyplot as plt
import numpy as np

from rx_model_funcs import *

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

start = time.time()

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

trange_list = [
['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
['2016-12-07 05:05:00', '2016-12-07 05:33:00'],
#['2015-09-08 11:05:00', '2015-09-08 11:15:00'],
#['2015-09-19 07:43:30'],
['2015-10-16 10:33:30'],
['2015-10-16 13:07:02'],
['2015-10-22 06:05:22'],
['2015-11-01 15:08:06'],
['2015-11-12 07:19:21'],
['2015-12-06 23:38:31'],
['2015-12-08 11:20:44'],
['2015-12-09 01:06:11'],
['2015-12-14 01:17:40'],
['2016-01-07 09:36:15'],
['2016-01-10 09:13:37'],
['2016-10-22 12:58:41'],
['2016-11-02 14:46:18'],
['2016-11-06 08:40:58'],
['2016-11-12 17:48:47'],
['2016-11-13 09:10:41'],
['2016-11-18 12:08:11'],
#['2016-11-23 07:49:33'],
#['2016-11-23 07:49:52'],
['2016-11-23 07:50:30'],
['2016-11-28 15:47:00'],
['2016-12-11 04:41:50'],
['2016-12-19 14:15:02'],
['2017-01-02 02:58:13'],
['2017-01-11 04:22:43'],
['2017-01-20 12:32:07'],
['2017-01-22 10:15:46'],
#['2017-01-22 10:15:58'],
['2017-01-22 10:47:33'],
['2017-01-27 12:05:43'],
]
count = 0
for trange in trange_list[0:]:

    mms_probe_num = '1'
    min_max_val = 15
    dr = 0.25
    y_min = - min_max_val
    y_max = min_max_val
    z_min = - min_max_val
    z_max = min_max_val

    model_inputs = {
        'trange': trange,
        "probe" : None,
        "omni_level" : 'hro',
        "mms_probe_num" : mms_probe_num,
        "model_type" : 't96',
        "m_p" : 0.5,
        "dr" : dr,
        "min_max_val" : min_max_val,
        "y_min" : y_min,
        "y_max" : y_max,
        "z_min" : z_min,
        "z_max" : z_max,
        "save_data" : False,
        }
    bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh, sw_params = rx_model(**model_inputs)


    figure_inputs = {
        "image" : [shear, rx_en/np.nanmax(rx_en), va_cs, bisec_msp],
        "t_range" : trange,
        "xrange" : [y_min, y_max],
        "yrange" : [z_min, z_max],
        "mms_probe_num" : mms_probe_num,
        "mms_sc_pos" : [np.nanmean(sw_params['mms_sc_pos'][:,1]),
                        np.nanmean(sw_params['mms_sc_pos'][:,2])],
        "dr" : dr,
        "dipole_tilt_angle" : sw_params['ps'],
        "imf_clock_angle" : sw_params['imf_clock_angle'],
        "sigma" : [2, 2, 2, 2],
        "mode" : "nearest",
        "alpha" : 1,
        "vmin" : [0, 0, None, None],
        "vmax" : [180, 1, None, None],
        "cmap_list" : ["viridis", "cividis", "plasma", "magma"],
        "draw_patch" : [True, True, True, True],
        "draw_ridge" : [True, True, True, True],
        "save_fig" : True,
        "fig_name" : f'all_ridge_plots',
        #"fig_format" : 'png',
        "c_label" : ['Shear', 'Reconnection Energy', 'Exhaust Velocity', 'Bisection Field'],
        "c_unit" : [r'${}^\circ$', 'nPa', 'km/s', 'nT'],
        "wspace" : 0.0,
        "hspace" : 0.17,
        "fig_size" : (8.775, 10),
        "box_style": dict(boxstyle='round', facecolor='black', alpha=0.8),
        "title_y_pos" : 1.07,
        "interpolation" : 'gaussian',
    }

    ridge_finder_multiple(**figure_inputs, fig_format='png')
    #print(f"Model run for date {trange[0]} to {trange[1]}")
    #ridge_finder_multiple(**figure_inputs, fig_format='pdf')
    count += 1
    '''
    # Check if 'plot_type' has length attribute. If it has length attribute then plot the ridge plot
    # for each of the plot type in the list. If it does not have length attribute then plot the
    # ridge plot for the specified plot type.
    plot_type = 'all'
    types_of_plot = ['shear', 'rx_en', 'va_cs', 'Bisec-msp', 'Bisec-msh', 'all']

    common_figure_inputs = {
        "t_range" : trange,
        "xrange" : [y_min, y_max],
        "yrange" : [z_min, z_max],
        "sigma" : 2.2,
        "dr" : dr,
        "dipole_tilt_angle" : sw_params['ps'],
        "imf_clock_angle" : sw_params['imf_clock_angle'],
        "draw_patch" : True,
        "draw_ridge" : True,
        "mms_probe_num" : mms_probe_num,
        "mms_sc_pos" : [np.nanmean(sw_params['mms_sc_pos'][:,1]),
                        np.nanmean(sw_params['mms_sc_pos'][:,2])],
        "fig_format" : 'png',
    }

    shear_figure_inputs = {
        "image" : shear,
        "fig_name" : "shear",
        "c_label" : "Shear",
        "c_unit" : r'${}^\circ$',
        "cmap" : "cividis",
        }

    rx_en_figure_inputs = {
        "image" : rx_en/np.nanmax(rx_en),
        "fig_name" : "rx-en_nPa",
        "c_label" : "Reconnection Energy",
        "c_unit" : "nPa",
        "cmap" : "viridis"
    }
    va_cs_figure_inputs = {
        "image" : va_cs,
        "fig_name" : "va-cs",
        "c_label" : "Exhaust Velocity",
        "c_unit" : "Km/s",
        "cmap" : "plasma"
    }
    bisec_msp_figure_inputs = {
        "image" : bisec_msp,
        "fig_name" : "bisec_msp",
        "c_label" : "Bisection Field",
        "c_unit" : "nT",
        "cmap" : "inferno"
    }
    bisec_msh_figure_inputs = {
        "image" : bisec_msh,
        "fig_name" : "bisec_msh",
        "c_label" : "Bisection Field",
        "c_unit" : "nT",
        "cmap" : "magma"
    }
    if isinstance(plot_type, list):
        for xx in plot_type:
            if xx not in types_of_plot:
                raise ValueError(
                    f'{xx} is not a valid plot type. Please choose from {types_of_plot}'
                    )
        if 'shear' in plot_type:
            print('Plotting shear')
            _ = ridge_finder(**shear_figure_inputs, **common_figure_inputs)
        if 'rx_en' in plot_type:
            print('Plotting rx_en')
            _ = ridge_finder(**rx_en_figure_inputs, **common_figure_inputs)
        if 'va_cs' in plot_type:
            print('Plotting va_cs')
            _ = ridge_finder(**va_cs_figure_inputs, **common_figure_inputs)
        if 'Bisec-msp' in plot_type:
            print('Plotting Bisec-msp')
            _ = ridge_finder(**bisec_msp_figure_inputs, **common_figure_inputs)
        if 'Bisec-msh' in plot_type:
            print('Plotting Bisec-msh')
            _ = ridge_finder(**bisec_msh_figure_inputs, **common_figure_inputs)
    elif plot_type == 'all':
        print('Plotting for all')
        _ = ridge_finder(**shear_figure_inputs, **common_figure_inputs)

        _ = ridge_finder(**rx_en_figure_inputs, **common_figure_inputs)

        _ = ridge_finder(**va_cs_figure_inputs, **common_figure_inputs)

        _ = ridge_finder(**bisec_msp_figure_inputs, **common_figure_inputs)

        _ = ridge_finder(**bisec_msh_figure_inputs, **common_figure_inputs)
    elif plot_type == 'shear':
        print(trange)
        print('Plotting shear')
        _ = ridge_finder(**shear_figure_inputs, **common_figure_inputs)
    elif plot_type == 'rx_en':
        print('Plotting rx_en')
        _ = ridge_finder(**rx_en_figure_inputs, **common_figure_inputs)
    elif plot_type == 'va_cs':
        print('Plotting va_cs')
        _ = ridge_finder(**va_cs_figure_inputs, **common_figure_inputs)
    elif plot_type == 'Bisec-msp':
        print('Plotting Bisec-msp')
        _ = ridge_finder(**bisec_msp_figure_inputs, **common_figure_inputs)

    elif plot_type == 'Bisec-msh':
        print('Plotting Bisec-msh')
        _ = ridge_finder(**bisec_msh_figure_inputs, **common_figure_inputs)
    else:
        raise KeyError(
            'plot_type must be one or a list of: all, shear, rx-en, va-cs, Bisec-msp, Bisec-msh')

    #_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=b_msy, bz=b_msz, save_fig=True,
    #                      scale=40, fig_name="magnetosheath")
    #_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=by, bz=bz, save_fig=True, scale=120,
    #                 fig_name="magnetosphere")
    '''
print(f'Took {round(time.time() - start, 3)} seconds')
