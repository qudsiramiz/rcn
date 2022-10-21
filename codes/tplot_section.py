    if jet_detection:

    tplot_globat_options = {"show_all_axes": True,
                            "black_background": True,
                            "crosshair": True,
                            "vertical_spacing": 0,
                            }
    for key in tplot_globat_options:
        ptt.tplot_options(key, tplot_globat_options[key])


    # ptt.timebar(ptt.get_data(mms_fpi_varnames[0])[0].min() + 200, color='red', dash=True, thick=2)
    # ptt.timespan(ptt.get_data(mms_fpi_varnames[0])[0].min() + 160, 100, keyword="seconds")
    energy_spectr_dict_option = {'Colormap': "Spectral_r",
                                 'ylog': True,
                                 'zlog': True}
    temp_dict_option = {'ylog': True,
                        'color': ['red', 'blue'],
                        'linestyle': '-',
                        'legend_names': ['para', 'perp']
                        }
    delta_v_min_max_dict_option = {'color': ['k', 'k'],
                                   'linestyle': '-',
                                   'legend_names': ['delta_v_min', 'delta_v_max'],
                                    'yrange': [-200, 200]
                                   }
    keys_to_plot = [f'mms{probe}_dis_numberdensity_{data_rate}',
                    'Tp',
                    f'mms{probe}_dis_energyspectr_omni_{data_rate}',
                    f'mms{probe}_des_energyspectr_omni_{data_rate}',
                    f'mms{probe}_dis_bulkv_lmn_{data_rate}',
                    'delta_v_min_max',
                    'theta_w_deg',
                    'R_w',
                    ]
    #ptt.options(f'mms{probe}_dis_numberdensity_{data_rate}', 'ylog', True)
    ptt.options(f'mms{probe}_dis_energyspectr_omni_{data_rate}', opt_dict=energy_spectr_dict_option)
    ptt.options(f'mms{probe}_des_energyspectr_omni_{data_rate}', opt_dict=energy_spectr_dict_option)
    ptt.options('Tp', opt_dict=temp_dict_option)
    ptt.options('delta_v_min_max', opt_dict=delta_v_min_max_dict_option)

    ptt.tplot(keys_to_plot, save_png="test_plot", display=False)
