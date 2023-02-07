import importlib
import glob
import numpy as np

import rc_stats_fncs as rcsf
importlib.reload(rcsf)

data_folder = '../data/rx_d'
fnames = np.sort(glob.glob(f"{data_folder}/reconnection_line_data_mms3_20221109.csv"))
# cut_type_list = ["jet", "walen1", "walen2", "walen_jet"]
cut_type_list = ["bz_neg", "bz_pos", 'bz', 'cone_angle', 'cone_and_bz_neg']
for file_name in fnames:
    for cut_type in cut_type_list[:]:
        mms_probe_num = file_name.split('/')[-1].split('_')[-1].split('.')[0]
        dark_mode = True
        fig_inputs = {
            'bins': np.linspace(0, 20, 25),
            'file_name': file_name,
            'dark_mode': dark_mode,
            'fig_name':  f"rx_hist_{mms_probe_num}_{dark_mode}",
            'fig_format': 'jpg',
            'fig_folder': '../figures/rx_hist/rx_hist_v16',
            'fig_size': (8, 8),
            'histtype': 'step',
            'linewidth': 3,
            'cut_type': cut_type,
            'r_lim': [0.01, 20],
            'density': True,
        }

        df_shear, df_rx_en, df_va_cs, df_bisec = rcsf.plot_hist(**fig_inputs)
