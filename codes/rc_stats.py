import importlib
import glob
import numpy as np

import rc_stats_fncs as rcsf
importlib.reload(rcsf)

data_folder = '../data/rx_d'
fnames = np.sort(glob.glob(f"{data_folder}/*_20220612.csv"))
cut_type_list = ["jet", "walen1", "walen2", "walen_jet"]
for file_name in fnames:
    for cut_type in cut_type_list[:]:
        mms_probe_num = file_name.split('/')[-1].split('_')[-1].split('.')[0]

        fig_inputs ={
            'nbins': 15,
            'file_name': file_name,
            'dark_mode': False,
            'fig_name':  f"rx_hist_{mms_probe_num}_",
            'fig_format': 'pdf',
            'fig_folder': '../figures/rx_hist_v07',
            'fig_size': (8, 8),
            'histtype': 'step',
            'linewidth': 3,
            'cut_type': cut_type,
        }

        df_shear, df_rx_en, df_va_cs, df_bisec = rcsf.plot_hist(**fig_inputs)
