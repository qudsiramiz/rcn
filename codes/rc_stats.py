import importlib
import glob
import numpy as np
import pandas as pd

import rc_stats_fncs as rcsf
importlib.reload(rcsf)

'''
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

'''

'''
data_folder = '../data/rx_d'
fname = np.sort(glob.glob(f"{data_folder}/reconnection_line_data_mms3_20221109.csv"))[0]
df = pd.read_csv(fname, index_col=False)
# Set date_from as index
df = df.set_index("date_from")

df_shear = df[df.method_used == "shear"]
df_rx_en = df[df.method_used == "rx_en"]
df_va_cs = df[df.method_used == "va_cs"]
df_bisec = df[df.method_used == "bisection"]

cone_angle_shear = np.arccos(df_shear.b_imf_x / np.sqrt(
                    df_shear.b_imf_x**2 +df_shear.b_imf_y**2+ df_shear.b_imf_z**2)) * 180 / np.pi
cone_angle_rx_en = np.arccos(df_rx_en.b_imf_x / np.sqrt(
                        df_rx_en.b_imf_x**2 +df_rx_en.b_imf_y**2+ df_rx_en.b_imf_z**2)) * 180 / np.pi
cone_angle_va_cs = np.arccos(df_va_cs.b_imf_x / np.sqrt(
                        df_va_cs.b_imf_x**2 +df_va_cs.b_imf_y**2+ df_va_cs.b_imf_z**2)) * 180 / np.pi
cone_angle_bisec = np.arccos(df_bisec.b_imf_x / np.sqrt(
                        df_bisec.b_imf_x**2 +df_bisec.b_imf_y**2+ df_bisec.b_imf_z**2)) * 180 / np.pi

b_imf_mag_shear = np.sqrt(df_shear.b_imf_x**2 +df_shear.b_imf_y**2+ df_shear.b_imf_z**2)
b_imf_mag_rx_en = np.sqrt(df_rx_en.b_imf_x**2 +df_rx_en.b_imf_y**2+ df_rx_en.b_imf_z**2)
b_imf_mag_va_cs = np.sqrt(df_va_cs.b_imf_x**2 +df_va_cs.b_imf_y**2+ df_va_cs.b_imf_z**2)
b_imf_mag_bisec = np.sqrt(df_bisec.b_imf_x**2 +df_bisec.b_imf_y**2+ df_bisec.b_imf_z**2)

df_shear["cone_angle"] = cone_angle_shear
df_rx_en["cone_angle"] = cone_angle_rx_en
df_va_cs["cone_angle"] = cone_angle_va_cs
df_bisec["cone_angle"] = cone_angle_bisec

df_shear["b_imf_mag"] = b_imf_mag_shear
df_rx_en["b_imf_mag"] = b_imf_mag_rx_en
df_va_cs["b_imf_mag"] = b_imf_mag_va_cs
df_bisec["b_imf_mag"] = b_imf_mag_bisec

'''

# Select all values with 40 < cone_angle < 120 and |Bz| > 0 nT
# df_shear_lim = df_shear[(df_shear.cone_angle > 40) & (df_shear.cone_angle < 120) & (df_shear.b_imf_z > 0)]
# df_rx_en_lim = df_rx_en[(df_rx_en.cone_angle > 40) & (df_rx_en.cone_angle < 120) & (df_rx_en.b_imf_z > 0)]
# df_va_cs_lim = df_va_cs[(df_va_cs.cone_angle > 40) & (df_va_cs.cone_angle < 120) & (df_va_cs.b_imf_z > 0)]
# df_bisec_lim = df_bisec[(df_bisec.cone_angle > 40) & (df_bisec.cone_angle < 120) & (df_bisec.b_imf_z > 0)]

# Select all values with cone angle > 120
df_shear_lim = df_shear[(df_shear.b_imf_y.abs() / df_shear.b_imf_mag > 0.7) & (df_shear.b_imf_z > 0)]
df_rx_en_lim = df_rx_en[(df_rx_en.b_imf_y.abs() / df_rx_en.b_imf_mag > 0.7) & (df_rx_en.b_imf_z > 0)]
df_va_cs_lim = df_va_cs[(df_va_cs.b_imf_y.abs() / df_va_cs.b_imf_mag > 0.7) & (df_va_cs.b_imf_z > 0)]
df_bisec_lim = df_bisec[(df_bisec.b_imf_y.abs() / df_bisec.b_imf_mag > 0.7) & (df_bisec.b_imf_z > 0)]

# Print the mean and median of the 'r_rc' column in a table
print(f"{'Method':<15}{'Mean':<15}{'Median':<15}")
print(f"{'Shear':<15}{df_shear_lim.r_rc.mean():<15.2f}{df_shear_lim.r_rc.median():<15.2f}")
print(f"{'Rx_en':<15}{df_rx_en_lim.r_rc.mean():<15.2f}{df_rx_en_lim.r_rc.median():<15.2f}")
print(f"{'Va_cs':<15}{df_va_cs_lim.r_rc.mean():<15.2f}{df_va_cs_lim.r_rc.median():<15.2f}")
print(f"{'Bisec':<15}{df_bisec_lim.r_rc.mean():<15.2f}{df_bisec_lim.r_rc.median():<15.2f}")

