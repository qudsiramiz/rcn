import importlib
from tkinter.tix import Tree
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import rx_model_funcs as rxmf
import seaborn_plots_fncs as spf
import SeabornFig2Grid as sfg

importlib.reload(rxmf)
importlib.reload(spf)
importlib.reload(sfg)

'''
file_name = "../data/rx_d/reconnection_line_data_mms3_20221027.csv"
cut_type_list = ["jet", "walen1", "walen2", "walen_jet"]

df = pd.read_csv(file_name, index_col=False)

df_shear = df[df.method_used == "shear"]
df_rx_en = df[df.method_used == "rx_en"]
df_va_cs = df[df.method_used == "va_cs"]
df_bisec = df[df.method_used == "bisection"]

df_shear_imfz = df_shear[df_shear["b_imf_z"] < 0]
df_rx_en_imfz = df_rx_en[df_rx_en["b_imf_z"] < 0]
df_va_cs_imfz = df_va_cs[df_va_cs["b_imf_z"] < 0]
df_bisec_imfz = df_bisec[df_bisec["b_imf_z"] < 0]

# If clock angle is greater than 180, then subtract 360 to get the angle
# between the IMF and the Bz field.
df_shear["imf_clock_angle"] = df_shear["imf_clock_angle"].apply(
    lambda x: 360 - x if x > 180 else x)
df_rx_en["imf_clock_angle"] = df_rx_en["imf_clock_angle"].apply(
    lambda x: 360 - x if x > 180 else x)
df_va_cs["imf_clock_angle"] = df_va_cs["imf_clock_angle"].apply(
    lambda x: 360 - x if x > 180 else x)
df_bisec["imf_clock_angle"] = df_bisec["imf_clock_angle"].apply(
    lambda x: 360 - x if x > 180 else x)
# plt.figure(figsize=(8,8))

"""
['Date', 'Probe', 'angle_b_lmn_vec_msp_msh_median', 'b_imf_x',
       'b_imf_y', 'b_imf_z', 'b_lmn_vec_msh_mean_l',
       'b_lmn_vec_msh_mean_m', 'b_lmn_vec_msh_mean_n',
       'b_lmn_vec_msh_median_l', 'b_lmn_vec_msh_median_m',
       'b_lmn_vec_msh_median_n', 'b_lmn_vec_msp_mean_l',
       'b_lmn_vec_msp_mean_m', 'b_lmn_vec_msp_mean_n',
       'b_lmn_vec_msp_median_l', 'b_lmn_vec_msp_median_m',
       'b_lmn_vec_msp_median_n', 'b_msh_x', 'b_msh_y', 'b_msh_z',
       'beta_msh_mean', 'date_from', 'date_to', 'dipole',
       'imf_clock_angle', 'ind_max_msh', 'ind_max_msp', 'ind_min_msh',
       'ind_min_msp', 'jet_detection', 'method_used', 'mms_spc_num',
       'np_median_msh', 'np_median_msp', 'np_msh_mean', 'np_msh_median',
       'np_msp_mean', 'np_msp_median', 'p_dyn', 'r_W', 'r_rc', 'r_spc',
       'spc_pos_x', 'spc_pos_y', 'spc_pos_z', 'theta_W', 'theta_w',
       'tp_para_msh_mean', 'tp_para_msh_median', 'tp_para_msp_mean',
       'tp_para_msp_median', 'tp_perp_msh_mean', 'tp_perp_msh_median',
       'tp_perp_msp_mean', 'tp_perp_msp_median', 'vp_lmn_vec_msh_mean_l',
       'vp_lmn_vec_msh_mean_m', 'vp_lmn_vec_msh_mean_n',
       'vp_lmn_vec_msh_median_l', 'vp_lmn_vec_msh_median_m',
       'vp_lmn_vec_msh_median_n', 'vp_lmn_vec_msp_mean_l',
       'vp_lmn_vec_msp_mean_m', 'vp_lmn_vec_msp_mean_n',
       'vp_lmn_vec_msp_median_l', 'vp_lmn_vec_msp_median_m',
       'vp_lmn_vec_msp_median_n', 'walen1', 'walen2', 'x_gsm', 'y_gsm',
       'z_gsm']
"""

b_lmn_vec_msh = np.full((df_shear.shape[0], 3), np.nan)

b_lmn_vec_msh[:, 0] = df_shear["b_lmn_vec_msh_median_l"]
b_lmn_vec_msh[:, 1] = df_shear["b_lmn_vec_msh_median_m"]
b_lmn_vec_msh[:, 2] = df_shear["b_lmn_vec_msh_median_n"]

b_lmn_vec_msp = np.full((df_shear.shape[0], 3), np.nan)
b_lmn_vec_msp[:, 0] = df_shear["b_lmn_vec_msp_median_l"]
b_lmn_vec_msp[:, 1] = df_shear["b_lmn_vec_msp_median_m"]
b_lmn_vec_msp[:, 2] = df_shear["b_lmn_vec_msp_median_n"]

msh_msp_shear = []
for xx, yy in zip(b_lmn_vec_msh, b_lmn_vec_msp):
    try:
        msh_msp_shear.append(rxmf.get_shear(xx, yy, angle_unit="degrees"))
    except Exception:
        pass


# Define the mass of proton in kg
m_p = 1.6726219e-27

# Define the absolute permeability of free space in m^2 kg^-1 s^-1
mu_0 = 4 * np.pi * 1e-7

# Define the Boltzmann constant in J K^-1
k_B = 1.38064852e-23
b_mag_msh = np.sqrt(df_shear["b_lmn_vec_msh_mean_l"]**2 + df_shear["b_lmn_vec_msh_mean_m"]**2 +
                    df_shear["b_lmn_vec_msh_mean_n"]**2)
b_mag_msp = np.sqrt(df_shear["b_lmn_vec_msp_mean_l"]**2 + df_shear["b_lmn_vec_msp_mean_m"]**2 +
                    df_shear["b_lmn_vec_msp_mean_n"]**2)

# Compute the cone angle
cone_angle = np.arccos(df_shear.b_imf_x / np.sqrt(
                       df_shear.b_imf_x**2 +df_shear.b_imf_y**2+ df_shear.b_imf_z**2)) * 180 / np.pi

# Compute the magnetosheath beta value
beta_msp_mean = 2 * mu_0 * df_shear.np_msp_mean.values * 1e6 * k_B * (
                2 * df_shear.tp_para_msp_mean.values +
                df_shear.tp_perp_msp_mean.values) / (3 * b_mag_msp ** 2 * 1e-18)

beta_msh_mean = 2 * mu_0 * df_shear.np_msh_mean.values * 1e6 * k_B * (
                2 * df_shear.tp_para_msh_mean.values +
                df_shear.tp_perp_msh_mean.values) / (3 * b_mag_msh ** 2 * 1e-18)

delta_beta = beta_msh_mean - beta_msp_mean

# df_shear["msh_msp_shear"] = msh_msp_shear
# df_shear["delta_beta"] = delta_beta

df_list = [df_shear, df_rx_en, df_va_cs, df_bisec]
for dfn in df_list:
    dfn["msh_msp_shear"] = msh_msp_shear
    dfn["delta_beta"] = delta_beta.values
    # Set values of beta_shear to nan if it is smaller than 0 or larger than 100
    dfn.loc[dfn["delta_beta"] < 0, "delta_beta"] = np.nan
    dfn.loc[dfn["delta_beta"] > 100, "delta_beta"] = np.nan
    # Multiply the delta_beta by 10
    dfn["delta_beta"] = dfn["delta_beta"]
    dfn["cone_angle"] = cone_angle.values



# For each dataframe in the list, divide the temperature by 1e6
for df in df_list:
    df["tp_para_msh_mean"] = df["tp_para_msh_mean"] / 1e6
    df["tp_para_msh_median"] = df["tp_para_msh_median"] / 1e6
    df["tp_para_msp_mean"] = df["tp_para_msp_mean"] / 1e6
    df["tp_para_msp_median"] = df["tp_para_msp_median"] / 1e6
    df["tp_perp_msh_mean"] = df["tp_perp_msh_mean"] / 1e6
    df["tp_perp_msh_median"] = df["tp_perp_msh_median"] / 1e6
    df["tp_perp_msp_mean"] = df["tp_perp_msp_mean"] / 1e6
    df["tp_perp_msp_median"] = df["tp_perp_msp_median"] / 1e6
label = ["Shear", "Rx En", "Va Cs", "Bisec"]
key_list = ["b_imf_z", "b_imf_x", "b_imf_y", "imf_clock_angle", "beta_msh_mean", "np_msp_median",
            "tp_para_msp_median", "tp_perp_msp_median", "msh_msp_shear", "cone_angle", "delta_beta"]
key2_list = ["IMF $B_{\\rm z}$ (nT)", "IMF $B_{\\rm x}$ (nT)", "IMF $B_{\\rm y} (nT)$",
             "IMF Clock Angle (${~}^{0}$)", "$\\beta_{\\rm p}$", "$N_p$ (MSP) (cm$^{-3}$)",
             "$Tp_{\\parallel} (10^6 K)$", "$Tp_{\\perp} (10^6 K)$", "Shear Angle (${~}^{0}$)", "Cone Angle ($B_x/|B|$) (${~}^{0}$)", "$\\Delta \\beta$"]

x_scale_list = [False, False, False, False, False, False, False, False, True, False, False]
y_scale_list = [False, False, False, False, True, True, True, True, False, False, True]

color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
dark_mode = True

if dark_mode:
    plt.style.use('dark_background')
    tick_color = 'w'  # color of the tick lines
    mtick_color = 'w'  # color of the minor tick lines
    label_color = 'w'  # color of the tick labels
    clabel_color = 'w'  # color of the colorbar label
else:
    plt.style.use('default')
    tick_color = 'k'  # color of the tick lines
    mtick_color = 'k'  # color of the minor tick lines
    label_color = 'k'  # color of the tick labels
    clabel_color = 'k'  # color of the colorbar label

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

label_fontsize = 15
tick_fontsize = 12
data_type = ["shear", "rx_en", "va_cs", "bisec"]
ind1 = 0
ind2 = 10
for i, (key, key2) in enumerate(zip(key_list[ind1:ind2], key2_list[ind1:ind2])):
    if key == "msh_msp_shear":
        axs_list = spf.seaborn_subplots(df_list=df_list, keys=["delta_beta", "msh_msp_shear"],
                                        labels=[r"$\Delta \beta$", key2],
                                        data_type=data_type, color_list=color_list, log_scale=False,
                                        x_log_scale=True, y_log_scale=False,
                                        fig_name=None, fig_format="png", nbins=[20, 20],
                                        dark_mode=dark_mode)
    else:
        axs_list = spf.seaborn_subplots(df_list=df_list, keys=["r_rc", key],
                                       labels=[r"Reconnection Distance $(R_\oplus)$", key2],
                                       data_type=data_type, color_list=color_list, log_scale=False,
                                       x_log_scale=x_scale_list[i], y_log_scale=y_scale_list[i],
                                       fig_name=None, fig_format="png", nbins=[20, 20],
                                       dark_mode=dark_mode)

# plt.show()
plt.close('all')
'''

for i, df in enumerate(df_list):
    spf.kde_plots(df=df, x='b_imf_z', y='b_imf_y', log_scale=False, y_log_scale=False,
                  xlim=[-10, 8], ylim=[-10, 15],
                  marker_size=20*df.r_rc.values, alpha=0.7, color=color_list[i],
                  data_type=data_type[i], x_label=r"$B_{\rm z}$ [nT]", y_label=r"$B_{\rm y}$ [nT]")
