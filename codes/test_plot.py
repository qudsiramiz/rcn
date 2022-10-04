import importlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import rx_model_funcs as rxmf
import seaborn_plots as sp
import matplotlib.gridspec as gridspec
import seaborn as sns; sns.set()
import SeabornFig2Grid as sfg

importlib.reload(rxmf)
importlib.reload(sp)

file_name = "../data/rx_d/reconnection_line_data_mms3_20220927.csv"
cut_type_list = ["jet", "walen1", "walen2", "walen_jet"]

df = pd.read_csv(file_name, index_col=False)

df_shear = df[df.method_used=="shear"]
df_rx_en = df[df.method_used=="rx_en"]
df_va_cs = df[df.method_used=="va_cs"]
df_bisec = df[df.method_used=="bisection"]

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
for xx,yy in zip(b_lmn_vec_msh, b_lmn_vec_msp):
    try:
        msh_msp_shear.append(rxmf.get_shear(xx,yy, angle_unit="degrees"))
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

# Compute the magnetosheath beta value
beta_msp_mean = 2 * mu_0 * df_shear.np_msp_mean.values * 1e6 * k_B * (2 *
                df_shear.tp_para_msp_mean.values +
                df_shear.tp_perp_msp_mean.values) / (3 * b_mag_msp ** 2 * 1e-18)

beta_msh_mean = 2 * mu_0 * df_shear.np_msh_mean.values * 1e6 * k_B * (2 *
                df_shear.tp_para_msh_mean.values +
                df_shear.tp_perp_msh_mean.values) / (3 * b_mag_msh ** 2 * 1e-18)

delta_beta = beta_msh_mean - beta_msp_mean

df_shear["msh_msp_shear"] = msh_msp_shear
df_shear["delta_beta"] = delta_beta

df_list = [df_shear, df_rx_en, df_va_cs, df_bisec]
label = ["Shear", "Rx En", "Va Cs", "Bisec"]
key_list = ["b_imf_z", "b_imf_x", "b_imf_y", "imf_clock_angle", "beta_msh_mean", "np_msp_median",
            "tp_para_msp_median", "tp_perp_msp_median"]
key2_list = ["IMF $B_{\\rm z}$ (nT)", "IMF $B_{\\rm x}$ (nT)", "IMF $B_{\\rm y} (nT)$",
             "IMF Clock Angle (${~}^{0}$)", "$\\beta_{\\rm p}$", "$N_p$ (MSP) (cm$^{-3}$)",
             "$Tp_{\\parallel} (K)$", "$Tp_{\\perp} (K)$"]

color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
dark_mode = True

if dark_mode:
    plt.style.use('dark_background')
    tick_color = 'w' # color of the tick lines
    mtick_color = 'w' # color of the minor tick lines
    label_color = 'w' # color of the tick labels
    clabel_color = 'w' # color of the colorbar label
else:
    plt.style.use('default')
    tick_color = 'k' # color of the tick lines
    mtick_color = 'k' # color of the minor tick lines
    label_color = 'k' # color of the tick labels
    clabel_color = 'k' # color of the colorbar label

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

label_fontsize = 15
tick_fontsize = 12

for key, key2 in zip(key_list, key2_list):
    plt.figure(figsize=(8,8))
    for i, df in enumerate(df_list):
        # Find the spearman and pearson correlation between key and "r_rc"
        spearman = df[key].corr(df["r_rc"], method="spearman")
        pearson = df[key].corr(df["r_rc"], method="pearson")

        plt.subplot(2, 2, i+1)
        plt.plot(df.r_rc, df[key], c=color_list[i], marker=".", ls=None, lw=0, label=label[i])
        # Make a line with slope of spearman correlation coefficient
        x = np.linspace(0, 25, 100)
        y = spearman*x + np.mean(df[key]) - spearman*np.mean(df.r_rc)
        plt.plot(x, y, c="w", ls="--", lw=2)

        # Get the r^2 value
        r2 = pearson**2

        # Print the correlation coefficient on the plot in a box
        plt.text(0.02, 0.02, f"$\\rho_{{\\rm {'s'}}}$ = {spearman:.2f}\n"
                             f"$\\rho_{{\\rm {'p'}}}$ = {pearson:.2f}",
                 transform=plt.gca().transAxes, va="bottom", ha="left",
                 bbox=dict(facecolor='k', alpha=1, edgecolor='k', boxstyle='round,pad=0.2'))
        plt.legend(loc=4, frameon=False, fontsize=10, ncol=1, handlelength=0.1)
        plt.ylabel(f"{key2}")
        plt.xlabel("Reconnection Distance $(R_\oplus)$")
        if key == "beta_msh_mean" or key=="tp_para_msp_median" or key=="tp_perp_msp_median":
            plt.yscale("log")
        # if key == "b_imf_z":
        #     plt.xlim(-5, 5)
        plt.xlim(0, 20)
        # Reverse the x-axis
        if key == "imf_clock_angle":
            plt.gca().invert_yaxis()
    plt.suptitle(f"Reconnection Distance vs {key2}")
    plt.tight_layout()

#plt.plot(df_shear.imf_clock_angle, df_shear.r_rc, c='b', marker=".", ls=None, lw=0, label="Shear")
#plt.plot(df_rx_en.imf_clock_angle, df_rx_en.r_rc, c='orange', marker=".", ls=None, lw=0)
#plt.plot(df_va_cs.imf_clock_angle, df_va_cs.r_rc, c='g', marker=".", ls=None, lw=0)
#plt.plot(df_bisec_imfz.imf_clock_angle, df_bisec_imfz.r_rc, c='r', marker=".", ls=None, lw=0, ms=5, alpha=0.5)
#plt.plot(df_bisec.beta_msh_mean, df_bisec.r_rc, c='w', marker=".", ls=None, lw=0, ms=5, alpha=0.5)
#plt.xscale('log')
    plt.savefig(f"../figures/rc_v_{key}.png")

shear_angle_theory = np.logspace(-1, np.log10(180), 100)
delta_beta_theory_half = np.tan(np.deg2rad(shear_angle_theory/2))
delta_beta_theory_one = 2 * np.tan(np.deg2rad(shear_angle_theory/2))
delta_beta_theory_two = 4 * np.tan(np.deg2rad(shear_angle_theory/2))

# Plot the delta beta vs shear angle
plt.figure(figsize=(8,8))
plt.scatter(df_shear.delta_beta, df_shear.msh_msp_shear, c=color_list[1], marker=".",
            s=25*df_shear.r_rc.values, alpha=0.7)
plt.plot(delta_beta_theory_half, shear_angle_theory, marker='.', c="w", ls=None, lw=0, label="Half")
plt.plot(delta_beta_theory_one, shear_angle_theory, marker='.', c="b", ls=None, lw=0, label="One")
plt.plot(delta_beta_theory_two, shear_angle_theory, marker='.', c="g", ls=None, lw=0, label="Two")

lgnd = plt.legend(loc=4, frameon=False, fontsize=10, ncol=1, handlelength=0.1)
for handle in lgnd.legendHandles:
    handle.size = [1]

plt.xlabel("$\\Delta \\beta$", fontsize=label_fontsize)
plt.ylabel("Shear Angle (${~}^{0}$)", fontsize=label_fontsize)
plt.title("Shear Angle vs $\\Delta \\beta$", fontsize=1.2 * label_fontsize)
plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
plt.xlim(1e-3, 2e1)
plt.ylim(1e-1, 2e2)
plt.xscale('log')
plt.yscale('log')
plt.savefig("../figures/delta_beta_v_shear_angle_mean.png")

plt.close('all')

# sp.kde_plots(df_shear.delta_beta.values, df_shear.msh_msp_shear.values)


axs1 = sp.kde_plots(df=df_shear, x="delta_beta", y="msh_msp_shear",
                    x_label=r"$\Delta \beta_{\rm msh - msp}$", y_label=r"Shear Angle (${~}^{0}$)",
                    log_scale=True)

axs2 = sp.kde_plots(df=df_shear, x="delta_beta", y="msh_msp_shear",
                    x_label=r"$\Delta \beta_{\rm msh - msp}$", y_label=r"Shear Angle (${~}^{0}$)",
                    log_scale=True)


fig = plt.figure(figsize=(40,80))
gs = gridspec.GridSpec(2, 1)

mg0 = sfg.SeabornFig2Grid(axs1, fig, gs[0])
mg1 = sfg.SeabornFig2Grid(axs2, fig, gs[1])
#mg2 = sfg.SeabornFig2Grid(g2, fig, gs[3])
#mg3 = sfg.SeabornFig2Grid(g3, fig, gs[2])

gs.tight_layout(fig)
#gs.update(top=0.7)

plt.savefig("../figures/test.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close('all')