
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
plt.figure(figsize=(8,8))

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
df_list = [df_shear, df_rx_en, df_va_cs, df_bisec]
label = ["Shear", "Rx En", "Va Cs", "Bisec"]
key = "b_imf_y"
key2 = "IMF By"

plt.figure(figsize=(8,8))
for i, df in enumerate(df_list):
    plt.subplot(2, 2, i+1)
    plt.plot(df[key], df.r_rc, c='aqua', marker=".", ls=None, lw=0, label=label[i])
    plt.legend()
    plt.xlabel(f"{key2}")
    plt.ylabel("Reconnection Distance (R_E)")
    # plt.xscale("log")

plt.suptitle(f"Reconnection Distance vs {key2}")
plt.tight_layout()

#plt.plot(df_shear.imf_clock_angle, df_shear.r_rc, c='b', marker=".", ls=None, lw=0, label="Shear")
#plt.plot(df_rx_en.imf_clock_angle, df_rx_en.r_rc, c='orange', marker=".", ls=None, lw=0)
#plt.plot(df_va_cs.imf_clock_angle, df_va_cs.r_rc, c='g', marker=".", ls=None, lw=0)
#plt.plot(df_bisec_imfz.imf_clock_angle, df_bisec_imfz.r_rc, c='r', marker=".", ls=None, lw=0, ms=5, alpha=0.5)
#plt.plot(df_bisec.beta_msh_mean, df_bisec.r_rc, c='w', marker=".", ls=None, lw=0, ms=5, alpha=0.5)
#plt.xscale('log')
plt.savefig(f"../figures/rc_v_{key}.png")
