import importlib
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import rc_stats_fncs as rcsf
importlib.reload(rcsf)

# Set usetex to True
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)


data_folder = '../data/rx_d'

fnames = np.sort(glob.glob(f"{data_folder}/reconnection_line_data_mms3_20221109.csv"))

df = pd.read_csv(fnames[0], index_col=False)

# Set date_from as index
df = df.set_index("date_from")

df_shear = df[df.method_used == "shear"]
df_rx_en = df[df.method_used == "rx_en"]
df_va_cs = df[df.method_used == "va_cs"]
df_bisec = df[df.method_used == "bisection"]

color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#606060']

tick_len = 12
tick_wid = 1.5
label_size = 20
label_pad = 5

# Create a figure and axis for two subplots
fig, ax = plt.subplots(1, 2, figsize=(17, 7.5))


# Set wsapce between subplots
fig.subplots_adjust(wspace=0.0, hspace=0.0)

# Make sure that the aspect ratio is
# ax[0].set_aspect('auto')

circle = plt.Circle((0, 0), 1, color='k', alpha=0.5, fill=False)

# Add circle to plot
ax[0].add_artist(circle)

# Create a rage of x and y values within the circle's radius
x = np.linspace(-1, 1, 1000)
y = np.sqrt(1 - x**2)

# Shade the lower h-af the circle
ax[0].fill_between(x, -y, where=y >= 0, color='k', alpha=0.5)


# Get all the reconnection locations where df_shear.b_imf_z > 0
df_shear_pos_z = df_shear[df_shear.b_imf_z > 0]

# Get all the reconnection locations where df_shear.b_imf_z < 0
df_shear_neg_z = df_shear[df_shear.b_imf_z < 0]

# Plot the reconnection locations for shear method
marker_size = 10 * df_shear_pos_z.r_rc.values
ax[0].scatter(df_shear_pos_z.spc_pos_y, df_shear_pos_z.spc_pos_x, marker='o', color=color_list[1],
              alpha=0.5, s=marker_size, label=r'$B_{\rm Z} > 0$', facecolors=color_list[1],
              edgecolors=color_list[1])

marker_size = 10 * df_shear_neg_z.r_rc.values

ax[0].scatter(df_shear_neg_z.spc_pos_y, df_shear_neg_z.spc_pos_x, marker='o', color=color_list[2],
              s=marker_size, label=r'$B_{\rm Z} < 0$', facecolors=color_list[2],
              edgecolors=color_list[2], alpha=0.5)
ax[0].set_xlabel(r'Y [GSM, $R_{\rm E}$]', fontsize=label_size * 1.2)
ax[0].set_ylabel(r'X [GSM, $R_{\rm E}$]', fontsize=label_size * 1.2)

# Change tickmarks properties
ax[0].tick_params(axis='both', which='both', direction='in', labelsize=label_size, top=True,
                  labeltop=False, right=True, labelright=False, bottom=True, labelbottom=True,
                  left=True, labelleft=True, length=tick_len, width=tick_wid, pad=label_pad)

MPangles = np.arange(0, 100 + 1, 1)

Bz = -5  # (nt)
Pdyn = 3  # (nPa)

R0 = 10.22+1.29*(np.tanh(0.184*(Bz + 8.14)))*(Pdyn**(-1/6.6))

alpha = ((0.58-0.007*Bz)*(1+(0.024*np.log(Pdyn))))
Tplus = MPangles*np.pi/180

MPx = (R0 * (2 / (1 + np.cos(Tplus))) ** alpha) * np.cos(Tplus)
MPy = (R0 * (2 / (1 + np.cos(Tplus))) ** alpha) * np.sin(Tplus)
ax[0].plot(MPy, MPx, linestyle='--', color=color_list[0])
ax[0].plot(-MPy, MPx, linestyle='--', color=color_list[0])

ax[0].legend(loc='upper right', fontsize=label_size, frameon=True)

# Add a horozontal arrow at (10, 0) pointing at right at the parabola and add text at the start of
# the arrow pointing at the parabola with the text "parabola"
ax[0].annotate("Magnetopause", xy=(-16.4, 1), xytext=(-5, -3.5), arrowprops=dict(arrowstyle="->",
               color=color_list[3], lw=2), fontsize=label_size, color=color_list[3])

# Write dawn and dusk labels
ax[0].text(13, -4, 'Dusk', fontsize=label_size, color=color_list[-1], ha='left', va='center')
ax[0].text(-13, -4, 'Dawn', fontsize=label_size, color=color_list[0], ha='right', va='center')

# Write panel (a) at the top right corner of the plot
ax[0].text(0.1, 0.95, '(a)', fontsize=1.2 * label_size, color='k', ha='right', va='top',
           transform=ax[0].transAxes)

# Flip the x axis
ax[0].invert_xaxis()

# Plot the reconnection locations for shear method in the yz plane

ax[1].set_aspect('equal')
ax[0].set_aspect('equal')
circle2 = plt.Circle((0, 0), 1, color='k', alpha=0.5, fill=False)

# Add circle to plot
ax[1].add_artist(circle2)

# Set x and y limits to create a square plot

# Create a rage of x and y values within the circle's radius
x = np.linspace(-1, 1, 1000)
y = np.sqrt(1 - x**2)

# Shade the lower h-af the circle
ax[1].fill_between(x, y, -y, where=y >= 0, color='k', alpha=0.5)


marker_size = 10 * df_shear_pos_z.r_rc.values
ax[1].scatter(df_shear_pos_z.spc_pos_y, df_shear_pos_z.spc_pos_z, marker='o', color=color_list[1],
              alpha=0.5, s=marker_size, label=r'$B_{\rm Z} > 0$', facecolors=color_list[1],
              edgecolors=color_list[1])

marker_size = 10 * df_shear_neg_z.r_rc.values
ax[1].scatter(df_shear_neg_z.spc_pos_y, df_shear_neg_z.spc_pos_z, marker='o', color=color_list[2],
              s=marker_size, label=r'$B_{\rm Z} < 0$', facecolors=color_list[2],
              edgecolors=color_list[2], alpha=0.5)

ax[1].set_xlabel(r'Y [GSM, $R_{\rm E}$]', fontsize=label_size * 1.2)
ax[1].set_ylabel(r'Z [GSM, $R_{\rm E}$]', fontsize=label_size * 1.2)

# Add y-label on the right side of the plot
ax[1].yaxis.set_label_position("right")


ax[1].tick_params(axis='both', which='both', direction='in', labelsize=label_size, top=True,
                  labeltop=False, right=True, labelright=True, bottom=True, labelbottom=True,
                  left=True, labelleft=False, length=tick_len, width=tick_wid, pad=label_pad)

ax[1].legend(loc='upper right', fontsize=label_size, frameon=True)


# Write dawn and dusk labels
ax[1].text(7, -9, 'Dusk', fontsize=label_size, color=color_list[-1], ha='left', va='center')
ax[1].text(-7, -9, 'Dawn', fontsize=label_size, color=color_list[0], ha='right', va='center')

# Write panel (b) at the top left corner of the plot
ax[1].text(0.1, 0.95, '(b)', fontsize=1.2 * label_size, color='k', ha='right', va='top',
           transform=ax[1].transAxes)

# Invert the x axis
ax[1].invert_xaxis()

ax[0].set_xlim(-18, 18)
ax[0].set_ylim(-5, 18)

ax[1].set_xlim(-12, 12)
ax[1].set_ylim(-11, 11)


plt.savefig('../figures/mms_rxn_location.pdf', dpi=300, bbox_inches='tight', transparent=True,
            pad_inches=0.05)

plt.close()
