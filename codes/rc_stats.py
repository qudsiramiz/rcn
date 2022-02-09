from turtle import position
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib.pyplot import MaxNLocator

file_name = "../data/reconnection_line_data.csv"

df = pd.read_csv(file_name)

df_shear = df[df.method_used=="shear"]
df_rx_en = df[df.method_used=="rx_en"]
df_va_cs = df[df.method_used=="va_cs"]
df_bisec = df[df.method_used=="bisection"]


# Set the font size for the axes
label_size = 20  # fontsize for x and y labels
t_label_size = 18  # fontsize for tick label
c_label_size = 18  # fontsize for colorbar label
ct_tick_size = 14  # fontsize for colorbar tick labels
l_label_size = 14  # fontsize for legend label

tick_len = 12  # length of the tick lines
mtick_len = 7  # length of the minor tick lines
tick_width = 1  # tick width in points
mtick_width = 0.7  # minor tick width in points

label_pad = 5  # padding between label and axis

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

nbins = 8

plt.close("all")

fig = plt.figure(num=None, figsize=(8,8), dpi=200, facecolor='k', edgecolor='w')
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0., hspace=0.1)
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1])

# Plot the histogram of the shear data
axs1 = plt.subplot(gs[0, 0])
axs1.hist(df_shear.r_rc, bins=nbins, range=(0, 15),color='#1f77b4', alpha=0.5, label='shear')
axs1.set_xlim(0, 15)
#axs1.set_xlabel(r'$r_{rc}$', fontsize=label_size, color=label_color, labelpad=label_pad)
axs1.set_ylabel('Count', fontsize=label_size, color=label_color, labelpad=label_pad)

# Plot the histogram of the rx_en data
axs2 = plt.subplot(gs[0, 1])
axs2.hist(df_rx_en.r_rc, bins=nbins, range=(0, 15), color='#ff7f0e', alpha=0.5, label='rx-en')
axs2.set_xlim(0, 15)
#axs2.set_xlabel(r'$r_{rc}$', fontsize=label_size, color=label_color, labelpad=label_pad)
axs2.set_ylabel('Count', fontsize=label_size, color=label_color, labelpad=label_pad)
axs2.yaxis.set_label_position("right")

# Plot the histogram of the va_cs data
axs3 = plt.subplot(gs[1, 0])
axs3.hist(df_va_cs.r_rc, bins=nbins, range=(0, 15), color='#2ca02c', alpha=0.5, label='va-cs')
axs3.set_xlim(0, 15)
axs3.set_xlabel(r'$R_{\rm {rc}} (R_\oplus)$', fontsize=label_size, color=label_color,
                labelpad=label_pad)
axs3.set_ylabel('Count', fontsize=label_size, color=label_color, labelpad=label_pad)

# Plot the histogram of the bisection data
axs4 = plt.subplot(gs[1, 1])
axs4.hist(df_bisec.r_rc, bins=nbins, range=(0, 15), color='#d62728', alpha=0.5, label='bisection')
axs4.set_xlim(0, 15)
axs4.set_xlabel(r'$R_{\rm {rc}} (R_\oplus)$', fontsize=label_size, color=label_color,
                labelpad=label_pad)
axs4.set_ylabel('Count', fontsize=label_size, color=label_color, labelpad=label_pad)
axs4.yaxis.set_label_position("right")

# Set the tick parameters
axs1.tick_params(axis='both', direction='in', which='major', left=True, right=True,
                 top=True, bottom=True, labelleft=True, labelright=False,
                 labeltop=False, labelbottom=True, labelsize=t_label_size,
                 length=tick_len, width=tick_width, labelcolor=label_color)

axs2.tick_params(axis='both', direction='in', which='major', left=True, right=True,
            top=True, bottom=True, labelleft=False, labelright=True,
            labeltop=False, labelbottom=True, labelsize=t_label_size,
            length=tick_len, width=tick_width, labelcolor=label_color)

axs3.tick_params(axis='both', direction='in', which='major', left=True, right=True,
                 top=True, bottom=True, labelleft=True, labelright=False,
                 labeltop=False, labelbottom=True, labelsize=t_label_size,
                 length=tick_len, width=tick_width, labelcolor=label_color)

axs4.tick_params(axis='both', direction='in', which='major', left=True, right=True,
            top=True, bottom=True, labelleft=False, labelright=True,
            labeltop=False, labelbottom=True, labelsize=t_label_size,
            length=tick_len, width=tick_width, labelcolor=label_color)

axs1.text(0.95, .95, 'Shear', ha='right', va='top',
            transform=axs1.transAxes, fontsize=c_label_size, color=label_color)
axs2.text(0.95, .95, 'Reconnection\n Energy', ha='right', va='top',
            transform=axs2.transAxes, fontsize=c_label_size, color=label_color)
axs3.text(0.95, .95, 'Exhaust\n Velocity', ha='right', va='top',
            transform=axs3.transAxes, fontsize=c_label_size, color=label_color)
axs4.text(0.95, .95, 'Bisection', ha='right', va='top',
            transform=axs4.transAxes, fontsize=c_label_size, color=label_color)

# Show minor ticks
axs1.minorticks_on()
axs1.tick_params(axis='both', which='minor', direction='in', length=mtick_len, left=True,
                 right=True, top=True, bottom=True, color=mtick_color, width=mtick_width)
axs2.minorticks_on()
axs2.tick_params(axis='both', which='minor', direction='in', length=mtick_len, left=True,
                 right=True, top=True, bottom=True, color=mtick_color, width=mtick_width)
axs3.minorticks_on()
axs3.tick_params(axis='both', which='minor', direction='in', length=mtick_len, left=True,
                 right=True, top=True, bottom=True, color=mtick_color, width=mtick_width)
axs4.minorticks_on()
axs4.tick_params(axis='both', which='minor', direction='in', length=mtick_len, left=True,
                 right=True, top=True, bottom=True, color=mtick_color, width=mtick_width)

# Set the number of ticks on the x-and y-axis
axs1.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
axs1.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
axs2.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
axs2.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
axs3.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
axs3.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
axs4.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
axs4.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))

# Setting the tickmarks labels in such a way that they don't overlap
plt.setp(axs1.get_xticklabels(), rotation=0, ha='right', va='top', visible=True)
plt.setp(axs1.get_yticklabels(), rotation=0, va='center', visible=True)

title_y_pos = 1.03

if dark_mode:
    fig.suptitle(f'Histogram of the distance of reconnection site from spacecraft location',
                 fontsize=label_size, color='w', y=title_y_pos, alpha=0.65)
else:
    fig.suptitle(f'Histogram of the distance of reconnection site from spacecraft location',
                 fontsize=label_size, color='crimson', y=title_y_pos, alpha=1)

if dark_mode:
    transparent=False
else:
    transparent=True

plt.savefig("../figures/reconnection_line_histogram.pdf", bbox_inches='tight', pad_inches=0.05,
            dpi=200, transparent=transparent, format='pdf')