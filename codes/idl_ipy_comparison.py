import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as hf
import datetime

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

# Set the fontstyle to Times New Roman
font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

# Read both the datafiles
dat_idl = hf.File('idl_data_edited.h5')

dat_ipy = hf.File('python_data.h5')

xx = np.linspace(-15.1, 15, 61)
x_ind = 10  #  10, 29, 32, 50
x_val = round(xx[x_ind], 2)

cmap = plt.cm.bwr

# Define the plot parameters
pad = 0.02
clabelpad = 10
labelsize = 18
ticklabelsize = 18
cticklabelsize = 18
clabelsize = 18
ticklength = 3
tickwidth = 2.5
ticklength = 6
mticklength = 4
cticklength = 5
mcticklength = 4
labelrotation = 0

# Define the edges of the pixels
ex = np.linspace(-15.1, 15, 61)
ey = np.linspace(-15.1, 15, 61)

key_type = input('Type of field to be plotted (dipole(1), external(2), total(3)) ==> ')
if key_type == '1':
    key_list = ['bx_igrf', 'by_igrf', 'bz_igrf']
elif key_type == '2':
    key_list = ['bx_t96', 'bx_t01', 'by_t96', 'by_t01', 'bz_t96', 'bz_t01']
elif key_type =='3':
    key_list = ['bx_t96_total', 'bx_t01_total', 'by_t96_total', 'by_t01_total', 'bz_t96_total',
                'bz_t01_total']
else:
    raise ValueError('Key type value must be 1 (for dipole field), 2 (for external field) or 3 (\
                     for total field)')

# Plot graphs for each key in the datafile
for key in key_list[0:]:

    fig = plt.figure(num=None, figsize=(6, 6), dpi=200, facecolor='w', edgecolor='gray')
    fig.subplots_adjust(left=0.01, right=0.95, top=0.99, bottom=0.01, wspace=0.08, hspace=0.)

    diff = dat_idl[key][x_ind, :, :] - dat_ipy[key][x_ind, :, :]

    key = key.replace('_', '-')

    # Compute the min/max value of colorbar
    max_val = np.nanmax([abs(np.nanmin(diff)), abs(np.nanmax(diff))])

    axs1 = fig.add_subplot(1, 1, 1)
    norm = mpl.colors.Normalize(vmin=-max_val, vmax=max_val)

    im1 = axs1.pcolormesh(ex, ey, np.transpose(diff), alpha=0.9, shading='auto', cmap=cmap,
                          norm=norm)
    axs1.text(0.85, 0.95, f'dst = {dst}', horizontalalignment='center', verticalalignment='center',
              transform=axs1.transAxes, fontsize=18)

    axs1.text(0.12, 0.95, f'X = {x_val}$R_E$', horizontalalignment='center',
              verticalalignment='center', transform=axs1.transAxes, fontsize=18)

    axs1.tick_params(axis='both', direction='in', top=True, labeltop=False, bottom=True, 
                     labelbottom=True, left=True, labelleft=True, right=True, labelright=False,
                     color='k', labelsize=16)

    axs1.set_xlim(-15, 15)
    axs1.set_ylim(-15, 15)

    axs1.set_xlabel( r'$\mathrm{y_{GSM}} ( R_E )$', fontsize=20 )
    axs1.set_ylabel( r'$\mathrm{z_{GSM}} ( R_E )$', fontsize=20 )
    # Create a new axis for colorbar
    divider1 = make_axes_locatable(axs1)

    # Define the location of the colorbar, it's size relative to main figure and the padding
    # between the colorbar and the figure, the orientation the colorbar
    cax1 = divider1.append_axes("top", size="5%", pad=0.01 )
    cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)

    cbar1.ax.xaxis.set_label_position('top')

    cbar1.ax.tick_params(axis='x', direction='in', labeltop=True, labelbottom=False, color='k',
                        top=True, bottom=True, labelsize=15)
    cax1.text(0.5, 0.5, f'{key}', horizontalalignment='center', verticalalignment='center',
              transform=cax1.transAxes, fontsize=16)

    fig_name = f'../figures/{key}_comparison_{today_date}_{x_ind}xind_v01.png'
    plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format='png', dpi=300)
    #plt.show()
    plt.close()

    print(f'Figures saved as {fig_name}')
