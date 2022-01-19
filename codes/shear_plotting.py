#from rx_model_batch import import pylab as pl
import pylab as pl
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import h5py as hf
import datetime

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

# Set the fontstyle to Times New Roman
font = { 'family' : 'sans-serif', 'weight' : 'normal', 'size' : 10 }
plt.rc( 'font', **font )
plt.rc('text', usetex=True)

trange = ['2016-12-24 15:08', '2016-12-24 15:12']

#y, z, rx_en, shear, va_cs, bisec = rx_model_batch(trange=trange)
model_type = 't01'  # t01
dr = 0.5  # Resolution of model run in R_E units
mp = 0.5  # Magnetopause thichkness
fn = f'../data/all_data_rx_model_{dr}re_{mp}mp_{model_type}_2021-10-20.h5'
fn = "../data/total_fields_0.5_15_15_[2.661813864565469, -18.0, -5.79, 0.37, 0, 0, 0, 0, 0, 0]_2021-11-12_v01.h5"
dat = hf.File(fn)

cmap = plt.cm.viridis

pad = 0.02
clabelpad = 10
labelsize = 18
ticklabelsize = 15
cticklabelsize = 15
clabelsize = 15
ticklength = 3
tickwidth = 1.5
ticklength = 6
mticklength = 4
cticklength = 5
mcticklength = 4
labelrotation = 0

try:
    plt.close('all')
except:
    pass
ex = np.linspace(-15.1, 15, int(30/dr) + 1)
ey = np.linspace(-15.1, 15, int(30/dr) + 1)
dat_ipy = dat
for key in list(dat.keys())[:]:
    '''

    # Define the figure
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    min_val = np.nanmin(dat_ipy[key][:,:])
    max_val = np.nanmax(dat_ipy[key][:,:])
    levels=np.linspace(min_val, max_val, 30)
    c1 = axes.contourf(dat_ipy['y_shu'][:,:], dat_ipy['z_shu'][:,:], dat_ipy[key][:,:],
                       levels=levels, cmap=plt.cm.viridis, norm=mpl.colors.Normalize())
    plt.colorbar(c1, ax=axes, ticks=np.linspace(min_val, max_val, 11))
    patch = patches.Circle((0, 0), radius=15, transform=axs1.transData)
    axes.set_clip_path(patch)

    axes.set_xlim(-15, 15)
    axes.set_ylim(-15, 15)
    axes.set_xlabel( r'$\mathrm{y_{GSM}} ( R_E )$', fontsize=20 )
    axes.set_ylabel( r'$\mathrm{z_{GSM}} ( R_E )$', fontsize=20 )
    axes.set_title(key, fontsize=20)
    '''
    fig = plt.figure(num=None, figsize=(6, 6), dpi=200, facecolor='w', edgecolor='gray')
    fig.subplots_adjust(left=0.01, right=0.95, top=0.99, bottom=0.01, wspace=0.02, hspace=0.)

    # Define the axes in the figure
    axs1 = fig.add_subplot(1, 1, 1)
    norm = mpl.colors.Normalize()

    im1 = axs1.pcolormesh(ex, ey, np.transpose(dat[key]))
    patch = patches.Circle((0, 0), radius=15, transform=axs1.transData)
    im1.set_clip_path(patch)

    axs1.set_xlim(-15, 15)
    axs1.set_ylim(-15, 15)

    axs1.set_xlabel( r'$\mathrm{y_{GSM}} ( R_\oplus )$', fontsize=20 )
    axs1.set_ylabel( r'$\mathrm{z_{GSM}} ( R_\oplus )$', fontsize=20 )

    axs1.text(0.98, 0.95, f'{model_type}', horizontalalignment='right',
              verticalalignment='center', transform=axs1.transAxes, fontsize=18)

    axs1.text(0.98, 0.02, f'$R_{{mp}}$ = {mp} $R_\oplus$', horizontalalignment='right',
              verticalalignment='center', transform=axs1.transAxes, fontsize=18)

    axs1.text(0.02, 0.95, f'$\Delta R$ = {dr}$R_\oplus$', horizontalalignment='left',
              verticalalignment='center', transform=axs1.transAxes, fontsize=18)


    # Create a new axis for colorbar
    divider1 = make_axes_locatable(axs1)

    # Define the location of the colorbar, it's size relative to main figure and the padding
    # between the colorbar and the figure, the orientation the colorbar
    cax1 = divider1.append_axes("top", size="5%", pad=0.01 )
    cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)

    cbar1.ax.xaxis.set_label_position('top')

    try :

        new_key = key.replace('_', '-')
        cbar1.set_label(f'{new_key}', fontsize=labelsize, labelpad=clabelpad)
    except :
        pass

    cbar1.ax.tick_params(axis='x', direction='in', labeltop=True, labelbottom=False, color='k',
                        top=True, bottom=True)
    '''
    '''
    fig_name = f'../figures/test_contourf_{new_key}_{dr}re_{mp}mp_{today_date}_{model_type}_v1.png'
    plt.savefig( fig_name, bbox_inches='tight', pad_inches = 0.05, format='png', dpi=300)
    #plt.close()
    #fig.tight_layout()
    plt.show()

    print(f'Figures saved as {fig_name}')
dat.close()
