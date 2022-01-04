import h5py as hf
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches


# Set the fontstyle to Times New Roman
font = { 'family' : 'Times New Roman', 'weight' : 'normal', 'size' : 10 }
plt.rc( 'font', **font )
plt.rc('text', usetex=True)

dat_idl = hf.File('param_total_edited_shear_v03.h5')

dat_ipy = hf.File('../data/all_data_rx_model_0.5re_20211005_v03.h5')
'''
dat_idl = hf.File('params_total_v01.h5')

bx_igrf_idl = dat_idl['bx_igrf'][:,:]
bx_igrf_idl[bx_igrf_idl<-1000] = np.nan

by_igrf_idl = dat_idl['by_igrf'][:,:]
by_igrf_idl[by_igrf_idl<-1000] = np.nan

bz_igrf_idl = dat_idl['bz_igrf'][:,:]
bz_igrf_idl[bz_igrf_idl<-1000] = np.nan

bx_t96_idl = dat_idl['bx_t96'][:,:]
bx_t96_idl[bx_t96_idl<-1000] = np.nan

by_t96_idl = dat_idl['by_t96'][:,:]
by_t96_idl[by_t96_idl<-1000] = np.nan

bz_t96_idl = dat_idl['bz_t96'][:,:]
bz_t96_idl[bz_t96_idl<-1000] = np.nan

bx_idl = dat_idl['bx'][:,:]
bx_idl[bx_idl<-1000] = np.nan

by_idl = dat_idl['by'][:,:]
by_idl[by_idl<-1000] = np.nan

bz_idl = dat_idl['bz'][:,:]
bz_idl[bz_idl<-1000] = np.nan

ys_idl = dat_idl['y_shu'][:,:]
ys_idl[ys_idl<-1000] = np.nan

zs_idl = dat_idl['z_shu'][:,:]
zs_idl[zs_idl<-1000] = np.nan

sr_idl = dat_idl['shear'][:,:]
sr_idl[sr_idl<-1000] = np.nan
sr_idl[sr_idl>200] = np.nan


fn = 'param_total_edited_shear_v3.h5'
new_dat = hf.File(fn, 'w')
new_dat.create_dataset('bx', data=np.transpose(bx_idl))
new_dat.create_dataset('by', data=np.transpose(by_idl))
new_dat.create_dataset('bz', data=np.transpose(bz_idl))

new_dat.create_dataset('bx_igrf', data=np.transpose(bx_igrf_idl))
new_dat.create_dataset('by_igrf', data=np.transpose(by_igrf_idl))
new_dat.create_dataset('bz_igrf', data=np.transpose(bz_igrf_idl))

new_dat.create_dataset('bx_t96', data=np.transpose(bx_t96_idl))
new_dat.create_dataset('by_t96', data=np.transpose(by_t96_idl))
new_dat.create_dataset('bz_t96', data=np.transpose(bz_t96_idl))

new_dat.create_dataset('y_shu', data=np.transpose(ys_idl))
new_dat.create_dataset('z_shu', data=np.transpose(zs_idl))
new_dat.create_dataset('shear', data=np.transpose(sr_idl))

new_dat.close()
'''
bx_idl_arr = dat_idl['bx'][:,:].flatten()
by_idl_arr = dat_idl['by'][:,:].flatten()
bz_idl_arr = dat_idl['bz'][:,:].flatten()


bx_igrf_idl_arr = dat_idl['bx_igrf'][:,:].flatten()
by_igrf_idl_arr = dat_idl['by_igrf'][:,:].flatten()
bz_igrf_idl_arr = dat_idl['bz_igrf'][:,:].flatten()

bx_t96_idl_arr = dat_idl['bx_t96'][:,:].flatten()
by_t96_idl_arr = dat_idl['by_t96'][:,:].flatten()
bz_t96_idl_arr = dat_idl['bz_t96'][:,:].flatten()

ys_idl_arr = dat_idl['y_shu'][:,:].flatten()
zs_idl_arr = dat_idl['z_shu'][:,:].flatten()
sr_idl_arr = dat_idl['shear'][:,:].flatten()

bx_ipy_arr = dat_ipy['bx'][:,:].flatten()
by_ipy_arr = dat_ipy['by'][:,:].flatten()
bz_ipy_arr = dat_ipy['bz'][:,:].flatten()

bx_igrf_ipy_arr = dat_ipy['bx_igrf'][:,:].flatten()
by_igrf_ipy_arr = dat_ipy['by_igrf'][:,:].flatten()
bz_igrf_ipy_arr = dat_ipy['bz_igrf'][:,:].flatten()

bx_t96_ipy_arr = dat_ipy['bx_t96'][:,:].flatten()
by_t96_ipy_arr = dat_ipy['by_t96'][:,:].flatten()
bz_t96_ipy_arr = dat_ipy['bz_t96'][:,:].flatten()

ys_ipy_arr = dat_ipy['y_shu'][:,:].flatten()
zs_ipy_arr = dat_ipy['z_shu'][:,:].flatten()
sr_ipy_arr = dat_ipy['shear'][:,:].flatten()

diff_bx = bx_ipy_arr - bx_idl_arr
diff_by = by_ipy_arr - by_idl_arr
diff_bz = bz_ipy_arr - bz_idl_arr

diff_bx_igrf = bx_igrf_ipy_arr - bx_igrf_idl_arr
diff_by_igrf = by_igrf_ipy_arr - by_igrf_idl_arr
diff_bz_igrf = bz_igrf_ipy_arr - bz_igrf_idl_arr

diff_bx_t96 = bx_t96_ipy_arr - bx_t96_idl_arr
diff_by_t96 = by_t96_ipy_arr - by_t96_idl_arr
diff_bz_t96 = bz_t96_ipy_arr - bz_t96_idl_arr

diff_ys = ys_ipy_arr - ys_idl_arr
diff_zs = zs_ipy_arr - zs_idl_arr

diff_sr = sr_ipy_arr - sr_idl_arr


'''
plt.figure()

plt.plot(bz_idl_arr, bz_ipy_arr, 'r.', ms=2)
plt.xlabel('IDL')
plt.ylabel('python')

plt.savefig('../figures/data_comp_zs.png')

plt.close()
'''

key_list = ['bx', 'by', 'bz', 'y_shu', 'z_shu', 'shear']

data_type = 'idl'
if (data_type=='idl'):
    dat = dat_idl
else:
    dat = dat_ipy

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

ex = np.linspace(-40, 40, 150)
ey = np.linspace(-40, 40, 150)
'''
for key in key_list :
    # Define the figure
    fig = plt.figure(num=None, figsize=(6, 6), dpi=200, facecolor='w', edgecolor='gray')
    fig.subplots_adjust(left=0.01, right=0.95, top=0.99, bottom=0.01, wspace=0.02, hspace=0.)

    # Define the axes in the figure
    axs1 = fig.add_subplot(1, 1, 1)
    norm = mpl.colors.Normalize(vmin=np.nanmin(dat[key]), vmax=np.nanmax(dat[key]))

    im1 = axs1.pcolormesh(ex, ey, np.transpose(dat[key]), alpha=0.9, shading='auto', cmap=cmap,
                          norm=norm)
    patch = patches.Circle((0, 0), radius=15, transform=axs1.transData)
    im1.set_clip_path(patch)

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

    try :
        cbar1.set_label(f'{key}', fontsize=labelsize, labelpad=clabelpad)
    except :
        pass

    cbar1.ax.tick_params(axis='x', direction='in', labeltop=True, labelbottom=False, color='k',
                        top=True, bottom=True)

    fig_name = f'../figures/{key}_0.5re_20211004_v3.png'
    plt.savefig( fig_name, bbox_inches='tight', pad_inches = 0.05, format='png', dpi=300)
    plt.close()

    print(f'Figures saved as {fig_name}')
dat.close()
plt.close('all')

'''
# Define the figure
fig = plt.figure(num=None, figsize=(15, 8), dpi=200, facecolor='w', edgecolor='gray')
fig.subplots_adjust(left=0.01, right=0.95, top=0.99, bottom=0.01, wspace=0.02, hspace=0.)

# Define the axes in the figure
axs1 = fig.add_subplot(1, 2, 1)

axs1.plot(ys_idl_arr, diff_bx_t96, 'r.', ms=3, label=r'$\delta b_x$')
axs1.plot(ys_idl_arr, diff_by_t96, 'b.', ms=3, label=r'$\delta b_y$')
axs1.plot(ys_idl_arr, diff_bz_t96, 'c.', ms=3, label=r'$\delta b_z$')
axs1.legend(fontsize=15)
axs1.set_xlabel(r'$Y_{shu}$', fontsize=18)
axs1.set_ylabel(r'$B_{ipy} - B_{idl} $', fontsize=18)


axs2 = fig.add_subplot(1, 2, 2)

axs2.plot(bx_t96_idl_arr, bx_t96_ipy_arr, 'r.', ms=3, label=r'$b_x$')
axs2.plot(by_t96_idl_arr, by_t96_ipy_arr, 'b.', ms=3, label=r'$b_y$')
axs2.plot(bz_t96_idl_arr, bz_t96_ipy_arr, 'c.', ms=3, label=r'$b_z$')
axs2.legend(fontsize=15)
axs2.set_xlabel(r'$B_{idl}$', fontsize=18)
axs2.set_ylabel(r'$B_{ipy}$', fontsize=18)
axs2.tick_params(axis='y', direction='out', left=False, labelleft=False, right=True,
                    labelright=True)

axs2.set_xlim(-50, 50)
axs2.set_ylim(-50, 50)
axs2.set_title('T96 field')
fig_name = f'../figures/data_comparison_t96_v03.pdf'

plt.savefig( fig_name, bbox_inches='tight', pad_inches = 0.05, format='pdf', dpi=300)
plt.close()

'''
import pylab as pl
fig, axes = pl.subplots(1, 1, figsize=(5, 5))
c1 = axes.contourf(dat_ipy['y_shu'][:,:], dat_ipy['z_shu'][:,:], dat_ipy['bx'][:,:], extent=[-40, 40, -40, 40], levels=30, vmin=-30, vmax=30)
pl.colorbar(c1, ax=axes, cmin=-30, cmax=30)
axes.set_xlim(-15, 15)
axes.set_ylim(-15, 15);

'''
