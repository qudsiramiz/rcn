from mayavi import mlab
import numpy as np
import h5py as hf
import time as tm
start = tm.time()

mlab.options.offscreen = False

dat = hf.File('../data/all_data_rx_model_0.25re_t01_20211013_v01.h5')
#dat = hf.File(fn)

xx = []
count = 0
j = 0
num = 121
for i in range(num*(num+1)):
    if(count >= num):
        count = 0
        j += 1
    else :
        xx.append(j)
        count += 1

yy = np.zeros_like(xx)

for i in range(num):
    for j in range(num):
        yy[num*i+j] = j

bx_0 = np.array(dat['b_msx'][:, :]).flatten()
by_0 = np.array(dat['b_msy'][:, :]).flatten()
bz_0 = np.array(dat['b_msz'][:, :]).flatten()

xx = np.array(xx)
zz = np.zeros_like(xx)
zz[:] = 0

# mayavi.mlab.close(scene=None, all=True)
try:
#mlab.clf()
    mlab.close()
except:
    pass

#q = ax.quiver(xx, yy, zz, bx_24_0, by_24_0, bz_24_0, length=0.01, cmap='gray', lw=3)
#q.set_array(np.random.rand(np.prod(xx.shape)))

#plt.show()
#‘2darrow’ or ‘2dcircle’ or ‘2dcross’ or ‘2ddash’ or ‘2ddiamond’ or 
#‘2dhooked_arrow’ or ‘2dsquare’ or ‘2dthick_arrow’ or ‘2dthick_cross’ 
#or ‘2dtriangle’ or ‘2dvertex’ or ‘arrow’ or ‘axes’ or ‘cone’ or ‘cube’ 
#or ‘cylinder’ or ‘point’ or ‘sphere’. Default: 2darrow

mode = 'arrow'
n_mask = 10
scale_mode = 'vector'
colormap1 = 'bone'
colormap2 = 'Reds'

figure = mlab.figure(fgcolor=(0., 0., 0.), bgcolor=(0.5, 0.5, 0.5))

ax1 = mlab.quiver3d(xx, yy, zz, bx_0, by_0, bz_0, vmin=0, vmax=1, extent=[-15, 15, -15, 15, 0, 0],
                    scale_mode=scale_mode, colormap=colormap1, mode=mode, mask_points=n_mask,
                    figure=figure)

ax1.module_manager.vector_lut_manager.reverse_lut = False

#lut = ax1.module_manager.scalar_lut_manager.lut.table.to_array()
#ilut = lut[::-1]
#ax1.module_manager.scalar_lut_manager.lut.table = ilut
#mlab.draw()
#ax2 = mlab.quiver3d(xx, yy, zz, bx_r_0, by_r_0, bz_r_0, vmin=0, vmax=1, extent=[8, 36, 8, 36, 0, 0],
#scale_mode=scale_mode, colormap=colormap2, mode=mode, mask_points=n_mask, figure=figure)
#ax2.module_manager.vector_lut_manager.reverse_lut = False

#mlab.draw()
#mlab.colorbar(orientation='vertical', nb_labels=5)

mlab.outline()
mlab.axes()

mlab.xlabel('Y')
mlab.ylabel('Z')
# mlab.xlim(-15, 15)
# mlab.ylim(-15, 15)
mlab.view(azimuth=0, elevation=0, distance=75)
    
#mlab.savefig(f'../figures/mayavi_figures/{n_spc}spc/comparison_{n_spc}_spc_{str(z_ind).zfill(3)}.png', size=(2000, 2000))
#mlab.clf()
#mlab.colorbar(orientation='horizontal', nb_labels=5, title='simulation')

# linetype = line, ribbon, tube
# integration_direction = forward, backward, both
# seedtype = line, plane, point, sphere
    
#mayavi.mlab.flow(xx, yy, zz, bx_r_0, by_r_0, bz_r_0, vmin=-1, vmax=1, extent=[8, 36, 8, 36, 0, 9],
#                 linetype='line', integration_direction='forward', seed_visible=True,
#                 seedtype='point', colormap='seismic')
    
#mayavi.mlab.savefig(f'comparison_{n_spc}_spc.png', magnification='auto')

#plt.colorbar()

#print(f'figure saved for idx = {z_ind} for {n_spc} spacecraft')

#mlab.close()
mlab.show()

print(f'It took {round(tm.time() - start, 2)} seconds')
#from numpy import pi, sin, cos, mgrid
#dphi, dtheta = pi/250.0, pi/250.0
#[phi,theta] = mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
#m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
#r = sin(m0*phi)**m1 + cos(m2*phi)**m3 + sin(m4*theta)**m5 + cos(m6*theta)**m7
#x = r*sin(phi)*cos(theta)
#y = r*cos(phi)
#z = r*sin(phi)*sin(theta)
#
## View it.
#from mayavi import mlab
#s = mlab.mesh(x, y, z)
#mlab.show()
#
