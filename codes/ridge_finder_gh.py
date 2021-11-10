import time
import numpy as np
import h5py as hf
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import signal, interpolate
from astropy.convolution import convolve
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Set the fontstyle to Times New Roman
font = { 'family' : 'serif', 'weight' : 'normal', 'size' : 10 }
plt.rc( 'font', **font )
plt.rc('text', usetex=True)

start = time.time()

# Load data
dat = hf.File('../data/all_data_rx_model_0.25re_0.5mp_t96_2021-10-20.h5')

image = dat['shear'][:, :]

y_shu = dat['y_shu'][:, :]
z_shu = dat['z_shu'][:, :]

y_shu_min = np.nanmin(y_shu)
z_shu_min = np.nanmin(z_shu)

y_shu_max = np.nanmax(y_shu)
z_shu_max = np.nanmax(z_shu)

y_arr = np.linspace(y_shu_min, y_shu_max, np.shape(y_shu)[0])
z_arr = np.linspace(z_shu_min, z_shu_max, np.shape(z_shu)[0])

print('Creating interpolation function...\n')

# Create an interpolation function for the dataset
interp_func = interpolate.RectBivariateSpline(y_arr, z_arr, image)

print('Interpolation done \n')

# interp = input('Carry out interpolation process? (y/n)==>  ')
interp=''
print('')

if(interp):
    print('Creating interpolation function...\n')
    interp_func = interpolate.interp2d(y_shu, z_shu, image)
    print('Interpolation done \n')

# Define a kernel for smoothing the image
kernel = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]])
#kernel = np.array([[1, 3, 3, 1], [3, 9, 9, 3], [3, 9, 9, 3], [1, 3, 3, 1]])
#kernel = np.array([[1, 3, 5, 3, 1], [1, 3, 9, 3, 1], [3, 9, 25, 9, 3], [1, 3, 9, 3, 1],
#                   [1, 3, 5, 3, 1]])
#kernel = np.array([[ -3-3j, 0-10j,  +3 -3j], [-10+0j, 0+ 0j, +10 +0j], [ -3+3j, 0+10j,  +3 +3j]])

smooth_image = signal.convolve2d(image, kernel, mode='same', boundary='symm', fillvalue=np.nan)
image = smooth_image/64
#smooth_image = convolve(image, kernel, boundary='fill', fill_value=0.0, nan_treatment='interpolate',
#        normalize_kernel=True, mask=None, preserve_nan=True)
#image = smooth_image

# Compute the Hessian matrix
x_deriv_kernel = np.atleast_2d([-1, 0, 1])
y_deriv_kernel = x_deriv_kernel.T


Lx = signal.convolve2d(image, x_deriv_kernel, mode='same', boundary='symm', fillvalue=np.nan)
Ly = signal.convolve2d(image, y_deriv_kernel, mode='same', boundary='symm', fillvalue=np.nan)

Lxy=signal.convolve2d(Lx, y_deriv_kernel, mode='same', boundary='symm', fillvalue=np.nan)
Lxx=signal.convolve2d(Lx, x_deriv_kernel, mode='same', boundary='symm', fillvalue=np.nan)
Lyy=signal.convolve2d(Ly, y_deriv_kernel, mode='same', boundary='symm', fillvalue=np.nan)

"""
Lx = convolve(image, x_deriv_kernel,boundary='fill', fill_value=0.0,
        nan_treatment='interpolate', normalize_kernel=False, mask=None, preserve_nan=True)

Ly = convolve(image, y_deriv_kernel,boundary='fill', fill_value=0.0,
        nan_treatment='interpolate', normalize_kernel=False, mask=None, preserve_nan=True)

Lxy = convolve(Lx, y_deriv_kernel, boundary='fill', fill_value=0.0, nan_treatment='interpolate',
        normalize_kernel=False, mask=None, preserve_nan=True)
Lxx = convolve(Lx, x_deriv_kernel, boundary='fill', fill_value=0.0, nan_treatment='interpolate',
        normalize_kernel=False, mask=None, preserve_nan=True)
Lyy = convolve(Ly, y_deriv_kernel, boundary='fill', fill_value=0.0, nan_treatment='interpolate',
        normalize_kernel=False, mask=None, preserve_nan=True)
"""
first_deri = Lxy
second_deri = Lxy

eig_val = np.full((np.shape(image)[0], np.shape(image)[1], 2), np.nan)
eig_vec = np.full((np.shape(image)[0], np.shape(image)[1], 2, 2), np.nan)

ridge = np.full(np.shape(image), np.nan)
x1_val = np.full(np.shape(image), np.nan)
y1_val = np.full(np.shape(image), np.nan)
x2_val = np.full(np.shape(image), np.nan)
y2_val = np.full(np.shape(image), np.nan)
diff = np.full(np.shape(image), np.nan)

count = 0
for i in range(np.shape(image)[0]):
    for j in range(np.shape(image)[1]):

        result = np.linalg.eig([[Lxx[i, j], Lxy[i, j]], [Lxy[i, j], Lyy[i, j]]])

        eig_val[i, j] = result[0]

        eig_vec[i, j] = result[1]
        second_deri[i, j] = sorted(result[0], key=abs)[0]

        step = 1.00
        x1_val[i, j] = i + eig_vec[i, j][0, 0] * step
        y1_val[i, j] = j + eig_vec[i, j][0, 1] * step

        x2_val[i, j] = i - eig_vec[i, j][0, 0] * step
        y2_val[i, j] = j - eig_vec[i, j][0, 1] * step

        xx1 = interp_func(x1_val[i, j], y1_val[i, j])
        xx2 = interp_func(x2_val[i, j], y2_val[i, j])
        diff[i, j] = (xx1 - image[i, j]) * (image[i, j] - xx2)

        if(diff[i, j] < 0 and second_deri[i, j] < 0):
            ridge[i, j] = (second_deri[i, j])**2

#lambda_max_idx = np.argmax(abs(lambda_0), axis=0)

idx = np.full(len(y_arr), np.nan)

for i in range(len(y_arr)):
    try:
        temp = np.where(ridge[i, :] == np.nanmax(ridge[i, :]))[0]
        idx[i] = temp[0]
    except:
        pass
idx = idx.astype(int)
idx[np.where(idx < 0)] = 0

cmap=plt.cm.Spectral

# plt.close('all')
fig, axs1 = plt.subplots(1, 1, figsize=(8, 6))
im1 = axs1.imshow(abs(image.T), extent=[-15, 15, -15, 15], origin='lower', cmap=cmap)
divider1 = make_axes_locatable(axs1)

axs1.plot(y_arr, z_arr[idx], 'kd', ms=1, lw=2)
patch = patches.Circle((0, 0), radius=15, transform=axs1.transData)
im1.set_clip_path(patch)


axs1.set_xlabel(r'Y [$R_\oplus$]', fontsize=18)
axs1.set_ylabel(r'Z [$R_\oplus$]', fontsize=18)

# Define the location of the colorbar, it's size relative to main figure and the padding
# between the colorbar and the figure, the orientation the colorbar
cax1 = divider1.append_axes("top", size="5%", pad=0.01 )
cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)
cbar1.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                     labelbottom=False, pad=0.01)
cbar1.ax.xaxis.set_label_position('top')
cbar1.ax.set_xlabel('Shear Angle', fontsize=18)
#fig.show()

fig_name = f'../figures/ridge_plot.pdf'
plt.savefig( fig_name, bbox_inches='tight', pad_inches = 0.05, format='pdf', dpi=300)
plt.close()

print(f'Took {round(time.time() - start, 3)} seconds')
