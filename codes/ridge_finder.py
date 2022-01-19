import time
import numpy as np
import h5py as hf
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import signal, interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
kernel = np.array([[1, 3, 3, 1], [3, 9, 9, 3], [3, 9, 9, 3], [1, 3, 3, 1]])
# kernel = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]])

smooth_image = signal.convolve2d(image, kernel, mode='full', boundary='symm', fillvalue=np.nan)

# Compute the Hessian matrix
x_deriv_kernel = np.atleast_2d([-1, 0, 1])
y_deriv_kernel = x_deriv_kernel.T

Lx = signal.convolve2d(smooth_image, x_deriv_kernel, mode='same', boundary='symm', fillvalue=np.nan)
Ly = signal.convolve2d(smooth_image, y_deriv_kernel, mode='same', boundary='symm', fillvalue=np.nan)

Lxy=signal.convolve2d(Lx, y_deriv_kernel)
Lxx=signal.convolve2d(Lx, x_deriv_kernel)
Lyy=signal.convolve2d(Ly, y_deriv_kernel)

first_deri = Lxy
second_deri = Lxy

eig_val = np.full((np.shape(image)[0], np.shape(image)[1], 2), np.nan)
eig_vec =  np.full((np.shape(image)[0], np.shape(image)[1], 2, 2), np.nan)

ridge = np.full(np.shape(image), np.nan)
x1_val = np.full(np.shape(image), np.nan)
y1_val = np.full(np.shape(image), np.nan)
x2_val = np.full(np.shape(image), np.nan)
y2_val = np.full(np.shape(image), np.nan)

for i in range(np.shape(image)[0]):
    for j in range(np.shape(image)[1]):

        result = np.linalg.eig([[Lxx[i, j], Lxy[i, j]], [Lxy[i, j], Lyy[i, j]]])

        eig_val[i, j] = result[0]
        eig_vec[i, j] = result[1]
        second_deri[i, j] = np.min(result[0])
        step = 1.00
        x1_val[i, j] = i + eig_vec[i, j][0, 0] * step
        y1_val[i, j] = j + eig_vec[i, j][0, 1] * step

        x2_val[i, j] = i - eig_vec[i, j][0, 0] * step
        y2_val[i, j] = j - eig_vec[i, j][0, 1] * step

        xx1 = interp_func(x1_val[i, j], y1_val[i, j])
        xx2 = interp_func(x2_val[i, j], y2_val[i, j])
        diff = (xx1 - image[i, j]) * (image[i, j] - xx2)

        if(diff < 0 and second_deri[i, j] < 0):
            ridge[i, j] = (second_deri[i, j])**2

idx = np.full(len(y_arr), np.nan)
for i in range(len(y_arr)):
    temp = np.where(ridge[:, i] == np.nanmax(ridge[:, i]))[0]
    idx[i] = temp[0]
idx = idx.astype(int)

# plt.close('all')
fig, (axs1, axs2) = plt.subplots(1, 1, figsize=(8, 6))
im1 = axs1.imshow(image.T, extent=[-15, 15, -15, 15], origin='lower', cmap=plt.cm.Spectral)
divider1 = make_axes_locatable(axs1)


axs1.plot(y_arr, z_arr[idx], 'k-', ms=1, lw=2)
# Define the location of the colorbar, it's size relative to main figure and the padding
# between the colorbar and the figure, the orientation the colorbar
cax1 = divider1.append_axes("top", size="5%", pad=0.01 )
cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)

cbar1.ax.xaxis.set_label_position('top')

im2 = axs2.imshow(ridge.T, origin='lower', cmap=plt.cm.Spectral,
        norm=mpl.colors.LogNorm(vmin=100, vmax=1e3))
divider2 = make_axes_locatable(axs2)

# Define the location of the colorbar, it's size relative to main figure and the padding
# between the colorbar and the figure, the orientation the colorbar
cax2 = divider2.append_axes("top", size="5%", pad=0.01 )
cbar2 = plt.colorbar(im2, cax=cax2, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)

cba
cbar2.ax.xaxis.set_label_position('top')


fig.show()

print(f'Took {round(time.time() - start, 3)} seconds')

