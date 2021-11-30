import time

import h5py as hf
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import meijering  # sato, frangi, hessian

# Set the fontstyle to Times New Roman
font = {'family': 'sans-serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=False)

start = time.time()

# Load data
model_type = 't96'  # t01
dr = 0.25  # Resolution of model run in R_E units
mp = 0.5  # Magnetopause thichkness
fn = f'../data/all_data_rx_model_{dr}re_{mp}mp_{model_type}_2021-10-20.h5'
dat = hf.File(fn)
image = dat["shear"][:]

def ridge_finder(image=None, sigma=2.2, mode="nearest", alpha=0.5):
        r"""
        Finds ridges in an image and plot the points with maximum ridge value on the given image.

        Parameters
        ----------
        image : ndarray
                The image to find ridges in. Default is None.
        sigma : float
                The size of the filter. Default is 2.2.
        mode : str
                The mode of the filter. Can be 'nearest', 'reflect', 'constant', 'mirror', 'wrap' or
                'linear'. Default is 'nearest'.
        alpha : float
                The alpha value for the filter. Default is 0.5.

        Raises
        ------
        ValueError: If the image is not a numpy array.

        Returns
        -------
        ridge_points : ndarray
        """
        if image is None:
            raise ValueError("No image given")
        image = np.transpose(image)

        cmap = plt.cm.Spectral

        kwargs = {'sigmas': [sigma], 'black_ridges': False, 'mode': mode, 'alpha': alpha}

        result = meijering(image, **kwargs)

        y_val = np.full(121, np.nan)
        for i in range(121):
                y_val[i] = np.argmax(result[:, i])/4 - 15.1


        # plt.close('all')
        fig, axs1 = plt.subplots(1, 1, figsize=(8, 6))

        im1 = axs1.imshow(abs(image), extent=[-15.1, 15, -15.1, 15], origin='lower', cmap=cmap)
        divider1 = make_axes_locatable(axs1)

        axs1.plot(np.linspace(-15.1, 15, 121), y_val, 'k-', ms=2)

        patch = patches.Circle((0, 0), radius=15, transform=axs1.transData)
        im1.set_clip_path(patch)

        axs1.set_xlabel(r'Y [$R_\oplus$]', fontsize=18)
        axs1.set_ylabel(r'Z [$R_\oplus$]', fontsize=18)

        # Define the location of the colorbar, it's size relative to main figure and the padding
        # between the colorbar and the figure, the orientation the colorbar
        cax1 = divider1.append_axes("top", size="5%", pad=0.01 )
        cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05,
                             pad=0.01)
        cbar1.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                             labelbottom=False, pad=0.01)
        cbar1.ax.xaxis.set_label_position('top')
        cbar1.ax.set_xlabel('Shear Angle', fontsize=18)

        fig.show()

        fig_name = f'../figures/ridge_plot_new.pdf'
        plt.savefig( fig_name, bbox_inches='tight', pad_inches = 0.05, format='pdf', dpi=300)
        #plt.close()
        return y_val
print(f'Took {round(time.time() - start, 3)} seconds')

y_val = ridge_finder()
