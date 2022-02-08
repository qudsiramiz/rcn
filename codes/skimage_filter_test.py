from skimage import data
from skimage import color
from skimage.filters import meijering, sato, frangi, hessian
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np

shear = shear
rx_en = rx_en
va_cs = va_cs
bisec_msp = bisec_msp

font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=False)


def identity(image, **kwargs):
    """Return the original image, ignoring any kwargs."""
    return image


#image = color.rgb2gray(shear)
image_rotated = np.transpose(bisec_msp)
image_smooth = sp.ndimage.filters.gaussian_filter(image_rotated, order=0, sigma=[8, 8], mode="nearest")
image = image_smooth
cmap = plt.cm.plasma

kwargs = {'sigmas': [2], 'mode': 'reflect'}

fig = plt.figure(num=None, figsize=(25, 25), dpi=200, facecolor='w', edgecolor='gray')
#fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

fig, axes = plt.subplots(2, 5)
for i, black_ridges in enumerate([1, 0]):
    for j, func in enumerate([identity, meijering, sato, frangi, hessian]):
        kwargs['black_ridges'] = black_ridges
        result = func(image, **kwargs)
        axes[i, j].imshow(result, cmap=cmap, aspect='equal', origin='lower')
        if i == 0:
            axes[i, j].set_title(['Original\nimage', 'Meijering\nneuriteness',
                                  'Sato\ntubeness', 'Frangi\nvesselness',
                                  'Hessian\nvesselness'][j])
        if j == 0:
            axes[i, j].set_ylabel('black_ridges = ' + str(bool(black_ridges)))
        axes[i, j].set_xticks([])
        axes[i, j].set_yticks([])

plt.tight_layout()
plt.savefig('../figures/skimage_filter_test.pdf', bbox_inches='tight', dpi=300, transparent=True, pad_inches=0.05, bbox_extra_artists=[])
plt.close()