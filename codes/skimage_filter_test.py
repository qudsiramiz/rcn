from skimage import data
from skimage import color
from skimage.filters import meijering, sato, frangi, hessian
import matplotlib.pyplot as plt

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
image = np.transpose(va_cs)
cmap = plt.cm.plasma

kwargs = {'sigmas': [1], 'mode': 'reflect'}

fig, axes = plt.subplots(2, 5)
for i, black_ridges in enumerate([1, 0]):
    for j, func in enumerate([identity, meijering, sato, frangi, hessian]):
        kwargs['black_ridges'] = black_ridges
        result = func(image, **kwargs)
        axes[i, j].imshow(result, cmap=cmap, aspect='auto')
        if i == 0:
            axes[i, j].set_title(['Original\nimage', 'Meijering\nneuriteness',
                                  'Sato\ntubeness', 'Frangi\nvesselness',
                                  'Hessian\nvesselness'][j])
        if j == 0:
            axes[i, j].set_ylabel('black_ridges = ' + str(bool(black_ridges)))
        axes[i, j].set_xticks([])
        axes[i, j].set_yticks([])

plt.tight_layout()
plt.savefig('../figures/skimage_filter_test.pdf', bbox_inches='tight', dpi=300, transparent=True, pad_inches=0.05, frameon=False, bbox_extra_artists=[])
plt.close()