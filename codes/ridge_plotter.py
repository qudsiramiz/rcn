from scipy import signal, interpolate
from skimage import data
from skimage import color
from skimage.filters import meijering, sato, frangi, hessian
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import h5py as hf
import datetime
import time

start = time.time()

# Set the fontstyle to Times New Roman
font = { 'family' : 'sans-serif', 'weight' : 'normal', 'size' : 10 }
plt.rc( 'font', **font )
plt.rc('text', usetex=False)


def identity(image, **kwargs):
    """Return the original image, ignoring any kwargs."""
    return image

#y, z, rx_en, shear, va_cs, bisec = rx_model_batch(trange=trange)
model_type = 't96'  # t01
dr = 0.25  # Resolution of model run in R_E units
mp = 0.5  # Magnetopause thichkness
fn = f'../data/all_data_rx_model_{dr}re_{mp}mp_{model_type}_2021-10-20.h5'
dat = hf.File(fn)

image = np.transpose(dat["shear"][:])

cmap = plt.cm.Spectral

kwargs = {'sigmas': [2.2], 'mode': 'nearest', 'alpha': 1.1}

fig, axes = plt.subplots(1, 2)
for i, black_ridges in enumerate([0]):
    for j, func in enumerate([identity, meijering]):
        kwargs['black_ridges'] = black_ridges
        result = func(image, **kwargs)

        axes[j].imshow(result, extent=[0, 120, 0, 120], cmap=cmap, aspect='auto', origin="lower")
        if i == 0:
            axes[j].set_title(['Original\nimage', 'Frangi\nvesselness'
                                  ][j])
        if j == 0:
            axes[j].set_ylabel('black-ridges = ' + str(bool(black_ridges)))
        axes[j].set_xticks([])
        axes[j].set_yticks([])

y_val = np.full(121, np.nan)
for i in range(121):
    y_val[i] = np.argmax(result[:,i])
axes[0].plot(np.linspace(0, 120, 121), y_val, 'k-', ms=2 )

plt.tight_layout()
plt.show()

print(f'Took {round(time.time() - start, 3)} seconds')
