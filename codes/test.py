from scipy import signal
from scipy import misc
ascent = misc.ascent()
scharr = np.array([[ -3-3j, 0-10j,  +3 -3j],
                   [-10+0j, 0+ 0j, +10 +0j],
                   [ -3+3j, 0+10j,  +3 +3j]]) # Gx + j*Gy
grad = signal.convolve2d(ascent, scharr, boundary='symm', mode='same')


import matplotlib.pyplot as plt
fig, (ax_orig, ax_mag, ax_ang) = plt.subplots(3, 1, figsize=(6, 15))
im1 = ax_orig.imshow(ascent, cmap='gray')
ax_orig.set_title('Original')

im2 = ax_mag.imshow(np.absolute(grad), cmap='gray')
ax_mag.set_title('Gradient magnitude')

im3 = ax_ang.imshow(np.angle(grad), cmap='hsv') # hsv is cyclic, like angles
ax_ang.set_title('Gradient orientation')

cax1 = divider1.append_axes("top", size="5%", pad=0.01 )
cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)
cbar1.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                     labelbottom=False, pad=0.01)
cbar1.ax.xaxis.set_label_position('top')

cax2 = divider1.append_axes("top", size="5%", pad=0.01 )
cbar2 = plt.colorbar(im2, cax=cax2, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)
cbar2.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                     labelbottom=False, pad=0.01)
cbar2.ax.xaxis.set_label_position('top')

cax3 = divider1.append_axes("top", size="5%", pad=0.01 )
cbar3 = plt.colorbar(im2, cax=cax3, orientation='horizontal', ticks=None, fraction=0.05, pad=0.01)
cbar3.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                     labelbottom=False, pad=0.01)
cbar3.ax.xaxis.set_label_position('top')

fig.show()
