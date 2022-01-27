def ridge_finder(
    image=None,
    t_range=['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
    xrange=[-15.1, 15],
    yrange=[-15.1, 15],
    dr=0.5,
    sigma=2.2,
    mode="nearest",
    alpha=1.,
    vmin=None,
    vmax=None,
    cmap="viridis",
    draw_patch=False,
    draw_ridge=False,
    save_fig=True,
    fig_name="new",
    fig_format="pdf",
    c_label="none",
    c_unit="none",
    ):
    r"""
    Finds ridges in an image and plot the points with maximum ridge value on the given image.

    Parameters
    ----------
    image : ndarray
            The image to find ridges in. Default is None.
    xrange : list of floats, optional
            The range of x-values for image. Default is [-15.1, 15].
    yrange : list of floats, optional
            The range of y-values for image. Default is [-15.1, 15].
    dr : float, optional
            The step size for the grid. Default is 0.5.
    sigma : float
            The size of the filter. Default is 2.2.
    mode : str
            The mode of the filter. Can be 'nearest', 'reflect', 'constant', 'mirror', 'wrap' or
            'linear'. Default is 'nearest'.
    alpha : float
            The alpha value for the filter. Default is 0.5.
    vmin : float, optional
            The minimum value of the colorbar. Default is None.
    vmax : float, optional
            The maximum value of the colorbar. Default is None.
    cmap : str, optional
            The colormap to use. Default is 'viridis'.
    draw_patch : bool, optional
            Whether to draw the circular patch. Default is False.
    draw_ridge : bool, optional
            Whether to draw the ridge line. Default is False.
    save_fig : bool, optional
            Whether to save the figure. Default is True.
    fig_name : str, optional
            The name of the figure. Default is "new".
    fig_format : str, optional
            The format of the figure. Default is "pdf".
    c_label : str, optional
            The label for the colorbar. Default is "none".
    c_unit : str, optional
            The units for the colorbar label. Default is "none".

    Raises
    ------
    ValueError: If the image is not a numpy array.

    Returns
    -------
    ridge_points : ndarray
    """
    if image is None:
        raise ValueError("No image given")

    # NOTE: This is a hack to ensure that the output of shear angle, reconnection energy etc. agrees
    # with what has been reported in literature. Plot for shear angle seems to agree reasonably well
    # (for "trange = ['2016-12-07 05:11:00', '2016-12-07 05:21:00']") with the one reported by
    # FuselierJGR2019 (doi:10.1029/2019JA027143, see fig. 4).
    image_rotated = np.transpose(np.flipud(np.fliplr(image)))

    if cmap is None:
        cmap = "viridis"
    if(vmin is not None and vmax is not None):
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = plt.Normalize()

    #cmap = plt.cm.jet

    kwargs = {'sigmas': [3], 'black_ridges': False, 'mode': mode, 'alpha': 1}

    # Smoothen the image
    image_smooth = sp.ndimage.filters.gaussian_filter(image_rotated, sigma=[5, 5], mode=mode)
    result = meijering(image_smooth, **kwargs)

    x_len = image_rotated.shape[0]
    y_len = image_rotated.shape[1]

    y_val = np.full(y_len, np.nan)
    im_max_val = np.full(y_len, np.nan)
    for i in range(y_len):
        y_val[i] = np.argmax(result[:, i]) * dr + yrange[0]
        im_max_val[i] = np.argmax(image_rotated[:, i]) * dr + yrange[0]

    # plt.close('all')
    fig, axs1 = plt.subplots(1, 1, figsize=(8, 6))

    im1 = axs1.imshow(image_smooth, extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                      origin='lower', cmap=cmap, norm=norm)
    divider1 = make_axes_locatable(axs1)

    if draw_ridge:
        axs1.plot(np.linspace(xrange[0], xrange[1], x_len), y_val, 'k-', alpha=0.9)
        axs1.plot(np.linspace(xrange[0], xrange[1], x_len), im_max_val, 'k*', ms=1, alpha=0.5)

    # Plot a horizontal line at x=0 and a vertical line at y=0
    axs1.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    axs1.axvline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)

    if(draw_patch):
        patch = patches.Circle((0, 0), radius=(xrange[1] - xrange[0])/2., transform=axs1.transData,
                            fc='none', ec='k', lw=0.1)
        axs1.add_patch(patch)
        im1.set_clip_path(patch)

    axs1.set_xlabel(r'Y [GSM, $R_\oplus$]', fontsize=18)
    axs1.set_ylabel(r'Z [GSM, $R_\oplus$]', fontsize=18)

    # Define the location of the colorbar, it's size relative to main figure and the padding
    # between the colorbar and the figure, the orientation the colorbar
    cax1 = divider1.append_axes("top", size="5%", pad=0.01)
    cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05,
                         pad=0.01)
    cbar1.ax.tick_params(axis="x", direction="in", top=True, labeltop=True, bottom=False,
                         labelbottom=False, pad=0.01)
    cbar1.ax.xaxis.set_label_position('top')
    cbar1.ax.set_xlabel(f'{c_label} ({c_unit})', fontsize=18)

    # Write the timme range on the plot
    axs1.text(1.0, 0.5, f'Time range: {trange[0]} - {trange[1]}', horizontalalignment='left',
              verticalalignment='center', transform=axs1.transAxes, rotation=270, color='r')

    # fig.show()

    if save_fig:
        try:
            fig_name = f'../figures/ridge_plot_vir_{fig_name}_{dr}dr_{m_p}mp_{t_range[0][:10]}_{t_range[0][-8:]}_{t_range[1][:10]}_{t_range[1][-8:]}.{fig_format}'
            plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=fig_format, dpi=300)
            print(f'Figure saved as {fig_name}')
        except  Exception as e:
            print(e)
            print(f'Figure not saved, folder does not exist. Create folder ../figures')
            #pass
        plt.close()
    return y_val
