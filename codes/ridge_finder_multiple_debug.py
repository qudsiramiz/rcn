
def ridge_finder_multiple_debug(
    image=[None, None, None, None],
    t_range=['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
    xrange=[-15.1, 15],
    yrange=[-15.1, 15],
    mms_probe_num='1',
    mms_sc_pos= [0, 0],
    dr=0.5,
    dipole_tilt_angle=None,
    imf_clock_angle=None,
    sigma=[2.2, 2.2, 2.2, 2.2],
    mode="nearest",
    alpha=1.,
    vmin=[None, None, None, None],
    vmax=[None, None, None, None],
    cmap=["viridis", "viridis", "viridis", "viridis"],
    draw_patch=[True, True, True, True],
    draw_ridge=[True, True, True, True],
    save_fig=True,
    fig_name="new",
    fig_format="png",
    c_label=[None, None, None, None],
    c_unit=[None, None, None, None],
    ):
    r"""
    Finds ridges in an image and plot the points with maximum ridge value on the given image.

    Parameters
    ----------
    image : ndarray
            The image to find ridges in. Default is None.
    trace_range : list of str
            The time range to find the ridge in. Default is ['2016-12-24 15:08:00',
            '2016-12-24 15:12:00'].
    xrange : list of floats, optional
            The range of x-values for image. Default is [-15.1, 15].
    yrange : list of floats, optional
            The range of y-values for image. Default is [-15.1, 15].
    mms_probe_num : str, optional
            The probe number of the MMS spacecraft. Default is 1.
    mms_sc_pos : list of floats, optional
            The position of the spacecraft in the image. Default is [0, 0].
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

    fig = plt.figure(num=None, figsize=(6, 13 ), dpi=200, facecolor='w', edgecolor='gray')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01,  wspace=0.02, hspace=0.)
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1])

    print(len(image))
    y_vals = []
    for i in range(len(image)):
        image_rotated = np.transpose(image[i])

        if cmap is None:
            cmap = "viridis"
        else:
            cmap = cmap[i]
        if(vmin is not None and vmax is not None):
            norm = plt.Normalize(vmin=vmin[i], vmax=vmax[i])
        else:
            norm = plt.Normalize()

        kwargs = {'sigmas': [3], 'black_ridges': False, 'mode': mode, 'alpha': 1}

        # Smoothen the image
        image_smooth = sp.ndimage.filters.gaussian_filter(image_rotated, sigma=[5, 5], mode=mode)
        result = meijering(image_smooth, **kwargs)

        x_len = image_rotated.shape[0]
        y_len = image_rotated.shape[1]

        y_val = np.full(y_len, np.nan)
        y_vals.append(y_val)
        im_max_val = np.full(y_len, np.nan)
        for xx in range(y_len):
            y_val[xx] = np.argmax(result[:, xx]) * dr + yrange[0]
            im_max_val[xx] = np.argmax(image_rotated[:, xx]) * dr + yrange[0]

        # plt.close('all')
        # TODO: Find a better way to do this
        if i==0:
            j = 0
            k = 0
        elif i==1:
            j = 0
            k = 1
        elif i==2:
            j = 1
            k = 0
        elif i==3:
            j = 1
            k = 1

        axs1 = plt.subplot(gs[j, k])
        im1 = axs1.imshow(image_smooth, extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                          origin='lower', cmap=cmap, norm=norm)
        divider1 = make_axes_locatable(axs1)

        # Take rolling average of the y_val array
        y_val_avg = np.full(len(y_val), np.nan)
        for i in range(len(y_val)):
            y_val_avg[i] = np.nanmean(y_val[max(0, i-5):min(len(y_val), i+5)])
    
        if draw_ridge:
            axs1.plot(np.linspace(xrange[0], xrange[1], x_len), y_val_avg, 'k-', alpha=0.9)
            axs1.plot(np.linspace(xrange[0], xrange[1], x_len), im_max_val, 'k*', ms=1, alpha=0.5)

        # Plot a horizontal line at x=0 and a vertical line at y=0
        axs1.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
        axs1.axvline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)

        if(draw_patch):
            patch = patches.Circle((0, 0), radius=(xrange[1] - xrange[0])/2.,
                                    transform=axs1.transData, fc='none', ec='k', lw=0.1)
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
        cbar1.ax.set_xlabel(f'{c_label[i]} ({c_unit[i]})', fontsize=18)

        # Draw the spacecraft position
        axs1.plot(mms_sc_pos[0], mms_sc_pos[1], 'k', marker=mms_probe_num, ms=10, alpha=1)

        # Write the timme range on the plot
        axs1.text(1.0, 0.5, f'Time range: {t_range[0]} - {t_range[1]}', horizontalalignment='left',
                  verticalalignment='center', transform=axs1.transAxes, rotation=270, color='r')
        axs1.text(0.01, 0.99, f'Clock Angle: {np.round(imf_clock_angle, 2)}$^\circ$',
        horizontalalignment='left', verticalalignment='top', transform=axs1.transAxes, rotation=0,
        color='r')
        axs1.text(0.99, 0.99, f'Dipole tilt: {np.round(dipole_tilt_angle * 180/np.pi, 2)} $^\circ$',
                  horizontalalignment='right', verticalalignment='top', transform=axs1.transAxes,
                rotation=0, color='r')
    # fig.show()

    if save_fig:
        try:
            fig_time_range = f"{parser.parse(t_range[0]).strftime('%Y-%m-%d_%H-%M-%S')}_{parser.parse(t_range[1]).strftime('%Y-%m-%d_%H-%M-%S')}"
            fig_name = f'../figures/{fig_name}/ridge_plot_{fig_name}_{fig_time_range}.{fig_format}'
            plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format=fig_format, dpi=300)
            print(f'Figure saved as {fig_name}')
        except  Exception as e:
            print(e)
            print(f'Figure not saved, folder does not exist. Create folder ../figures')
            #pass
        plt.close()
    return y_vals
