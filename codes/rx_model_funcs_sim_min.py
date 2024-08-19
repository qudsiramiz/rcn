import datetime
import multiprocessing as mp
import os
import matplotlib

import geopack.geopack as gp
import h5py as hf
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from matplotlib.pyplot import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import frangi


def ridge_finder_multiple(
    image=[None, None, None, None],
    convolution_order=[1, 1, 1, 1],
    t_range=["2016-12-24 15:08:00", "2016-12-24 15:12:00"],
    dt=5,
    b_imf=[-5, 0, 0],
    b_msh=[-5, 0, 0],
    v_msh=[-200, 50, 50],
    xrange=[-15.1, 15],
    yrange=[-15.1, 15],
    mms_probe_num="1",
    mms_sc_pos=[0, 0],
    dr=0.5,
    dipole_tilt_angle=None,
    p_dyn=None,
    imf_clock_angle=None,
    sigma=[2.2, 2.2, 2.2, 2.2],
    mode="nearest",
    alpha=1.0,
    vmin=[None, None, None, None],
    vmax=[None, None, None, None],
    cmap_list=["viridis", "viridis", "viridis", "viridis"],
    draw_patch=[True, True, True, True],
    draw_ridge=[False, False, False, False],
    save_fig=True,
    fig_name="new",
    fig_format="png",
    c_label=[None, None, None, None],
    c_unit=[None, None, None, None],
    wspace=0.1,
    hspace=0.1,
    fig_size=(10, 10),
    box_style=None,
    title_y_pos=0.95,
    interpolation="nearest",
    tsy_model="t96",
    dark_mode=True,
    rc_file_name="rc_file.csv",
    rc_folder="../data",
    save_rc_file=False,
    walen1=False,
    walen2=False,
    jet_detection=False,
    fig_version="v6",
):
    r"""
    Finds ridges in an image and plot the points with maximum ridge value on the given image.

    Parameters
    ----------
    image : list of numpy arrays
        List of images to be plotted.
    convolution_order : list of ints
        List of the order of the convolution to be used while smoothing each image. Values must be
        non-negative integers. Default is [1, 1, 1, 1].
    t_range : list of str
            The time range to find the ridge in. Default is ['2016-12-24 15:08:00',
            '2016-12-24 15:12:00'].
    dt : float, optional
        The time differential, in minutes, for observation if 't_range' has only one element.
        Default is 5 minutes.
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
    dipole_tilt_angle : float, optional
            The dipole tilt angle. Default is None.
    imf_clock_angle : float, optional
            The IMF clock angle. Default is None.
    sigma : list of floats, optional
            List of sigmas to be used for the ridge plot. Default is [2.2, 2.2, 2.2, 2.2].
    mode : str
            The mode of the filter. Can be 'nearest', 'reflect', 'constant', 'mirror', 'wrap' or
            'linear'. Default is 'nearest'.
    alpha : float
            The alpha value for the filter. Default is 0.5.
    vmin : list of floats, optional
            List of vmin values for the ridge plot. Default is [None, None, None, None].
    vmax : list of floats, optional
            List of vmax values for the ridge plot. Default is [None, None, None, None].
    cmap_list : list of str, optional
            List of colormaps to be used for the ridge plot. Default is ['viridis', 'viridis',
            'viridis', 'viridis'].
    draw_patch : list of bool, optional
            Whether to draw the circular patch. Default is [True, True, True, True].
    draw_ridge : list of bool, optional
            Whether to draw the ridge line. Default is [True, True, True, True].
    save_fig : bool, optional
            Whether to save the figure. Default is True.
    fig_name : str, optional
            The name of the figure. Default is "new".
    fig_format : str, optional
            The format of the figure. Default is "pdf".
    c_label : list of str, optional
            List of colorbar labels. Default is [None, None, None, None].
    c_unit : list of str, optional
            List of colorbar units. Default is [None, None, None, None].
    wspace : float, optional
            The width space between subplots. Default is 0.1.
    hspace : float, optional
            The height space between subplots. Default is 0.1.
    fig_size : tuple of floats, optional
            The size of the figure. Default is (10, 10).
    box_style : dict, optional
            The style of the box. Default is None.
    title_y_pos : float, optional
            The y-position of the title. Default is 0.95.
    interpolation : str, optional
            The interpolation method for imshow. Default is 'nearest'.
            Options are 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning',
            'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell'
    dark_mode : bool, optional
        Sets the dark mode for the plot and adjusts the color of labels and tickmarks accordingly.
        Default is True.

    Raises
    ------
    ValueError: If the image is not a numpy array.

    Returns
    -------
    ridge_points : ndarray
    """
    if image is None:
        raise ValueError("No image given")

    if len(t_range) == 1:
        # Check if t_range is a datetime object
        if isinstance(t_range[0], datetime.datetime):
            t_range_date = t_range[0]
        else:
            t_range_date = datetime.datetime.strptime(t_range[0], "%Y-%m-%d %H:%M:%S")
        t_range_date_min = t_range_date - datetime.timedelta(minutes=dt)
        t_range_date_max = t_range_date + datetime.timedelta(minutes=dt)
        t_range = [
            t_range_date_min.strftime("%Y-%m-%d %H:%M:%S"),
            t_range_date_max.strftime("%Y-%m-%d %H:%M:%S"),
        ]

    if dark_mode:
        plt.style.use("dark_background")
        # tick_color = 'w'  # color of the tick lines
        mtick_color = "w"  # color of the minor tick lines
        label_color = "w"  # color of the tick labels
        clabel_color = "w"  # color of the colorbar label
    else:
        # plt.style.use("default")
        # tick_color = 'k'  # color of the tick lines
        mtick_color = "k"  # color of the minor tick lines
        label_color = "k"  # color of the tick labels
        clabel_color = "k"  # color of the colorbar label

    # Set the fontstyle to Times New Roman

    font = {"family": "Helvetica", "weight": "normal", "size": 14}
    plt.rc("font", **font)
    plt.rc("text", usetex=True)

    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["font.family"] = "Helvetica"
    matplotlib.rcParams["font.size"] = 15

    fig = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor="w", edgecolor="k")
    fig.subplots_adjust(
        left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=wspace, hspace=hspace
    )
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1])

    # Set the font size for the axes
    label_size = 18  # fontsize for x and y labels
    t_label_size = 18  # fontsize for tick label
    c_label_size = 18  # fontsize for colorbar label
    ct_tick_size = 14  # fontsize for colorbar tick labels
    l_label_size = 14  # fontsize for legend label

    tick_len = 10  # length of the tick lines
    mtick_len = 7  # length of the minor tick lines
    tick_width = 1  # tick width in points
    mtick_width = 0.7  # minor tick width in points

    # box_style = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    if dark_mode:
        box_style = box_style
    else:
        box_style = dict(boxstyle="round", color="w", alpha=0.8, linewidth=1)
    y_vals = []
    x_intr_vals_list = []
    y_intr_vals_list = []
    for i in range(len(image)):
        image_rotated = np.transpose(image[i])

        # Create the masked image from result for all the new processings
        # Find the number of rows in the original image
        n_rows, n_cols = image_rotated.shape

        # Make a grid of the data based on mumber of rows and columns
        X, Y = np.ogrid[:n_rows, :n_cols]

        # Find the central row and column
        c_row = int(n_rows / 2)
        c_col = int(n_cols / 2)
        # Find the distance of each pixel from the central pixel in terms of pixels
        dist_pxl = np.sqrt((X - c_row) ** 2 + (Y - c_col) ** 2)
        mask_image = dist_pxl > 15 / dr

        if cmap_list is None:
            cmap_list = ["viridis", "viridis", "viridis", "viridis"]
        else:
            cmap_list = cmap_list
        if vmin is not None and vmax is not None:
            norm = plt.Normalize(vmin=vmin[i], vmax=vmax[i])
        else:
            norm = plt.Normalize()

        kwargs = {"sigmas": [sigma[i]], "black_ridges": False, "mode": mode, "alpha": 1}

        # Smoothen the image
        image_smooth = sp.ndimage.gaussian_filter(
            image_rotated, order=convolution_order[i], sigma=[5, 5], mode=mode
        )
        image_smooth_p = sp.ndimage.gaussian_filter(
            image_rotated, order=0, sigma=[5, 5], mode=mode
        )
        result = frangi(image_smooth, **kwargs)  # frangi, hessian, meijering, sato

        m_result = result.copy()
        m_result[mask_image] = np.nan
        new_image_rotated = image_rotated.copy()
        new_image_rotated[mask_image] = np.nan

        x_len = image_rotated.shape[0]
        y_len = image_rotated.shape[1]

        y_val = np.full(y_len, np.nan)
        y_vals.append(y_val)
        im_max_val = np.full(y_len, np.nan)
        for xx in range(y_len):
            try:
                y_val[xx] = np.nanargmax(m_result[:, xx]) * dr + yrange[0]
                im_max_val[xx] = np.nanargmax(new_image_rotated[:, xx]) * dr + yrange[0]
            except Exception:
                pass

        axs1 = plt.subplot(gs[0, i])
        im1 = axs1.imshow(
            image_smooth_p,
            extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
            origin="lower",
            cmap=cmap_list[i],
            norm=norm,
            interpolation=interpolation,
            alpha=1,
        )
        divider1 = make_axes_locatable(axs1)
        # Draw a circle of radius 10 around the center of the image
        # axs1.add_patch(plt.Circle((0, 0), radius=15, color='gray', fill=False, lw=0.5))

        # Take rolling average of the y_val array
        y_val_avg = np.full(len(y_val), np.nan)
        im_max_val_avg = np.full(len(y_val), np.nan)

        r_a_l = 5
        for xx in range(len(y_val)):
            y_val_avg[xx] = np.nanmean(
                y_val[max(0, xx - r_a_l) : min(len(y_val), xx + r_a_l)]
            )
            im_max_val_avg[xx] = np.nanmean(
                im_max_val[max(0, xx - r_a_l) : min(len(y_val), xx + r_a_l)]
            )

        axs1.axhline(0, color="k", linestyle="-", linewidth=0.5, alpha=0.5)
        axs1.axvline(0, color="k", linestyle="-", linewidth=0.5, alpha=0.5)

        if draw_ridge:
            # axs1.plot(np.linspace(xrange[0], xrange[1], x_len), y_val_avg, color='aqua', ls='-',
            #          alpha=0.9)
            x_intr_vals = np.linspace(xrange[0], xrange[1], x_len)
            y_intr_vals = im_max_val_avg
            # If the square root of the sum of squares of x_intr_vals and y_intr_vals is greater
            # than 15, then mask those values
            # r_intr_vals = np.sqrt(x_intr_vals ** 2 + y_intr_vals ** 2)
            # mask = r_intr_vals > 15
            # Mask the values of x_intr_vals and y_intr_vals
            # x_intr_vals[mask] = np.nan
            # y_intr_vals[mask] = np.nan
            # if z component of b_imf is negative, then the ridge is on the left side of the
            # image
            # if b_imf[2] <= 0:
            #    axs1.plot(x_intr_vals, y_intr_vals, color='aqua', ls='-', alpha=0.9)
        # Plot a horizontal line at x=0 and a vertical line at y=0
        if draw_patch:
            patch = patches.Circle(
                (0, 0),
                radius=xrange[1],
                transform=axs1.transData,
                fc="none",
                ec="k",
                lw=0.5,
            )
            # im1.set_clip_path(patch)
        # axs1.add_patch(patch)
        if i == 0 or i == 3:
            axs1.set_ylabel(
                r"Z [GSM, $R_{\rm E}$]", fontsize=label_size, color=label_color
            )
        if i == 3:
            axs1.yaxis.set_label_position("right")

        axs1.set_xlabel(r"Y [GSM, $R_{\rm E}$]", fontsize=label_size, color=label_color)
        if dark_mode:
            text_color = "white"
        else:
            text_color = "black"
        if i == 0:
            # axs1.text(-0.3, 1.16, f'Model: {tsy_model}', horizontalalignment='left',
            #           verticalalignment='bottom', transform=axs1.transAxes, rotation=0,
            #           color=text_color, fontsize=l_label_size, bbox=box_style)
            axs1.text(
                -0.3,
                1.16,
                f"Clock Angle: {np.round(imf_clock_angle, 2)}$^\\circ$",
                horizontalalignment="left",
                verticalalignment="bottom",
                transform=axs1.transAxes,
                rotation=0,
                color=text_color,
                fontsize=l_label_size,
                bbox=box_style,
            )

        if i == 3:
            axs1.text(
                1.3,
                1.16,
                f"$B_{{\\rm {{imf}}}}$ = [{b_imf[0]}, {b_imf[1]}, {b_imf[2]}]",
                horizontalalignment="right",
                verticalalignment="bottom",
                transform=axs1.transAxes,
                rotation=0,
                color=text_color,
                fontsize=l_label_size,
                bbox=box_style,
            )
        # elif i == 3:
        #     axs1.text(1.3, -0.15,
        #               f'Dipole tilt: {np.round(dipole_tilt_angle * 180/np.pi, 2)} ${{\\hspace{{-.2em}}}}^\\circ$',
        #               horizontalalignment='right', verticalalignment='top',
        #               transform=axs1.transAxes, rotation=0, color=text_color, fontsize=l_label_size,
        #               bbox=box_style)

        # Define the location of the colorbar, it's size relative to main figure and the padding
        # between the colorbar and the figure, the orientation the colorbar
        cax1 = divider1.append_axes("top", size="5%", pad=0.01)
        cbar1 = plt.colorbar(
            im1, cax=cax1, orientation="horizontal", ticks=None, fraction=0.05, pad=0.01
        )
        cbar1.ax.tick_params(
            axis="x",
            direction="in",
            top=True,
            labeltop=True,
            bottom=False,
            labelbottom=False,
            pad=0.01,
            labelsize=ct_tick_size,
            labelcolor=label_color,
        )
        # Get the location of all ticks on the colorbar
        cbar_ticks = cbar1.ax.get_xticks()
        # Remove the first tick
        cbar_ticks = cbar_ticks[1:]
        # Set the ticks to the new tick values
        cbar1.ax.set_xticks(cbar_ticks)

        cbar1.ax.xaxis.set_label_position("top")

        cbar1.ax.set_xlabel(f"{c_label[i]}", fontsize=c_label_size, color=clabel_color)

        # Set tick label parameters
        if i == 0:
            axs1.tick_params(
                axis="both",
                direction="in",
                which="major",
                left=True,
                right=True,
                top=True,
                bottom=True,
                labelleft=True,
                labelright=False,
                labeltop=False,
                labelbottom=True,
                labelsize=t_label_size,
                length=tick_len,
                width=tick_width,
                labelcolor=label_color,
            )
        elif i == 1 or i == 2:
            axs1.tick_params(
                axis="both",
                direction="in",
                which="major",
                left=True,
                right=True,
                top=True,
                bottom=True,
                labelleft=False,
                labelright=False,
                labeltop=False,
                labelbottom=True,
                labelsize=t_label_size,
                length=tick_len,
                width=tick_width,
                labelcolor=label_color,
            )
        else:
            axs1.tick_params(
                axis="both",
                direction="in",
                which="major",
                left=True,
                right=True,
                top=True,
                bottom=True,
                labelleft=False,
                labelright=True,
                labeltop=False,
                labelbottom=True,
                labelsize=t_label_size,
                length=tick_len,
                width=tick_width,
                labelcolor=label_color,
            )

        if i == 0:
            # Add a label '(a)' to the plot to indicate the panel number
            axs1.text(
                0.05,
                0.15,
                "(a)",
                horizontalalignment="left",
                verticalalignment="top",
                transform=axs1.transAxes,
                rotation=0,
                color=text_color,
                fontsize=1.2 * l_label_size,
            )
        elif i == 1:
            # Add a label '(b)' to the plot to indicate the panel number
            axs1.text(
                0.05,
                0.15,
                "(b)",
                horizontalalignment="left",
                verticalalignment="top",
                transform=axs1.transAxes,
                rotation=0,
                color=text_color,
                fontsize=1.2 * l_label_size,
            )
        elif i == 2:
            # Add a label '(c)' to the plot to indicate the panel number
            axs1.text(
                0.05,
                0.15,
                "(c)",
                horizontalalignment="left",
                verticalalignment="top",
                transform=axs1.transAxes,
                rotation=0,
                color=text_color,
                fontsize=1.2 * l_label_size,
            )
        elif i == 3:
            # Add a label '(d)' to the plot to indicate the panel number
            axs1.text(
                0.05,
                0.15,
                "(d)",
                horizontalalignment="left",
                verticalalignment="top",
                transform=axs1.transAxes,
                rotation=0,
                color=text_color,
                fontsize=1.2 * l_label_size,
            )

        # Show minor ticks
        axs1.minorticks_on()
        axs1.tick_params(
            axis="both",
            which="minor",
            direction="in",
            length=mtick_len,
            left=True,
            right=True,
            top=True,
            bottom=True,
            color=mtick_color,
            width=mtick_width,
        )
        # Set the number of ticks on the x-axis
        axs1.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="lower"))
        # Set the number of ticks on the y-axis
        axs1.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="lower"))

        # Setting the tickmarks labels in such a way that they don't overlap
        plt.setp(axs1.get_xticklabels(), rotation=0, ha="right", va="top", visible=True)
        plt.setp(axs1.get_yticklabels(), rotation=0, va="center", visible=True)
        # Set the title of the plot
        # fig.suptitle(f'$B_{{\\rm {{imf}}}}$ = {b_imf}',
        #              fontsize=label_size, color=text_color, y=title_y_pos, alpha=0.65)

    # plt.show()
    if save_fig:
        try:
            # TODO: Add folder name as one of the path and make sure that the code creates the
            # folder. Gives out error if the folder can't be created.
            fig_folder = (
                f"../figures/test/{tsy_model}/{interpolation}"
                + f"_interpolation_mms{mms_probe_num}/{fig_version}"
            )
            check_folder = os.path.isdir(fig_folder)
            # If folder doesn't exist, then create it.
            if not check_folder:
                os.makedirs(fig_folder)
                print("created folder : ", fig_folder)
            else:
                print(f"folder already exists: {fig_folder}\n")

            # fig_folder = "../figures/test"
            fig_name = f"{fig_folder}/ridge_plot_{int(b_imf[0])}_{int(b_imf[1])}_{int(b_imf[2])}.{fig_format}"
            plt.savefig(
                fig_name, bbox_inches="tight", pad_inches=0.05, format=fig_format, dpi=200
            )
            print(f"Figure saved as {fig_name}")
        except Exception as e:
            print(e)
            print("Figure not saved, folder does not exist. Create folder ../figures")
            # pass
        # plt.close()
    plt.close()
    return y_vals, x_intr_vals_list, y_intr_vals_list
