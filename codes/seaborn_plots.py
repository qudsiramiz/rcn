import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def kde_plots(x_data,
              y_data,
              data_type=None,
              cmap="Blues",
              fig_save=True,
              kind="scatter",
              nlevels=20,
              cbar_status=True,
              thresh= 0.01,
              fill_status=True,
              log_scale=True,
              hue_norm=mpl.colors.LogNorm(),
              xlim=[1e-3, 2e1],
              ylim=[1e-1, 2e2],
              height=8,
              ratio=10,
              space=0,
              vertical=False,
              ) :

    pad = 0.02
    labelsize = 45
    ticklabelsize = 38
    clabelsize = 30
    ticklength = 10
    
    plt.subplots_adjust(left=.1, right=.9, top=0.9, bottom=0.1)

    # Find indices of the data where the data is not nan
    ind_not_nan_x = np.where(np.logical_not(np.isnan(x_data)))[0]
    ind_not_nan_y = np.where(np.logical_not(np.isnan(y_data)))[0]
    # Find the intersection of the two sets of indices
    ind_not_nan = np.intersect1d(ind_not_nan_x, ind_not_nan_y)
    x_data = x_data[ind_not_nan]
    y_data = y_data[ind_not_nan]

    # axs1 = sns.jointplot(x=x_data, y=y_data, kind=kind, cbar=cbar_status, thresh=thresh,
    #                     fill=fill_status, levels=nlevels, log_scale=log_scale, hue_norm=hue_norm,
    #                     cmap=cmap, xlim=xlim, ylim=ylim, height=height, ratio=ratio, space=space,
    #                     vertical=vertical)
    axs1 = sns.jointplot(x=x_data, y=y_data, kind=kind, hue_norm=hue_norm,
                        cmap=cmap, xlim=xlim, ylim=ylim, height=height, ratio=ratio, space=space)

    axs1.fig.axes[0].set_xscale('log')
    axs1.fig.axes[0].set_yscale('log')
    pos_joint_ax = axs1.ax_joint.get_position()
    pos_marg_x_ax = axs1.ax_marg_x.get_position()
    axs1.ax_joint.set_position([pos_joint_ax.x0, pos_joint_ax.y0, pos_marg_x_ax.width,
                                pos_joint_ax.height])
    axs1.fig.axes[-1].set_position([1, pos_joint_ax.y0, .07, pos_joint_ax.height])

    # get the current colorbar ticks
    cbar_ticks = axs1.fig.axes[-1].get_yticks()
    # get the maximum value of the colorbar
    _, cbar_max = axs1.fig.axes[-1].get_ylim()
    # change the labels (not the ticks themselves) to a percentage
    #axs1.fig.axes[-1].set_yticklabels([f'{t / cbar_max * 1:.3f} %' for t in cbar_ticks], size=clabelsize)


    # axs1.fig.axes[-1].set_xlabel('Density', fontsize=clabelsize, labelpad=10)

    axs1.fig.axes[0].tick_params( axis='both', which='major', direction='in', labelbottom=True,
                                 bottom=True, labeltop=False, top=True, labelleft=True, left=True,
                                 labelright=False, right=True, width=1.5, length=ticklength,
                                 labelsize=ticklabelsize, labelrotation=0 )

    axs1.fig.axes[0].tick_params( axis='both', which='minor', direction='in', labelbottom=False,
                                 bottom=False, left=False, width=1.5, length=ticklength,
                                 labelsize=ticklabelsize, labelrotation=0 )

    axs1.fig.axes[1].tick_params( axis='both', which='both', direction='in', labelbottom=False,
                                 bottom=False, labelleft=False, left=False, width=1.5,
                                 length=ticklength, labelsize=ticklabelsize, labelrotation=0 )

    axs1.fig.axes[2].tick_params( axis='both', which='both', direction='in', labelbottom=False,
                                 bottom=False, labelleft=False, left=False, width=1.5,
                                 length=ticklength, labelsize=ticklabelsize, labelrotation=0 )

    #axs1.fig.axes[3].tick_params( axis='y', which='major', direction='in', labelbottom=False,
    #                             bottom=False, labelleft=False, left=False, labelright=True,
    #                             right=True, width=1.5, length=ticklength, labelsize=clabelsize,
    #                             labelrotation=0 )

    axs1.set_axis_labels(r"$\Delta \beta$", "Shear Angle (${~}^{0}$)", fontsize=labelsize)
    axs1.ax_joint.set_ylabel( r"Shear Angle (${~}^{0}$)", fontsize=labelsize )
    axs1.ax_joint.set_xlabel( r"$\Delta \beta$", fontsize=labelsize )

    if (fig_save) :
        fname = f'../figures/beta_r_kdeplot_{cmap}.pdf'
        axs1.savefig(fname, format='pdf', dpi=400)
    #plt.show()
    plt.close('all')
    return axs1