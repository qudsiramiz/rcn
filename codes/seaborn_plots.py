import seaborn as sns; sns.set()
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

import SeabornFig2Grid as sfg


def kde_plots(
              df=None,
              x=None,
              y=None,
              x_label="",
              y_label="",
              log_scale=True,
              x_scale="lin",
              y_scale="lin",
              xlim=[0, 2e1],
              ylim=[1e-1, 2e2],
              color="blue",
              spearman=None,
              pearson=None,
              cmap="Blues",
              fig_save=True,
              hue_norm=mpl.colors.LogNorm(),
              height=8,
              ratio=8,
              space=0,
              alpha=0.7,
              ):

    pad = 7
    labelsize = 45
    ticklabelsize = 38
    clabelsize = 30
    ticklength = 10
    
    #plt.subplots_adjust(left=.1, right=.9, top=0.9, bottom=0.1)

    # Find indices of the data where the data is not nan
    # ind_not_nan_x = np.where(np.logical_not(np.isnan(x_data)))[0]
    # ind_not_nan_y = np.where(np.logical_not(np.isnan(y_data)))[0]
    # # Find the intersection of the two sets of indices
    # ind_not_nan = np.intersect1d(ind_not_nan_x, ind_not_nan_y)
    # x_data = x_data[ind_not_nan]
    # y_data = y_data[ind_not_nan]

    # axs1 = sns.jointplot(x=x_data, y=y_data, kind=kind, cbar=cbar_status, thresh=thresh,
    #                     fill=fill_status, levels=nlevels, log_scale=log_scale, hue_norm=hue_norm,
    #                     cmap=cmap, xlim=xlim, ylim=ylim, height=height, ratio=ratio, space=space,
    #                     vertical=vertical)
    # axs1 = sns.jointplot(x=x_data, y=y_data, kind=kind, hue_norm=hue_norm,
    #                     cmap=cmap, xlim=xlim, ylim=ylim, height=height, ratio=ratio, space=space)

    # Remove all occurances of delta_beta where delta_beta is less than 0
    if log_scale:
        df = df[df[x] > 0]
        df = df[df[y] > 0]

    axs1 = sns.JointGrid(x=x, y=y, data=df, xlim=xlim, ylim=ylim, height=height, ratio=ratio,
                         space=space, hue_norm=hue_norm)
    axs1.plot_joint(sns.scatterplot, s=25 * df["r_rc"], alpha=alpha)
    axs1.plot_marginals(sns.histplot, kde=True, alpha=alpha, log_scale=log_scale)
    # hue='delta_beta', hue_norm=hue_norm, bins=20, kde=True, stat='density', common_norm=False, common_bins=False, multiple='stack', shrink=.8, alpha=0.8)
    #axs1.plot(sns.scatterplot, sns.histplot)

    if y == "msh_msp_shear":
        shear_angle_theory = np.logspace(-1, np.log10(180), 100)
        delta_beta_theory_half = np.tan(np.deg2rad(shear_angle_theory/2))
        delta_beta_theory_one = 2 * np.tan(np.deg2rad(shear_angle_theory/2))
        delta_beta_theory_two = 4 * np.tan(np.deg2rad(shear_angle_theory/2))

        axs1.fig.axes[0].plot(delta_beta_theory_half, shear_angle_theory, marker='.', c="w",
                              ls=None, lw=0, label=r"$\lambda$ = 0.5")
        axs1.fig.axes[0].plot(delta_beta_theory_one, shear_angle_theory, marker='.', c="b",
                              ls=None, lw=0, label=r"$\lambda$ = 1")
        axs1.fig.axes[0].plot(delta_beta_theory_two, shear_angle_theory, marker='.', c="g",
                              ls=None, lw=0, label=r"$\lambda$ = 2")
        lgnd = axs1.fig.axes[0].legend(loc=4, fontsize=30, frameon=False)
        for handle in lgnd.legendHandles:
            handle.size = [1]

    if spearman is not None:
        x_spearman = np.linspace(0, 25, 100)
        y_spearman = spearman * x_spearman + np.mean(df[y]) - spearman*np.mean(df[x])
        axs1.fig.axes[0].plot(x_spearman, y_spearman, c="w", ls="--", lw=2)
        axs1.fig.axes[0].text(0.02, 0.02, f"$\\rho_{{\\rm {'s'}}}$ = {spearman:.2f}\n"
                                          f"$\\rho_{{\\rm {'p'}}}$ = {pearson:.2f}",
                              transform=plt.gca().transAxes, va="bottom", ha="left",
                              bbox=dict(facecolor='k', alpha=1, edgecolor='k',
                              boxstyle='round,pad=0.2'))

    if ~log_scale and x_scale == "log":
        axs1.fig.axes[0].set_xscale('log')
    
    if ~log_scale and y_scale == "log":
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

    axs1.fig.axes[0].tick_params(axis='both', which='major', direction='in', labelbottom=True,
                                 bottom=True, labeltop=False, top=True, labelleft=True, left=True,
                                 labelright=False, right=True, width=1.5, length=ticklength,
                                 labelsize=ticklabelsize, labelrotation=0, pad=pad)

    axs1.fig.axes[0].tick_params(axis='both', which='minor', direction='in', labelbottom=False,
                                 bottom=False, left=False, width=1.5, length=ticklength,
                                 labelsize=ticklabelsize, labelrotation=0)

    axs1.fig.axes[1].tick_params(axis='both', which='both', direction='in', labelbottom=False,
                                 bottom=False, labelleft=False, left=False, width=1.5,
                                 length=ticklength, labelsize=ticklabelsize, labelrotation=0)

    axs1.fig.axes[2].tick_params(axis='both', which='both', direction='in', labelbottom=False,
                                 bottom=False, labelleft=False, left=False, width=1.5,
                                 length=ticklength, labelsize=ticklabelsize, labelrotation=0)

    #axs1.fig.axes[3].tick_params(axis='y', which='major', direction='in', labelbottom=False,
    #                             bottom=False, labelleft=False, left=False, labelright=True,
    #                             right=True, width=1.5, length=ticklength, labelsize=clabelsize,
    #                             labelrotation=0)

    axs1.set_axis_labels(x_label, y_label, fontsize=labelsize)


    if (fig_save) :
        fname = f'../figures/{x}_vs_{y}.png'
        axs1.savefig(fname, format='png', dpi=400)
    #plt.show()
    plt.close('all')
    return axs1


def seaborn_subplots(
                     df_list=None,
                     keys=[],
                     figsize=(40, 80),
                     labels=[],
                     color_list=[],
                     y_scale="log",
                     x_scale="log",
                     log_scale=True,
                     ):

    axs_list = []
    for i, df in enumerate(df_list):
        # Find the spearman and pearson correlation between key and "r_rc"
        spearman = df[keys[1]].corr(df["r_rc"], method="spearman")
        pearson = df[keys[1]].corr(df["r_rc"], method="pearson")
        
        x_lim = [df[keys[0]].min(), df[keys[0]].max()]
        y_lim = [df[keys[1]].min(), df[keys[1]].max()]
        axs = kde_plots(df=df, x=keys[0], y=keys[1], x_label=labels[0],
                        y_label=labels[1], log_scale=log_scale, x_scale=x_scale, y_scale=y_scale,
                        xlim=x_lim, ylim=y_lim, color=color_list[i],
                        spearman=spearman, pearson=pearson, fig_save=False)
        axs_list.append(axs)

    fig = plt.figure(figsize=(figsize[0], figsize[1]))
    gs = gridspec.GridSpec(2, 2)

    mg0 = sfg.SeabornFig2Grid(axs_list[0] fig, gs[0])
    mg1 = sfg.SeabornFig2Grid(axs_list[1] fig, gs[1])
    mg2 = sfg.SeabornFig2Grid(axs_list[2] fig, gs[2])
    mg3 = sfg.SeabornFig2Grid(axs_list[3] fig, gs[3])

    gs.tight_layout(fig)
    #gs.update(top=0.7)
    plt.savefig("../figures/test.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    #plt.show()
