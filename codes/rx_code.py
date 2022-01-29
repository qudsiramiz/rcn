# This the python version of IDL code named 'RX_model_batch.pro'
import datetime
import time
import h5py as hf
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from rx_model_funcs import *

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

start = time.time()

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

#def rx_model_batch(
#    probe=None,
#    omni_level='hro',
#    maximum_shear=True,
#    mms_probe='mms1',
#    movie=None,
#    trange=None,
#    model_type='t96',
#    m_p=0.5,
#    dr=0.25
#    y_min= -20,
#    y_max= 20,
#    z_min= -20,
#    z_max= 20,
#    draw_patch=False,
#    save_data=False,
#    plot_type="shear"
#):
trange_list = [
                ['2016-12-24 15:08:00', '2016-12-24 15:12:00'],
                ['2016-12-07 05:05:00', '2016-12-07 05:33:00'],
#                ['2015-09-08 11:05:00', '2015-09-08 11:15:00'],
                ['2015-10-16 10:28:00', '2015-10-16 10:38:00'],
                ['2015-10-16 13:02:00', '2015-10-16 13:12:00'],
                ['2015-10-22 06:00:00', '2015-10-22 06:10:00'],
                ['2015-11-01 15:03:00', '2015-11-01 15:13:00'],
                ['2015-11-12 07:14:00', '2015-11-12 07:24:00'],
                ['2015-12-06 23:35:00', '2015-12-06 23:45:00'],
                ['2015-12-08 11:15:00', '2015-12-08 11:25:00'],
                ['2015-12-14 01:12:00', '2015-12-14 01:22:00'],
                ['2016-01-10 09:08:00', '2016-01-10 09:18:00'],
#                ['2016-02-07 20:18:00', '2016-02-07 20:28:00']
    ]
code_run = '1'
#if(code_run):
for trange in trange_list[-1:]:
    print(trange)
    r"""
    RX model from the IDL code.

    Parameters
    ----------
    probe : TYPE, optional
        DESCRIPTION. The default is probe.
    maximum_shear : TYPE, optional
        DESCRIPTION. The default is maximum_shear.
    movie : TYPE, optional.
        The default is movie.
    mmsprobe : TYPE, optional
        DESCRIPTION. The default is mmsprobe.
    times : TYPE, optional
        DESCRIPTION. The default is times.
    trange : TYPE, optional
        Range of times to use. The default is trange.
    model_type : TYPE, optional
        Type of Tsyganenko model. The default is 't96'.
    m_p : float, optional
        Thickness of the magnetopause. The default is 0.5.
    dr : float, optional
        Resolution of the model. The default is 0.25.
    y_min : float, optional
        Minimum y value. The default is -20.
    y_max : float, optional
        Maximum y value. The default is 20.
    z_min : float, optional
        Minimum z value. The default is -20.
    z_max : float, optional
        Maximum z value. The default is 20.
    draw_patch : bool, optional
        Draw the patch of a given radius around the figure. The default is False. If it is set to
        true then the radius of the patch is automatically computed as the maximum value of the
        x-axis.
    save_data : bool, optional
        If True, the data will be saved in an HDF file. The default is False.
    plot_type : String or array of strings, optional
        Type of the plot one wants to get once the code has finished running. The default is
        "shear". Other options are "rx_en", "va_cs", "bisec_msp", "bisec_msh" and "all" for all of
        them.

    Raises
    ------
    KeyError:
        If "plot_type" is not in the list of plot types.
    Returns
    -------
    None.
    """
    probe = None
    omni_level = 'hro'
    maximum_shear = True
    mms_probe = None
    movie = None
    #trange = ['2016-12-24 15:08:00', '2016-12-24 15:12:00']
    #trange = ['2016-12-07 05:11:00', '2016-12-07 05:21:00']
    #trange = ['2016-12-29 03:53:00', '2016-12-29 04:03:00']
    #trange = ['2015-09-08 11:05:00', '2015-09-08 11:15:00']
    model_type = 't96'
    m_p = 0.5  # Magnetopause thichkness
    dr = 0.5  # Resolution of model run in R_E units 
    min_max_val = 15
    y_min = -min_max_val
    y_max = min_max_val
    z_min = -min_max_val
    z_max = min_max_val
    save_data = False
    plot_type = "ttt"
    draw_patch = True


    sw_params = get_sw_params(probe=probe, omni_level=omni_level, trange=trange, mms_probe=mms_probe, verbose=True)

    n_arr_y = int((y_max - y_min) / dr) + 1
    n_arr_z = int((z_max - z_min) / dr) + 1

    bx = np.full((n_arr_y, n_arr_z), np.nan)
    by = np.full((n_arr_y, n_arr_z), np.nan)
    bz = np.full((n_arr_y, n_arr_z), np.nan)

    bx_ext = np.full((n_arr_y, n_arr_z), np.nan)
    by_ext = np.full((n_arr_y, n_arr_z), np.nan)
    bz_ext = np.full((n_arr_y, n_arr_z), np.nan)

    bx_igrf = np.full((n_arr_y, n_arr_z), np.nan)
    by_igrf = np.full((n_arr_y, n_arr_z), np.nan)
    bz_igrf = np.full((n_arr_y, n_arr_z), np.nan)

    b_msx = np.full((n_arr_y, n_arr_z), np.nan)
    b_msy = np.full((n_arr_y, n_arr_z), np.nan)
    b_msz = np.full((n_arr_y, n_arr_z), np.nan)

    x_shu = np.full((n_arr_y, n_arr_z), np.nan)
    y_shu = np.full((n_arr_y, n_arr_z), np.nan)
    z_shu = np.full((n_arr_y, n_arr_z), np.nan)

    rho_sh = np.full((n_arr_y, n_arr_z), np.nan)

    # rp = np.full((n_arr, n_arr), np.nan)

    # r = np.full((n_arr, n_arr), np.nan)
    # zp = np.full((n_arr, n_arr), np.nan)
    # x0 = np.full((n_arr, n_arr), np.nan)

    shear = np.full((n_arr_y, n_arr_z), np.nan)
    rx_en = np.full((n_arr_y, n_arr_z), np.nan)
    gd1 = np.full((n_arr_y, n_arr_z), np.nan)
    gd2 = np.full((n_arr_y, n_arr_z), np.nan)
    va_cs = np.full((n_arr_y, n_arr_z), np.nan)
    bisec_msp = np.full((n_arr_y, n_arr_z), np.nan)
    bisec_msh = np.full((n_arr_y, n_arr_z), np.nan)
    y_coord = np.full((n_arr_y, n_arr_z), np.nan)
    z_coord = np.full((n_arr_y, n_arr_z), np.nan)
    b_sh_ca = np.full((n_arr_y, n_arr_z), np.nan)
    b_sh_mag = np.full((n_arr_y, n_arr_z), np.nan)
    n_sh = np.full((n_arr_y, n_arr_z), np.nan)

    d_theta = np.pi/100

    # Shue et al.,1998, equation 10
    # if (b_imf_z >=0):
    #     ro = (11.44 + 0.013 * b_imf_z) * (p_dyn)**(-1.0/6.6)
    #     #ro = (10.22 + 1.29 * np.tanh(0.184 * (b_imf_z + 8.14))) * (p_dyn)**(-1.0/6.6)
    # else:
    #    ro = (11.44 + 0.14 * b_imf_z) * (p_dyn)**(-1.0/6.6)
    ro = (10.22 + 1.29 * np.tanh(0.184 * (sw_params['b_imf'][2] + 8.14))) * (
                                          sw_params['p_dyn'])**(-1.0/6.6)

    # Shue et al.,1998, equation 11
    # alpha = (0.58 - 0.010 * b_imf_z) * (1 + 0.010 * p_dyn)
    alpha = (0.58 - 0.007 * sw_params['b_imf'][2]) * (1 + 0.024 * np.log(sw_params['p_dyn']))
    rmp = ro * (2/(1 + np.cos(0.0))) ** alpha  # Stand off position of the magnetopause

    A = 2
    len_y = int((y_max - y_min)/dr) + 1
    len_z = int((z_max - z_min)/dr) + 1

    p = mp.Pool()

    input = ((j, k, y_max, z_max, dr, m_p, ro, alpha, rmp, sw_params, 't96') for j in range(len_y) for k in range(len_z))

    print("Running the model\n")
    res = p.map(model_run, input)
    print("Model run complete")

    p.close()
    p.join()

    for r in res:
        j = r[0]
        k = r[1]
        bx[j, k] = r[2]
        by[j, k] = r[3]
        bz[j, k] = r[4]

        shear[j, k] = r[5]
        rx_en[j, k] = r[6]
        va_cs[j, k] = r[7]
        bisec_msp[j, k] = r[8]
        bisec_msh[j, k] = r[9]

    if save_data:
        try:
            fn = f'../data/all_data_rx_model_{dr}re_{m_p}mp_{model_type}_{today_date}.h5'
            data_file = hf.File(fn, 'w')

            data_file.create_dataset('bx', data=bx)
            data_file.create_dataset('by', data=by)
            data_file.create_dataset('bz', data=bz)

            data_file.create_dataset('shear', data=shear)
            data_file.create_dataset('rx_en', data=rx_en)
            data_file.create_dataset('va_cs', data=va_cs)
            data_file.create_dataset('bisec_msp', data=bisec_msp)
            data_file.create_dataset('bisec_msh', data=bisec_msh)

            data_file.close()
            print(f'Date saved to file {fn}')
        except Exception as e:
            print(e)
            print(f'Data not saved to file {fn}. Please make sure that file name is correctly assigned and that the directory exists and you have write permissions')

    # Check if 'plot_type' has length attribute. If it has length attribute then plot the ridge plot
    # for each of the plot type in the list. If it does not have length attribute then plot the
    # ridge plot for the specified plot type.
    ps = sw_params['ps']
    imf_phi = sw_params['imf_clock_angle']
    plot_type = 'all'
    types_of_plot = ['shear', 'rx_en', 'va_cs', 'bisec_msp', 'bisec_msh', 'all']
    if isinstance(plot_type, list):
        for xx in plot_type:
            if xx not in types_of_plot:
                raise ValueError(
                    f'{xx} is not a valid plot type. Please choose from {types_of_plot}'
                    )
        if 'shear' in plot_type:
            print('Plotting shear')
            ridge_finder(image=shear, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
            sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='shear', c_label='Shear',
            c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'rx_en' in plot_type:
            print('Plotting rx_en')
            ridge_finder(image=rx_en/np.nanmax(rx_en), t_range=trange, xrange=[y_min, y_max],
            yrange=[z_min, z_max], sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='rx_en',
            c_label='Rx_en', c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'va_cs' in plot_type:
            print('Plotting va_cs')
            ridge_finder(image=va_cs, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
            sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='va_cs', c_label='Va_cs',
            c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'bisec_msp' in plot_type:
            print('Plotting bisec_msp')
            ridge_finder(image=bisec_msp, t_range=trange, xrange=[y_min, y_max],
            yrange=[z_min, z_max], sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='bisec_msp',
            c_label='Bisec_msp', c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
        if 'bisec_msh' in plot_type:
            print('Plotting bisec_msh')
            ridge_finder(image=bisec_msh, t_range=trange, xrange=[y_min, y_max],
            yrange=[z_min, z_max], sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='bisec_msh',
            c_label='Bisec_msh', c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type == 'all':
        print('Plotting for all')
        ridge_finder(image=shear, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='shear', c_label='Shear',
        c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=rx_en/np.nanmax(rx_en), t_range=trange, xrange=[y_min, y_max],
        yrange=[z_min, z_max], sigma=2.8, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='rx-en_nPa',
        c_label='Reconnection Energy', c_unit='nPa', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=va_cs, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=3., dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='va-cs', c_label='Exhaust Velocity',
        c_unit='km/s', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=bisec_msp, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='bisec_msp', c_label='Bisection Field',
        c_unit='nT', draw_patch=draw_patch, draw_ridge=True)

        ridge_finder(image=bisec_msh, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='bisec_msh', c_label='Bisection Field',c_unit='nT', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='shear':
        print(trange)
        print('Plotting shear')
        ridge_finder(image=shear, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=1, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='shear_test', c_label='Shear',
        c_unit=r'${}^\circ$', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='rx_en':
        print('Plotting rx_en')
        ridge_finder(image=rx_en/np.nanmax(rx_en), t_range=trange, xrange=[y_min, y_max],
        yrange=[z_min, z_max], sigma=2.8, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='rx-en_nPa',
        c_label='Reconnection Energy', c_unit='nPa', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='va_cs':
        print('Plotting va_cs')
        y_val =ridge_finder(image=va_cs, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=3., dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='va-cs', c_label='Exhaust Velocity',
        c_unit='km/s', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='bisec_msp':
        print('Plotting bisec_msp')
        ridge_finder(image=bisec_msp, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='bisec_msp', c_label='Bisection Field',c_unit='nT', draw_patch=draw_patch, draw_ridge=True)
    elif plot_type=='bisec_msh':
        print('Plotting bisec_msh')
        ridge_finder(image=bisec_msh, t_range=trange, xrange=[y_min, y_max], yrange=[z_min, z_max],
        sigma=2.2, dr=dr, dipole_tilt_angle=ps, imf_clock_angle=imf_phi, fig_name='bisec_msh', c_label='Bisection Field',
        c_unit='nT', draw_patch=draw_patch, draw_ridge=True)
    else:
        raise KeyError('plot_type must be one or a list of: all, shear, rx-en, va-cs, bisec_msp, bisec_msh')

    #_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=b_msy, bz=b_msz, save_fig=True,
    #                      scale=40, fig_name="magnetosheath")
    #_ = draping_field_plot(x_coord=y_shu, y_coord=z_shu, by=by, bz=bz, save_fig=True, scale=120,
    #                 fig_name="magnetosphere")
print(f'Took {round(time.time() - start, 3)} seconds')