'''
import numpy as np
import scipy as sp
import scipy as sc
import scipy.optimize
from matplotlib import pyplot as plt

def f(x, y):
    """ Function of the surface"""
    # example equation

    z = x**2 + y**2 -10
    return z

def line_func(p0, t, b_msh):
    """Function of the straight line.
    :param t:     curve-parameter of the line

    :returns      xyz-value as array"""
    b_msh_norm = np.linalg.norm(b_msh)
    b_msh_dir = b_msh / b_msh_norm

    return p0 + t*b_msh_dir

def target_func(p0, t, b_msh, line_intrp):
    """Function that will be minimized by fmin
    :param t:      curve parameter of the straight line

    :returns:      (z_line(t) - z_surface(t))**2 – this is zero
                   at intersection points"""
    p_line = line_func(p0, t, b_msh)
    z_surface = line_intrp(p_line[0])
    print(z_surface)
    return np.sum((p_line[1] - z_surface)**2)


def line_fnc(
    r0=np.array([0, 0]),
    b_msh=np.array([0, 0]),
    r = 0,
    ):
    """
    Function to compute the equation of a line along the direction of the magnetic field.

    Parameters
    ----------
    r0 : array
        Starting position of the line
    b_msh : array, optional
        Magnetic field from magnetosheath at the magnetopause
    r : float, optional
        Position of the point where the line is computed
    """
    # Compute the direction of the magnetic field line
    b_msh_norm = np.linalg.norm(b_msh)
    b_msh_dir = b_msh / b_msh_norm

    # return the line function
    return r0 + r * b_msh_dir

def line_fnc_der(x, y):
    line_intrp = sp.interpolate.CubicSpline(x, y)
    return line_intrp
#np.random.seed(10)
x = np.linspace(-5, 5, 100)
y = (0.5 - np.random.rand(100)) * 1000
line_intrp = line_fnc_der(x, y)

def target_fnc(r, r0, b_msh, line_intrp):
    p_line = line_fnc(r0=r0, b_msh=b_msh, r=r)
    #print(p_line)
    z_surface = line_intrp(p_line[0])
    #print(z_surface)
    return np.sum((p_line[1] - z_surface)**2)

#r = 2
#r0 = np.array([-5, -0.14])
#b_msh = np.array([-5, -3])
#p0 = r0
#r_opt = sp.optimize.fmin(target_fnc, x0=r, args=(r0, b_msh, line_intrp))
#r_opt2 = sp.optimize.minimize(target_fnc, x0=r, args=(r0, b_msh, line_intrp))
#print(r_opt, r_opt2.x)
#r_intsc = line_fnc(r_opt)
#t_opt = sc.optimize.fmin(target_func, x0=-10)
#intersection_point = line_func(t_opt)


trange = ['2015-10-16 10:28:30', '2015-10-16 10:38:30']
time_clip = True

if (mms_probe_num is not None):
    mms_fgm_varnames = [f'mms{mms_probe_num}_fgm_b_gsm_srvy_l2_bvec']
    mms_fgm_vars = spd.mms.fgm(trange=trange, probe=mms_probe_num,
                               time_clip=time_clip, latest_version=True)
    print(mms_fgm_vars, 'foo \n')
    mms_fgm_time = ptt.get_data(mms_fgm_varnames[0])[0]
    mms_fgm_b_gsm = ptt.get_data(mms_fgm_varnames[0])[1:4][0]
'''
'''
import importlib
import rx_model_funcs as rmf
importlib.reload(rmf)

figure_inputs = {
    "image" : [shear, rx_en/np.nanmax(rx_en), va_cs, bisec_msp],
    "convolution_order" : [0, 1, 1, 1],
    "t_range" : trange,
    "b_imf" : np.round(sw_params['b_imf'],2),
    "b_msh" : np.round(sw_params['mms_b_gsm'],2),
    "xrange" : [y_min, y_max],
    "yrange" : [z_min, z_max],
    "mms_probe_num" : mms_probe_num,
    "mms_sc_pos" : np.round(np.nanmean(sw_params['mms_sc_pos'], axis=0), 2),
    "dr" : dr,
    "dipole_tilt_angle" : sw_params['ps'],
    "imf_clock_angle" : sw_params['imf_clock_angle'],
    "sigma" : [2, 2, 2, 2],
    "mode" : "nearest",
    "alpha" : 1,
    "vmin" : [0, 0, None, None],
    "vmax" : [180, 1, None, None],
    "cmap_list" : ["viridis", "cividis", "plasma", "magma"],
    "draw_patch" : [True, True, True, True],
    "draw_ridge" : [True, True, True, True],
    "save_fig" : True,
    "fig_name" : f'all_ridge_plots_v2',
    #"fig_format" : 'png',
    "c_label" : ['Shear', 'Reconnection Energy', 'Exhaust Velocity', 'Bisection Field'],
    "c_unit" : [r'${}^\\circ$', 'nPa', 'km/s', 'nT'],
    "wspace" : 0.0,
    "hspace" : 0.17,
    "fig_size" : (8.775, 10),
    "box_style": dict(boxstyle='round', facecolor='black', alpha=0.8),
    "title_y_pos" : 1.09,
    "interpolation" : 'gaussian',
    "tsy_model" : model_type,
    "dark_mode" : True,
    "rc_file_name" : f"reconnection_line_data_mms{mms_probe_num}.csv",
    "rc_folder" : "../data/rx_d/",
    "save_rc_file" : False
}

y, xx, yy = rmf.ridge_finder_multiple(**figure_inputs, fig_format='png')
'''
