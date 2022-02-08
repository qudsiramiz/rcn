import numpy as np
import scipy as sc
import scipy.optimize
from matplotlib import pyplot as plt

def f(x, y):
    """ Function of the surface"""
    # example equation
    z = x**2 + y**2 -10
    return z

p0 = np.array([1, 2, 1]) # starting point for the line
direction = np.array( [1, 3, -1]) # direction vector

def line_func(t):
    """Function of the straight line.
    :param t:     curve-parameter of the line

    :returns      xyz-value as array"""
    return p0 + t*direction

def target_func(t):
    """Function that will be minimized by fmin
    :param t:      curve parameter of the straight line

    :returns:      (z_line(t) - z_surface(t))**2 – this is zero
                   at intersection points"""
    p_line = line_func(t)
    z_surface = f(*p_line[:2])
    return np.sum((p_line[2] - z_surface)**2)


def line_fnc(
    r0=np.array([0, 0, 0]),
    b_msh=np.array([-5, 0, 0]),
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
    line_intrp2 = interpolate.CubicSpline(x, y)
    return line_intrp2

def target_fnc(r, b_msh, r0, x, y):
    p_line = line_fnc(r=r, b_msh=b_msh, r0=[0, 0, 0])
    line_intrp = line_fnc_der(x, y)
    z_surface = line_intrp(p_line[2])

    return np.sum((p_line[2] - z_surface)**2)

r_opt = sc.optimize.minimize(target_fnc, 0, np.array([-5, -0.14, -3.5]))
r_intsc = line_fnc(r_opt.x)
#t_opt = sc.optimize.fmin(target_func, x0=-10)
#intersection_point = line_func(t_opt)