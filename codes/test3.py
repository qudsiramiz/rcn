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

def target_func(t, a):
    """Function that will be minimized by fmin
    :param t:      curve parameter of the straight line

    :returns:      (z_line(t) - z_surface(t))**2 – this is zero
                   at intersection points"""
    p_line = line_func(t)
    z_surface = f(*p_line[:2]) + a
    return np.sum((p_line[2] - z_surface)**2)

t_opt = sc.optimize.fmin(target_func, x0=-10, args=0)
intersection_point = line_func(t_opt)