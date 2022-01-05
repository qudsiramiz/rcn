# Import the required packages
import geopack
from geopack import geopack
from geopack.geopack import igrf_gsm
from geopack.t96 import t96
from geopack.t01 import t01
import time
import numpy as np
import multiprocessing as mp
import itertools
import h5py as hf


# from func_t01 import func_t01

now = time.time()

# Define the length of array
v = 61

# Define the observation time
time_o = 1420070400

# Compute dipole tilt
ps = geopack.recalc(time_o)

# Define the param
param = [5, -100, 1, 1, 0, 0, 0, 0, 0, 0]

x_gsm = np.linspace(-15.1, 15, v)
y_gsm = np.linspace(-15.1, 15, v)
z_gsm = np.linspace(-15.1, 15, v)

# Define the function to compute magnetic fields from different models for different values of
# x_gsm, y_gsm and z_gsm


def func_all(*args):

    i = args[0][0]
    j = args[0][1]
    k = args[0][2]

    bx_t96, by_t96, bz_t96 = t96(param, ps, x_gsm[i], y_gsm[j], z_gsm[k])
    bx_t01, by_t01, bz_t01 = t01(param, ps, x_gsm[i], y_gsm[j], z_gsm[k])
    bx_igrf, by_igrf, bz_igrf = igrf_gsm(x_gsm[i], y_gsm[j], z_gsm[k])

    return i, j, k, bx_t96, by_t96, bz_t96, bx_t01, by_t01, bz_t01, bx_igrf, by_igrf, bz_igrf


p = mp.Pool()

# Define the dimension of input/output
dim1, dim2, dim3 = v, v, v

input = ((i, j, k) for i, j, k in itertools.product(range(dim3), repeat=3))

res = p.map(func_all, input)

p.close()
p.join()

bx_t96 = np.empty((dim1, dim2, dim3))
by_t96 = np.empty((dim1, dim2, dim3))
bz_t96 = np.empty((dim1, dim2, dim3))

bx_t01 = np.empty((dim1, dim2, dim3))
by_t01 = np.empty((dim1, dim2, dim3))
bz_t01 = np.empty((dim1, dim2, dim3))

bx_igrf = np.empty((dim1, dim2, dim3))
by_igrf = np.empty((dim1, dim2, dim3))
bz_igrf = np.empty((dim1, dim2, dim3))

for r in res:
    i = r[0]
    j = r[1]
    k = r[2]

    bx_t96[i, j, k] = r[3]
    by_t96[i, j, k] = r[4]
    bz_t96[i, j, k] = r[5]

    bx_t01[i, j, k] = r[6]
    by_t01[i, j, k] = r[7]
    bz_t01[i, j, k] = r[8]

    bx_igrf[i, j, k] = r[9]
    by_igrf[i, j, k] = r[10]
    bz_igrf[i, j, k] = r[11]

bx_t96_total = bx_igrf + bx_t96
by_t96_total = by_igrf + by_t96
bz_t96_total = bz_igrf + bz_t96

bx_t01_total = bx_igrf + bx_t01
by_t01_total = by_igrf + by_t01
bz_t01_total = bz_igrf + bz_t01

data_file = hf.File('python_data.h5', 'w')
data_file.create_dataset('bx_t96', data=bx_t96)
data_file.create_dataset('by_t96', data=by_t96)
data_file.create_dataset('bz_t96', data=bz_t96)

data_file.create_dataset('bx_t01', data=bx_t01)
data_file.create_dataset('by_t01', data=by_t01)
data_file.create_dataset('bz_t01', data=bz_t01)

data_file.create_dataset('bx_igrf', data=bx_igrf)
data_file.create_dataset('by_igrf', data=by_igrf)
data_file.create_dataset('bz_igrf', data=bz_igrf)

data_file.create_dataset('bx_t96_total', data=bx_t96_total)
data_file.create_dataset('by_t96_total', data=by_t96_total)
data_file.create_dataset('bz_t96_total', data=bz_t96_total)

data_file.create_dataset('bx_t01_total', data=bx_t01_total)
data_file.create_dataset('by_t01_total', data=by_t01_total)
data_file.create_dataset('bz_t01_total', data=bz_t01_total)

data_file.close()

print(f'took {round(time.time() - now, 3)} seconds')
