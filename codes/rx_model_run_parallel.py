import geopack.geopack as gp
import PyGeopack as gpp
#from geopack import geopack
#from geopack.geopack import igrf_gsm
#from geopack.t96 import t96
#from geopack.t01 import t01
#from geopack.t04 import t04
import datetime
from dateutil import parser
import time
#from func_t01 import func_t01
now = time.time()

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

v = 61
t1 = datetime.datetime(2015, 1, 1, 0, 0, 0)
t0 = datetime.datetime(1970, 1, 1)
time_o = (t1 - t0).total_seconds()
ps = gp.recalc(time_o)
par = [5, -100, 1, 1, 0, 0, 0, 0, 0, 0]

x_gsm = np.linspace(-15.1, 15, v)
y_gsm = np.linspace(-15.1, 15, v)
z_gsm = np.linspace(-15.1, 15, v)

def func_all(*args) :

    i = args[0][0]
    j = args[0][1]
    k = args[0][2]

    bx_t96, by_t96, bz_t96 = gp.t96.t96(par, ps, x_gsm[i], y_gsm[j], z_gsm[k])
    bx_t01, by_t01, bz_t01 = gp.t01.t01(par, ps, x_gsm[i], y_gsm[j], z_gsm[k])
    bx_igrf, by_igrf, bz_igrf = gp.igrf_gsm(x_gsm[i], y_gsm[j], z_gsm[k])

    Bx_t96, By_t96, Bz_t96 = gpp.ModelField(x_gsm, y_gsm, z_gsm, 20150101, 0, Model='T96',
                                            CoordIn='GSM', CoordOut='GSM', parmod=par)
    Bx_t01, By_t01, Bz_t01 = gpp.ModelField(x_gsm, y_gsm, z_gsm, 20150101, 0, Model='T01',
                                            CoordIn='GSM', CoordOut='GSM', parmod=par)

    return i, j, k, bx_t96, by_t96, bz_t96, bx_t01, by_t01, bz_t01, bx_igrf, by_igrf, bz_igrf, Bx_t96[0], By_t96[0], Bz_t96[0], Bx_t01[0], By_t01[0], Bz_t01[0]

import multiprocessing as mp

p = mp.Pool()

dim1, dim2, dim3 = v, v, v

import itertools

#input = ((i, j, k) for i,j,k in itertools.combinations_with_replacement(range(dim3), 3))
input = ((i, j, k) for i,j,k in itertools.product(range(dim3), repeat=3))

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

Bx_t96 = np.empty((dim1, dim2, dim3))
By_t96 = np.empty((dim1, dim2, dim3))
Bz_t96 = np.empty((dim1, dim2, dim3))

Bx_t01 = np.empty((dim1, dim2, dim3))
By_t01 = np.empty((dim1, dim2, dim3))
Bz_t01 = np.empty((dim1, dim2, dim3))

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

    Bx_t96[i, j, k] = r[12]
    By_t96[i, j, k] = r[13]
    Bz_t96[i, j, k] = r[14]

    Bx_t01[i, j, k] = r[15]
    By_t01[i, j, k] = r[16]
    Bz_t01[i, j, k] = r[17]

bx_t96_total = bx_igrf + bx_t96
by_t96_total = by_igrf + by_t96
bz_t96_total = bz_igrf + bz_t96

bx_t01_total = bx_igrf + bx_t01
by_t01_total = by_igrf + by_t01
bz_t01_total = bz_igrf + bz_t01

fn = f'../data/total_fields_ipy_{v}_15_15_15_{par}_{today_date}_geo.h5'

data_file = hf.File(fn, 'w')
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

data_file.create_dataset('Bx_t96', data=Bx_t96)
data_file.create_dataset('By_t96', data=By_t96)
data_file.create_dataset('Bz_t96', data=Bz_t96)

data_file.create_dataset('Bx_t01', data=Bx_t01)
data_file.create_dataset('By_t01', data=By_t01)
data_file.create_dataset('Bz_t01', data=Bz_t01)

data_file.close()

print(f'took {round(time.time() - now, 3)} seconds')
