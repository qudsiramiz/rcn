from geopack.t01 import t01
from geopack.t96 import t96
from geopack import geopack
import datetime

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

time = 1420070400
ps = geopack.recalc(time)
par = [5,-100,1,1,0,0,0,0,0,0]

input_r = [[5, 0.5, 0.5]]

v = 61

x_gsm = np.linspace(-15.1, 15, v)
y_gsm = np.linspace(-15.1, 15, v)
z_gsm = np.linspace(-15.1, 15, v)

bx_t96 = np.full((v,v,v), np.nan)
by_t96 = np.full((v,v,v), np.nan)
bz_t96 = np.full((v,v,v), np.nan)

bx_t01 = np.full((v,v,v), np.nan)
by_t01 = np.full((v,v,v), np.nan)
bz_t01 = np.full((v,v,v), np.nan)

bx_igrf = np.full((v,v, v), np.nan)
by_igrf = np.full((v,v, v), np.nan)
bz_igrf = np.full((v,v, v), np.nan)

bx_t96_total = np.full((v,v,v), np.nan)
by_t96_total = np.full((v,v,v), np.nan)
bz_t96_total = np.full((v,v,v), np.nan)

bx_t01_total = np.full((v,v,v), np.nan)
by_t01_total = np.full((v,v,v), np.nan)
bz_t01_total = np.full((v,v,v), np.nan)


for i in range(len(x_gsm)):
    print(i)
    for j in range(len(y_gsm)):
        for k in range(len(z_gsm)):

            bx_t96[i,j,k], by_t96[i,j,k], bz_t96[i,j,k] = t96(par, ps, x_gsm[i], y_gsm[j], z_gsm[k])
            bx_t01[i,j,k], by_t01[i,j,k], bz_t01[i,j,k] = t01(par, ps, x_gsm[i], y_gsm[j], z_gsm[k])

            bx_igrf[i,j,k], by_igrf[i,j,k], bz_igrf[i,j,k] = geopack.igrf_gsm(x_gsm[i], y_gsm[j],
                                                                              z_gsm[k])

            bx_t96_total[i, j, k] = bx_igrf[i, j, k] + bx_t96[i, j, k]
            by_t96_total[i, j, k] = by_igrf[i, j, k] + by_t96[i, j, k]
            bz_t96_total[i, j, k] = bz_igrf[i, j, k] + bz_t96[i, j, k]

            bx_t01_total[i, j, k] = bx_igrf[i, j, k] + bx_t01[i, j, k]
            by_t01_total[i, j, k] = by_igrf[i, j, k] + by_t01[i, j, k]
            bz_t01_total[i, j, k] = bz_igrf[i, j, k] + bz_t01[i, j, k]

data_file = hf.File(f'../data/total_fields_ipy_0.5_15_15_{par}_{today_date}.h5')
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
'''
for (xgsm, ygsm, zgsm) in zp(x_gsm, y_gsm, z_gsm):

    bx_t96, by_t96, bz_t96 = t96(par, ps, xgsm, ygsm, zgsm)
    bx_t01, by_t01, bz_t01 = t01(par, ps, xgsm, ygsm, zgsm)

    bx_igrf, by_igrf, bz_igrf = geopack.igrf_gsm(xgsm, ygsm, zgsm)

    bx = bx_igrf + bx_t96
    by = by_igrf + by_t96
    bz = bz_igrf + bz_t96

    bx01 = bx_igrf + bx_t01
    by01 = by_igrf + by_t01
    bz01 = bz_igrf + bz_t01

    print('t96-->', round(bx_t96,4), round(by_t96,4), round(bz_t96,4), round(ps, 4))
    print('t01-->', round(bx_t01,4), round(by_t01,4), round(bz_t01,4))
    print('igrf-->', round(bx_igrf,4), round(by_igrf,4), round(bz_igrf,4))
    print('total (t96)-->', round(bx,4), round(by,4), round(bz,4))
    print('total (t01)-->', round(bx01,4), round(by01,4), round(bz01,4))
'''
