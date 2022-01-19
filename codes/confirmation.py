import geopack.geopack as gp
import PyGeopack as gpp
import datetime
import time
from tabulate import tabulate
#from func_t01 import func_t01
now = time.time()

today_date = datetime.datetime.today().strftime('%Y-%m-%d')

t1 = datetime.datetime(2001, 1, 1, 0, 0, 0)
t0 = datetime.datetime(1970, 1, 1)
time_o = (t1 - t0).total_seconds()
ps = gp.recalc(time_o)
par = [5, 0, 1, 1, 0, 0, 0, 0, 0, 0]

x_gsm = -8.00
y_gsm = -2.00
z_gsm = -3.00


bx_t96, by_t96, bz_t96 = gp.t96.t96(par, ps, x_gsm, y_gsm, z_gsm)
bx_t01, by_t01, bz_t01 = gp.t01.t01(par, ps, x_gsm, y_gsm, z_gsm)
bx_igrf, by_igrf, bz_igrf = gp.igrf_gsm(x_gsm, y_gsm, z_gsm)

bx_total_t96 = bx_t96 + bx_igrf
by_total_t96 = by_t96 + by_igrf
bz_total_t96 = bz_t96 + bz_igrf

bx_total_t01 = bx_t01 + bx_igrf
by_total_t01 = by_t01 + by_igrf
bz_total_t01 = bz_t01 + bz_igrf

#b_t96_idl = [22.705237, -2.5203028, -52.475987]

#diff_t96 = [b_t96_idl[0] - bx_t96, b_t96_idl[1] - by_t96, b_t96_idl[2] - bz_t96]

#print('Geo T96==>    ', bx_igrf + bx_t96, by_igrf + by_t96, bz_igrf + bz_t96)

Bx_t96, By_t96, Bz_t96 = gpp.ModelField(x_gsm, y_gsm, z_gsm, 20010101, 0, Model='T96',
                                        CoordIn='GSM', CoordOut='GSM', parmod=par)

#print('PyGeo T96 ==> ', Bx_t96[0], By_t96[0], Bz_t96[0])

Bx_t01, By_t01, Bz_t01 = gpp.ModelField(x_gsm, y_gsm, z_gsm, 20010101, 0, Model='T01',
                                        CoordIn='GSM', CoordOut='GSM', parmod=par)

#print('Geo T01 ==>   ', bx_igrf + bx_t01, by_igrf + by_t01, bz_igrf + bz_t01)
#print('PyGeo T01 ==> ', Bx_t01[0], By_t01[0], Bz_t01[0])

print(tabulate([['Position (in GSM)', x_gsm, y_gsm, z_gsm]], tablefmt='simple', numalign='center'))
print(tabulate([['Geo T96', bx_total_t96, by_total_t96, bz_total_t96],
                ['PyGeo T96', Bx_t96, By_t96, Bz_t96],
                ['Geo T01', bx_total_t01, by_total_t01, bz_total_t01],
                ['PyGeo T01', Bx_t01, By_t01, Bz_t01]],
                headers=['Model', 'B_x', 'B_y', 'B_z'], tablefmt='grid', numalign='center',
                stralign='center'))
