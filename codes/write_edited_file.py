import h5py as hf
 
#dat_idl = hf.File('../data/total_fields_idl_0.501_15_15_[5,0,1,1,0,0,0,0,0,0]_2021-10-14.h5')
dat_idl = hf.File('../data/total_fields_ipy_0.5_15_15_[5, -100, 1, 1, 0, 0, 0, 0, 0, 0]_2021-10-19_v01.h5')

dat_new = {}
for key in list(dat_idl.keys()):
    temp_data = dat_idl[key][:,:]
    temp_data[temp_data<-1000] = np.nan
    temp_data[temp_data>1000] = np.nan
    dat_new[key] = temp_data

fn = '../data/total_fields_ipy_0.5_15_15_[5, -100, 1, 1, 0, 0, 0, 0, 0, 0]_2021-10-19_edited.h5'

new_dat = hf.File(fn, 'w')
for key in dat_new.keys():
    new_dat.create_dataset(key, data=dat_new[key])
new_dat.close()
'''
bx_igrf_idl = dat_idl['bx_igrf'][:,:]
bx_igrf_idl[bx_igrf_idl<-1000] = np.nan

by_igrf_idl = dat_idl['by_igrf'][:,:]
by_igrf_idl[by_igrf_idl<-1000] = np.nan

bz_igrf_idl = dat_idl['bz_igrf'][:,:]
bz_igrf_idl[bz_igrf_idl<-1000] = np.nan

bx_ext_idl = dat_idl['bx_ext'][:,:]
bx_ext_idl[bx_ext_idl<-1000] = np.nan

by_ext_idl = dat_idl['by_ext'][:,:]
by_ext_idl[by_ext_idl<-1000] = np.nan

bz_ext_idl = dat_idl['bz_ext'][:,:]
bz_ext_idl[bz_ext_idl<-1000] = np.nan

bx_idl = dat_idl['bx'][:,:]
bx_idl[bx_idl<-1000] = np.nan

by_idl = dat_idl['by'][:,:]
by_idl[by_idl<-1000] = np.nan

bz_idl = dat_idl['bz'][:,:]
bz_idl[bz_idl<-1000] = np.nan

xs_idl = dat_idl['x_shu'][:,:]
xs_idl[xs_idl<-2000] = np.nan

ys_idl = dat_idl['y_shu'][:,:]
ys_idl[ys_idl<-1000] = np.nan

zs_idl = dat_idl['z_shu'][:,:]
zs_idl[zs_idl<-1000] = np.nan

sr_idl = dat_idl['shear'][:,:]
sr_idl[sr_idl<-1000] = np.nan
sr_idl[sr_idl>200] = np.nan


fn = 'param_total_edited_shear_v07.h5'
new_dat = hf.File(fn, 'w')
new_dat.create_dataset('bx', data=np.transpose(bx_idl))
new_dat.create_dataset('by', data=np.transpose(by_idl))
new_dat.create_dataset('bz', data=np.transpose(bz_idl))

new_dat.create_dataset('bx_igrf', data=np.transpose(bx_igrf_idl))
new_dat.create_dataset('by_igrf', data=np.transpose(by_igrf_idl))
new_dat.create_dataset('bz_igrf', data=np.transpose(bz_igrf_idl))

new_dat.create_dataset('bx_ext', data=np.transpose(bx_ext_idl))
new_dat.create_dataset('by_ext', data=np.transpose(by_ext_idl))
new_dat.create_dataset('bz_ext', data=np.transpose(bz_ext_idl))

new_dat.create_dataset('x_shu', data=np.transpose(xs_idl))
new_dat.create_dataset('y_shu', data=np.transpose(ys_idl))
new_dat.create_dataset('z_shu', data=np.transpose(zs_idl))
new_dat.create_dataset('shear', data=np.transpose(sr_idl))

new_dat.close()
'''
