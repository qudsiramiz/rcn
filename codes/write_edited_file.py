import h5py as hf
import numpy as np
'''
The file removes the fill values (-9999 from IDL) and replaces them with
np.nan. It also rotates the matrix by taking its transpose, since the way data
is saved in python and IDL are transpose of each other.
'''

dat_idl = hf.File('idl_data.h5')

dat_new = {}
for key in list(dat_idl.keys()):
    temp_data = dat_idl[key][:,:]
    temp_data[temp_data<-1000] = np.nan
    temp_data[temp_data>1000] = np.nan
    dat_new[key] = temp_data

fn = 'idl_data_edited.h5'

new_dat = hf.File(fn, 'w')
for key in dat_new.keys():
    new_dat.create_dataset(key, data=dat_new[key])
new_dat.close()
