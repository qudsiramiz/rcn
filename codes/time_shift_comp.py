import matplotlib.pyplot as plt
import numpy as np
import pyspedas as spd
import pytplot as ptt
import pandas as pd

time = ['2016-11-23 07:50:30', '2016-11-23 07:55:30']

download_data = True
if(download_data):
    omni_varnames = ['BX_GSE', 'BY_GSM', 'BZ_GSM', 'proton_density', 'Vx', 'Vy', 'Vz', 'SYM_H',
                    'Timeshift', 'RMS_Timeshift']


    omni_vars = spd.omni.data(trange=time, varnames=omni_varnames, level='hro', time_clip=True)

    omni_time = ptt.get_data(omni_vars[0])[0]
    omni_bx_gse = ptt.get_data(omni_vars[0])[1]
    omni_by_gsm = ptt.get_data(omni_vars[1])[1]
    omni_bz_gsm = ptt.get_data(omni_vars[2])[1]
    omni_np = ptt.get_data(omni_vars[3])[1]
    omni_vx = ptt.get_data(omni_vars[4])[1]
    omni_vy = ptt.get_data(omni_vars[5])[1]
    omni_vz = ptt.get_data(omni_vars[6])[1]
    omni_sym_h = ptt.get_data(omni_vars[7])[1]
    omni_time_shift = ptt.get_data(omni_vars[8])[1]
    omni_rms_time_shift = ptt.get_data(omni_vars[9])[1]

    mms_fgm_varnames = [f'mms1_fgm_b_gsm_srvy_l2_bvec']

    mms_vars = spd.mms.fgm(trange=time, probe=1, data_rate='srvy', level='l2', time_clip=True,
                           latest_version=True)

    mms_time = ptt.get_data(mms_fgm_varnames[0])[0]

    mms_fgm_b_gsm = ptt.get_data(mms_fgm_varnames[0])[1:4][0]

# Create a dataframe from MMS data
mms_df = pd.DataFrame(data={'time': mms_time, 'b_x_gsm': mms_fgm_b_gsm[:, 0],
                     'b_y_gsm': mms_fgm_b_gsm[:, 1], 'b_z_gsm': mms_fgm_b_gsm[:, 2]})

# Compute 1 minute average of MMS data
bx_gsm_average = mms_df.b_x_gsm.rolling(window=16).mean()
by_gsm_average = mms_df.b_y_gsm.rolling(window=16).mean()
bz_gsm_average = mms_df.b_z_gsm.rolling(window=16).mean()

bx_means = []
by_means = []
bz_means = []

for count in range(len(omni_time)):

    bx_means.append(np.nanmean(mms_df.b_x_gsm[(count*960):(count + 1)*960]))
    by_means.append(np.nanmean(mms_df.b_y_gsm[(count*960):(count + 1)*960]))
    bz_means.append(np.nanmean(mms_df.b_z_gsm[(count*960):(count + 1)*960]))

#mms_time_edit = np.linspace(omni_time[0], omni_time[-1], len(mms_time))
plt.close('all')
plt.figure(figsize=(10, 10))
plt.plot(mms_time, mms_fgm_b_gsm[:,2], color='gray', marker='.', label='MMS1_X', alpha=0.05)
plt.plot(omni_time, bz_means, 'k*', label='MMS1_X_average', ms=5)
plt.plot(omni_time, omni_bz_gsm, 'b.', label='OMNI_X', ms=5, alpha=0.5)
plt.legend()
#plt.show()
# save figure as pdf
plt.savefig('mms_bz_gsm.png', format='png', dpi=250)