import pyspedas as spd
import pytplot as ptt

trange = ['2017-09-02 17:30:30', '2017-09-02 17:31:30']

mms_fpi_varnames = ['mms3_dis_bulkv_gse_brst']
data_rate = 'brst'
mms_probe_num = 3
time_clip = True

_ = spd.mms.fpi(trange=trange, probe=mms_probe_num, data_rate=data_rate, level='l2',
                datatype='dis-moms', varnames=mms_fpi_varnames, time_clip=time_clip,
                latest_version=True)

fpi_time_unix = ptt.get_data(mms_fpi_varnames[0])[0]
fpi_v_gse = ptt.get_data(mms_fpi_varnames[0])[1:][0]

fpi_time_utc = spd.time_datetime(fpi_time_unix)

# Print the first and last time in the data
print(f"\033[1;32m Start time: {fpi_time_utc[0].strftime('%Y-%m-%d %H:%M:%S')} \033[0m")
print(f"\033[1;32m End time: {fpi_time_utc[-1].strftime('%Y-%m-%d %H:%M:%S')} \033[0m")