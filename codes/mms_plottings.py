"""Code for plotting MMS data using PySpedas."""

from datetime import datetime
import pytz
import numpy as np

from pyspedas import time_string, time_double
from pyspedas.mms import fpi, fgm

from pytplot import tplot, options, get_data, store_data

'''
    tplot: plots the data
    options: for plotting options, like colormap, line properties, legend names etc etc.
    get_data: Get the data corresponding to a specific variable, returns a tuple
    time_string: converts time from unix format to a string
    store_data: stores data that can be later used with tplot to plot 'store_dat	('var_name',
    data={'x':var_times, 'y':var_values})', 'tplot('var_name')' will then plot the data
'''

# Load the electron moments data

# electron_vars = fpi(datatype='des-moms', trange=['2015-10-16', '2015-10-17'],
# center_measurement=True)

# Plot the electron data

# tplot(['mms1_des_energyspectr_omni_fast', 'mms1_des_bulkv_gse_fast',
# 'mms1_des_numberdensity_fast'])

# Load the ion moments data
start_dates = ['2015-09-19 07:40:30', '2015-10-16 10:30:30', '2015-10-16 13:04:02',
               '2015-10-22 06-02-22', '2015-11-01 15:05:06']
end_dates = ['2015-09-19 07:46:30', '2015-10-16 10:36:30', '2015-10-16 13:10:02',
             '2015-10-22 06-08-22', '2015-11-01 15:11:06']

for i in range(len(start_dates)):
    ion_vars = fpi(datatype='dis-moms', trange=[start_dates[i], end_dates[i]],
                   center_measurement=True, time_clip=True)

    mms_fgm = fgm(trange=[start_dates[i], end_dates[i]], data_rate='brst', time_clip=True)

    options('mms1_dis_bulkv_gse_fast', 'legend_names', ['vx_gse', 'vy_gse', 'vz_gse'])
    options('mms1_dis_bulkv_gse_fast', 'color', ['r', 'b', 'g'])

    # Plot the ion data
    tplot(['mms1_fgm_b_gsm_brst_l2', 'mms1_dis_energyspectr_omni_fast',
           'mms1_dis_numberdensity_fast', 'mms1_dis_temppara_fast', 'mms1_dis_bulkv_gse_fast'],
          qt=True, combine_axes=True,
          save_png=f'/home/ahmadr/SpacePlasFest/mms/rcn/figures/{start_dates[i]}.png')

    t1_v, m1_v = get_data('mms1_dis_bulkv_gse_fast')

    time1_array = np.array([datetime.fromtimestamp(xx, pytz.timezone("UTC")) for xx in t1_v])

# Load the burst mode data for both electrons and ions

# both_vars = fpi(data_rate='brst', trange=['2015-10-16/13:06', '2015-10-16/13:07'],
# center_measurement=True)

# Plot the burst mode data

# tplot(['mms1_des_bulkv_gse_brst', 'mms1_dis_bulkv_gse_brst', 'mms1_des_numberdensity_brst',
# 'mms1_dis_numberdensity_brst'])

# Find all the loaded variables

# from pytplot import tplot_names
# tplot_names()
