import time
import numpy as np
import pandas as pd
import pyspedas as spd
import pytplot as ptt
import datetime

start_time = time.time()
# Read the list of dates from the csv file
df = pd.read_csv("../data/mms_magnetopause_crossings.csv")

r_e = 6378.137  # Earth radius in km

mms_varnames = ['mms3_mec_r_gsm']
dt = 15

with open('../data/mms_magnetopause_crossings_limited_v2.csv', 'w') as f:
    f.write("time, x, y, z, r_yz, index\n")
    # Go through each date and find the position of MMS at those dates
    for i, time_crossing in enumerate(df.DateStart[:]):

        trange = [time_crossing]
        # If trange has only one element, then make it a list of length 2 with 'dt' minutes of 
        # padding
        if len(trange) == 1:
            trange_date = datetime.datetime.strptime(
                trange[0], '%Y-%m-%d %H:%M:%S.%f')
            trange_date_min = trange_date - datetime.timedelta(seconds=dt)
            trange_date_max = trange_date + datetime.timedelta(seconds=dt)
            trange = [trange_date_min.strftime('%Y-%m-%d %H:%M:%S'),
                      trange_date_max.strftime('%Y-%m-%d %H:%M:%S')]
        mms_vars = spd.mms.mec(trange=trange, varnames=mms_varnames, probe=3,
                               data_rate='srvy', level='l2', time_clip=True,
                               latest_version=True)
        mms_time = ptt.get_data(mms_vars[0])[0]
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        mms_sc_pos = ptt.get_data(mms_vars[0])[1:3][0][0] / r_e
        x = mms_sc_pos[0]
        y = mms_sc_pos[1]
        z = mms_sc_pos[2]
        r_yz = np.sqrt(y**2 + z**2)
        # Open a file to save the data to a csv file

        # Check if the crossing was close to the sub-solar point
        if x > -5 and x < 12 and r_yz < 12:
            # print("Crossing at: " + time_crossing)
            # print("Position of MMS in GSM coordinates: " +
            #       str(mms_sc_pos))
            # print("Distance from sub-solar point: " + str(r_yz))
            # print("\n")
            # Write the time of crossing and the position of MMS in GSM coordinates to the csv file
            f.write(time_crossing + "," + str(x) + "," + str(y) +
                    "," + str(z) + "," + str(r_yz) + "," + str(i) + "\n")
    f.close()
# print(mms_sc_pos)

# Steps
# 1. Get the time of the crossing
# 2. Get the position of MMS at that time
# 3. Check if mms_spc_y and mms_spc_z are within the range of the crossing (
# sqrt(mms_spc_y**2 + mms_spc_z**2) <= 10 )
# 4. Find the value of x_shu 
# 5. Check if mms_sc_x is within the range of x_shu +/- 1 r_e
# 6. If yes, then check the Walen relation
# 7. If the Walen relation is satisfied, then check if there was a jet reversal within 5 minutes of
#    magnetopause crossing
# 8. If yes, then add the time of the crossing to the list of magnetopause crossings with
#    reconnection

# Print amount of time taken to run the script
print("Time taken to run the script: " + str(time.time() - start_time) + " seconds")