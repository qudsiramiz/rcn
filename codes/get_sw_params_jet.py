import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import get_sw_fnc_w_mms as rmf
import importlib
import glob
import datetime

importlib.reload(rmf)

data_folder = '../data/rx_d'
'''
fname = np.sort(glob.glob(f"{data_folder}/reconnection_line_data_mms3_20221109.csv"))[0]
df = pd.read_csv(fname, index_col=False)
# Set date_from as index
df = df.set_index("date_from")

df_shear = df[df.method_used == "shear"].copy()

date_from_list = df_shear.index.tolist()
date_to_list = df_shear.date_to.tolist()

# Convert date_from_list and date_to_list to datetime
# date_from_list = [datetime.datetime.strptime(date_from, "%Y-%m-%d %H:%M:%S")
#                   for date_from in date_from_list]
# 
# date_to_list = [datetime.datetime.strptime(date_to, "%Y-%m-%d %H:%M:%S")
#                 for date_to in date_to_list]

key_list = ["time_imf", "b_imf", "b_imf_x", "b_imf_y", "b_imf_z", "rho", "ps", "p_dyn", "sym_h",
            "tp", "v_sw", "ca"]
# Open a file to write the output
fname_out = f"{data_folder}/sw_params_jet.csv"
f = open(fname_out, "w")
# Write the header from key_list
f.write(",".join(key_list) + "\n")

for date_from, date_to in zip(date_from_list[0:], date_to_list[0:]):
    sw_params = rmf.get_sw_params(trange=[date_from, date_to], omni_level='hro', verbose=True)
    # For each key in the dictionary, write the value to the file
    for key in sw_params.keys():
        f.write(str(sw_params[key]) + ",")
    f.write("\n")

f.close()
'''



fname = np.sort(glob.glob(f"{data_folder}/sw_params_jet.csv"))[0]
df = pd.read_csv(fname, index_col=False)
df = df.set_index("time_imf")

# Convert tp from K to eV
df["tp"] = df["tp"] * 8.617e-5

m_proton = 1.672e-27

# Convert rho to n
df["rho"] = df["rho"] / (m_proton * 1.15)

key_list = ["b_imf", "b_imf_x", "b_imf_y", "b_imf_z", "rho", "ps", "p_dyn", "sym_h",
            "tp", "v_sw", "ca"]

# For each key in the dataframe, print the 5th and 95th percentile and the median value in a table
print("Key\t5th\tMedian\t95th")
for key in df.keys():
    print(f"{key}\t{df[key].quantile(0.05):.2f}\t{df[key].quantile(0.5):.2f}\t{df[key].quantile(0.95):.2f}")


# Read column 16 of a txt file in a dataframe
# df = pd.read_csv(fname, index_col=False, usecols=[16], header=None)