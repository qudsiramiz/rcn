import datetime
import importlib
import multiprocessing as mp
import numpy as np
import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout

import pandas as pd
import pytz

import jet_reversal_quick_check_function_beta as jrcfb

importlib.reload(jrcfb)

# Read the data from csv files
df_crossings = pd.read_csv("../data/mms_magnetopause_crossings.csv", index_col=False)
# df_jet_reversal_times = pd.read_csv("../data/mms_jet_reversal_times_list_20221017_beta.csv",
#                                    index_col=False)
# Set the index to the date column
df_crossings.set_index("DateStart", inplace=True)
# df_jet_reversal_times.set_index("jet_time", inplace=True)
# Select only those times where index is unique
# df_jet_reversal_times = df_jet_reversal_times.loc[~df_jet_reversal_times.index.duplicated(
#                                                                                       keep='first')]

# for xx, crossing_time in enumerate(df_crossings.index[indx_number:indx_max], start=indx_number):

date_obs = "paper"


def check_jet_reversal(crossing_time):
    # Convert the crossing time to a datetime object
    # TODO: Something weird is happening with the timestamp. Check it later: crossing_time =
    # '2017-01-02 02:58:13.0+00:00'
    crossing_time = datetime.datetime.strptime(crossing_time.split('+')[0],
                                               "%Y-%m-%d %H:%M:%S.%f")
    # Set the timezone to UTC
    crossing_time = crossing_time.replace(tzinfo=pytz.utc)
    # Try with 'brst' data rate, if that fails then try with 'fast'
    inputs = {'crossing_time': crossing_time,
              'dt': 600,
              'probe': 3,
              'jet_len': 3,
              'level': 'l2',
              'coord_type': 'lmn',
              'data_type': ['dis-moms', 'des-moms'],
              'time_clip': True,
              'latest_version': True,
              'date_obs': date_obs,
              'figname': 'mms_jet_reversal_check_lmn_mean',
              'fname': f'../data/mms_jet_reversal_times_list_{date_obs}_beta_brst.csv',
              # 'fname': '../data/test.csv',
              'error_file_log_name': f"../data/mms_jet_reversal_check_err_log_{date_obs}_beta.csv",
              "verbose": True
              }
    # inputs["data_rate"] = 'fast'
    # # ptt, df_mms, ind_range_msp, ind_range_msh, t_jet_center = jrcfb.jet_reversal_check(**inputs)
    # _ =  jrcfb.jet_reversal_check(**inputs)
    # v1, v2, ind_walen = jrcfb.jet_reversal_check(**inputs)
    try:
        # try:
        inputs["data_rate"] = 'brst'
        _ = jrcfb.jet_reversal_check(**inputs)
        # except Exception:
        #     inputs["data_rate"] = 'brst'
        #     _ = jrcfb.jet_reversal_check(**inputs)
    except Exception as e:
        # print(f"\033[91;31m\n{e} for date {crossing_time}\n\033[0m")
        # Save the crossing time to a file
        # Check if the file exists
        if not os.path.isfile(inputs["error_file_log_name"]):
            # If it doesn't exist, create it
            with open(inputs["error_file_log_name"], 'w') as f:
                f.write("DateStart,Error\n")
                f.write(f"{crossing_time},{e}\n")
        else:
            # If it exists, append to it
            df_added_list = pd.read_csv(inputs["error_file_log_name"], sep=',', index_col=False)
            if not np.any(df_added_list['DateStart'].values == str(crossing_time)):
                with open(inputs["error_file_log_name"], 'a') as f:
                    f.write(f"{crossing_time},{e}\n")
                f.close()
        pass
    #return ptt, df_mms, ind_range_msp, ind_range_msh, t_jet_center
    return None


@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


use_parallel = False

#with suppress_stdout_stderr():
for foo in range(1):
    if use_parallel:
        # Set the number of processes to use
        # num_processes = 20
        # Ask the user for index number
        # indx_min = int(input("Enter the index number: "))
        indx_min = 18
        # indx_min = 400
        # Ask the user for the maximum index number
        # indx_max = int(input("Enter the maximum index number: "))
        indx_max = indx_min + 2

        # Run the jet reversal check in parallel
        # pool = mp.Pool()
        # pool.map(check_jet_reversal, df_crossings.index[indx_min:indx_max])
        # pool.close()
        # pool.join()
        # Create a list of processes that we want to run and run them in parallel
        processes = [mp.Process(target=check_jet_reversal, args=(crossing_time,))
                     for crossing_time in df_crossings.index[indx_min:indx_max]]
        # Run processes
        for p in processes:
            p.start()
        # Exit the completed processes
        for p in processes:
            p.join()
        # for p in processes:
        #     p.terminate()
        for p in processes:
            p.close()
        #pool = mp.Pool()
        ## create a list of processes to run
        #processes = [pool.apply_async(check_jet_reversal, args=(crossing_time,)) for
        #             crossing_time in df_crossings.index[indx_min:indx_max]]
        ## run the processes
        #for p in processes:
        #    p.get()
        ## close the pool and wait for the processes to finish
        #pool.close()
        #pool.join()
    else:
        indx_min = 16515
        # indx_min = 400
        # Ask the user for the maximum index number
        # indx_max = int(input("Enter the maximum index number: "))
        indx_max = indx_min + 1
        for xx, crossing_time in enumerate(df_crossings.index[indx_min:indx_max],
                                           start=indx_min):
            check_jet_reversal(crossing_time)
    # indx_number = 0
    # indx_max = indx_number + 200
    # if __name__ == '__main__':
    #     # Setup a list of processes that we want to run
    #     processes = [mp.Process(target=check_jet_reversal, args=(crossing_time,))
    #                  for crossing_time in df_crossings.index[indx_number:indx_max]]
    #     # Run processes
    #     for p in processes:
    #         p.start()
    #     # Exit the completed processes
    #     for p in processes:
    #         p.join()
