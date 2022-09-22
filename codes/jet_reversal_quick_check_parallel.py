import datetime
import importlib
import multiprocessing as mp
import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout

import pandas as pd
import pytz

import jet_reversal_quick_check_function as jrcf

importlib.reload(jrcf)

# Read the data from csv files
df_crossings = pd.read_csv("../data/mms_magnetopause_crossings.csv")
# Set the index to the date column
df_crossings.set_index("DateStart", inplace=True)

# for xx, crossing_time in enumerate(df_crossings.index[indx_number:indx_max], start=indx_number):

def check_jet_reversal(crossing_time):
    # Convert the crossing time to a datetime object
    # TODO: Something weird is happening with the timestamp. Check it later: crossing_time =
    # '2017-01-02 02:58:13.0+00:00'
    crossing_time = datetime.datetime.strptime(crossing_time.split('+')[0], "%Y-%m-%d %H:%M:%S.%f")
    # Set the timezone to UTC
    crossing_time = crossing_time.replace(tzinfo=pytz.utc)
    # Try with 'brst' data rate, if that fails then try with 'fast'
    inputs = {'crossing_time': crossing_time,
              'dt': 120,
              'probe': 3,
              'jet_len': 3,
              'level': 'l2',
              'coord_type': 'lmn',
              'data_type': 'dis-moms',
              'time_clip': True,
              'latest_version': True,
              'figname': 'mms_jet_reversal_check_lmn_mean',
              'fname': '../data/mms_jet_reversal_times_list_20220922.csv',
              #'fname': '../data/test.csv',
              'error_file_log_name': "../data/mms_jet_reversal_check_error_log_20220922.csv",
              "verbose": False
        }
    inputs["data_rate"] = 'brst'
    # df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)
    try:
        try:
            inputs["data_rate"] = 'brst'
            df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)
        except:
            inputs["data_rate"] = 'fast'
            df_fpi, df_fgm, df_mms = jrcf.jet_reversal_check(**inputs)
    except Exception as e:
        # print(f"\033[91;31m\n{e} for date {crossing_time}\n\033[0m")
        # Save the crossing time to a file
        # Check if the file exists
        if not os.path.isfile(inputs["error_file_log_name"]):
            # If it doesn't exist, create it
            with open(inputs["error_file_log_name"], 'w') as f:
                f.write(f"DateStart,Error\n")
                f.write(f"{crossing_time},{e}\n")
        else:
            # If it exists, append to it
            with open(inputs["error_file_log_name"], 'a') as f:
                f.write(f"{crossing_time},{e}\n")
        f.close()
        pass


@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

with suppress_stdout_stderr():
# for xxx in range(1):
    # Set the number of processes to use
    # num_processes = 20
    # Ask the user for index number
    indx_min = int(input("Enter the index number: "))
    #indx_min = 400
    # Ask the user for the maximum index number
    # indx_max = int(input("Enter the maximum index number: "))
    indx_max = indx_min + 5000
    # create a pool of processes
    pool = mp.Pool()
    # create a list of processes to run
    processes = [pool.apply_async(check_jet_reversal, args=(crossing_time,)) for crossing_time in df_crossings.index[indx_min:indx_max]]
    # run the processes
    for p in processes:
        p.get()
    # close the pool and wait for the processes to finish
    pool.close()
    pool.join()

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

