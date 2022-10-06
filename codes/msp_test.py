import more_itertools as mit
import numpy as np


a = np.random.random((50))

n_points_msp_msh = 3

np_msp_bool_array = a > 0.3

ind_np_msp_vals = np.flatnonzero(np.convolve(np_msp_bool_array > 0,
                                 np.ones(n_points_msp_msh, dtype=int),
                                 'valid') >= n_points_msp_msh)

result = list(mit.run_length.encode(np_msp_bool_array))


#Find the length of longest subsequence of True, and the location if that index in result
max_true_count = -1
max_true_idx  = -1
for idx, (val, count) in enumerate(result):
    if val and max_true_count < count:
        max_true_count = count
        max_true_idx = idx

#Find total elements before and after the longest subsequence tuple
elems_before_idx = sum((idx[1] for idx in result[:max_true_idx]))
elems_after_idx = sum((idx[1] for idx in result[max_true_idx+1:]))

#Create the output list using the information
output = [(False, elems_before_idx), (True, max_true_count), (False, elems_after_idx)]
print(output)