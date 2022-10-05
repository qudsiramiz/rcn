# Import the required packages
import datetime
import itertools
import multiprocessing as mp
import time
import warnings

import geopack.geopack as gp
import h5py as hf
import numpy as np
import pyspedas as spd
import pytplot as ptt

# Define the function to compute magnetic fields from different models for different values of
# x_gsm, y_gsm and z_gsm


def model_run(*args):
    """
    Returns the value of the magnetic field at a given point in the model grid using three different
    models
    """
    i = args[0][0]
    j = args[0][1]
    k = args[0][2]

    bx_t96, by_t96, bz_t96 = gp.t96.t96(param, ps, x_gsm[i], y_gsm[j], z_gsm[k])
    bx_t01, by_t01, bz_t01 = gp.t01.t01(param, ps, x_gsm[i], y_gsm[j], z_gsm[k])
    bx_igrf, by_igrf, bz_igrf = gp.igrf_gsm(x_gsm[i], y_gsm[j], z_gsm[k])

    return i, j, k, bx_t96, by_t96, bz_t96, bx_t01, by_t01, bz_t01, bx_igrf, by_igrf, bz_igrf


now = time.time()

today_date = datetime.datetime.today().strftime("%Y-%m-%d")

# Define the time interval for data acquisition
trange = ["2016-12-24 15:08:00", "2016-12-24 15:12:00"]

# Try to get the solar wind data from OMNI, if it fails, set the parameter to a default value
try:
    # Define the omni data level
    omni_level = "hro"  # Valid options: "hro", "hro2" (default options)

    # Download the OMNI data (default level of "hro_1min") for the specified timerange.
    if ("omni_vars" in locals()):
        # If the data has already been downloaded, skip the download step
        print("Data already downloaded for the present date. Skipping the download process\n")
        pass
    else:
        omni_vars = spd.omni.data(trange=trange, level=omni_level)

    # Extract the proton density, magnetic field and velocity components from the OMNI data
    omni_time = ptt.get_data("BX_GSE")[0]
    # omni_time_unix = np.vectorize(datetime.utcfromtimestamp)(omni_time[:]) # converting omni_time
    # from unixtime to utc datetime object array in python

    omni_bx_gse = ptt.get_data("BX_GSE")[1]
    omni_by_gse = ptt.get_data("BY_GSE")[1]
    omni_bz_gse = ptt.get_data("BZ_GSE")[1]
    omni_by_gsm = ptt.get_data("BY_GSM")[1]
    omni_bz_gsm = ptt.get_data("BZ_GSM")[1]
    omni_np = ptt.get_data("proton_density")[1]
    omni_vx = ptt.get_data("Vx")[1]
    omni_vy = ptt.get_data("Vy")[1]
    omni_vz = ptt.get_data("Vz")[1]
    omni_sym_h = ptt.get_data("SYM_H")[1]

    # Compute the average time and magnetic field components for the specified time interval
    time_imf = np.nanmedian(omni_time)
    b_imf_x = np.nanmedian(omni_bx_gse)
    b_imf_y = np.nanmedian(omni_by_gsm)
    b_imf_z = np.nanmedian(omni_bz_gsm)

    # Raise warning if the average z-component of the field is outside the range of the model
    # parameters
    if (b_imf_z > 15 or b_imf_z < -18):
        warnings.warn(
            f"The given parameters produced the z-component of IMF field (b_imf_z) {b_imf_z} nT,"
            f"which is out of range in which model is valid (-18 nT < b_imf_z < 15 nT)"
        )

    # Find median values of other parameters
    np_imf = np.nanmedian(omni_np)
    vx_imf = np.nanmedian(omni_vx)
    vy_imf = np.nanmedian(omni_vy)
    vz_imf = np.nanmedian(omni_vz)
    sym_h_imf = np.nanmedian(omni_sym_h)

    # If any of the parameters are not finite set them to a default value
    if ~(np.isfinite(np_imf)):
        np_imf = 5
    if ~(np.isfinite(vx_imf)):
        np_imf = 500
    if ~(np.isfinite(sym_h_imf)):
        np_imf = -1

    m_p = 1.672e-27  # Mass of proton in SI unit
    rho = np_imf * m_p * 1.15
    # Solar wind ram pressure# in nPa, including roughly 4% Helium++ contribution
    p_dyn = 1.6726e-6 * 1.15 * np_imf * (vx_imf**2 + vy_imf**2 + vz_imf**2)

    # Raise warning if the solar dynamic pressure is outside the range of the model
    if (p_dyn > 8.5 or p_dyn < 0.5):
        warnings.warn(
            f"The given parameters produced a dynamic pressure of {p_dyn} nPa which is out of range"
            f"in which model is valid (0.5 nPa < p_dyn < 8.5 nPa)",
        )
    param = [p_dyn, sym_h_imf, b_imf_y, b_imf_z, 0, 0, 0, 0, 0, 0]
except Exception as e:
    print(f"OMNI data not found because of following exception: {e}, setting IMF parameters to "
          "default values")
    param = [5, -30, -1, -1, 0, 0, 0, 0, 0, 0]

# Compute the dipole tilt angle
ps = gp.recalc(time_imf)

# Define the maximum and minimum values of the model grid in terms of earth radii
x_min = -15.1
x_max = 15
y_min = -15.1
y_max = 15
z_min = -15.1
z_max = 15

# Define the resolution of the model grid in terms of earth radii
dr = 0.5

# Define the length of array based on the resolution of the model
n_arr = int(30/dr) + 1

x_gsm = np.linspace(x_min, x_max, n_arr)
y_gsm = np.linspace(y_min, y_max, n_arr)
z_gsm = np.linspace(z_min, z_max, n_arr)

# Define the pooling process
p = mp.Pool()  # number of cores, cannot exceed the maximum number of cores available

# Define the dimension of input/output
dim1, dim2, dim3 = n_arr, n_arr, n_arr

input = ((i, j, k) for i, j, k in itertools.product(range(dim3), repeat=3))

print("Running the model\n")
res = p.map(model_run, input)
print("Model run complete")

p.close()
p.join()

bx_t96 = np.empty((dim1, dim2, dim3))
by_t96 = np.empty((dim1, dim2, dim3))
bz_t96 = np.empty((dim1, dim2, dim3))

bx_t01 = np.empty((dim1, dim2, dim3))
by_t01 = np.empty((dim1, dim2, dim3))
bz_t01 = np.empty((dim1, dim2, dim3))

bx_igrf = np.empty((dim1, dim2, dim3))
by_igrf = np.empty((dim1, dim2, dim3))
bz_igrf = np.empty((dim1, dim2, dim3))

# Extract the data from the result
for r in res:
    i = r[0]
    j = r[1]
    k = r[2]

    bx_t96[i, j, k] = r[3]
    by_t96[i, j, k] = r[4]
    bz_t96[i, j, k] = r[5]

    bx_t01[i, j, k] = r[6]
    by_t01[i, j, k] = r[7]
    bz_t01[i, j, k] = r[8]

    bx_igrf[i, j, k] = r[9]
    by_igrf[i, j, k] = r[10]
    bz_igrf[i, j, k] = r[11]

bx_t96_total = bx_igrf + bx_t96
by_t96_total = by_igrf + by_t96
bz_t96_total = bz_igrf + bz_t96

bx_t01_total = bx_igrf + bx_t01
by_t01_total = by_igrf + by_t01
bz_t01_total = bz_igrf + bz_t01

# Save the data to an HDF file
version_number = 1
# Defien the filename
fname = f"../data/total_fields_0.5_15_15_{param}_{today_date}_v{str(version_number).zfill(2)}.h5"

data_file = hf.File(fname, "w")
data_file.create_dataset("bx_t96", data=bx_t96)
data_file.create_dataset("by_t96", data=by_t96)
data_file.create_dataset("bz_t96", data=bz_t96)

data_file.create_dataset("bx_t01", data=bx_t01)
data_file.create_dataset("by_t01", data=by_t01)
data_file.create_dataset("bz_t01", data=bz_t01)

data_file.create_dataset("bx_igrf", data=bx_igrf)
data_file.create_dataset("by_igrf", data=by_igrf)
data_file.create_dataset("bz_igrf", data=bz_igrf)


data_file.create_dataset("bx_t96_total", data=bx_t96_total)
data_file.create_dataset("by_t96_total", data=by_t96_total)
data_file.create_dataset("bz_t96_total", data=bz_t96_total)

data_file.create_dataset("bx_t01_total", data=bx_t01_total)
data_file.create_dataset("by_t01_total", data=by_t01_total)
data_file.create_dataset("bz_t01_total", data=bz_t01_total)

data_file.close()

print(f"Data saved to file {fname}")

print(f"took {round(time.time() - now, 3)} seconds")
