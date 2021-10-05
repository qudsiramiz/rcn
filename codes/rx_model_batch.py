# This the python version of IDL code named 'RX_model_batch.pro'
import geopack
from geopack import geopack
from geopack.t96 import t96
import numpy as np
import pyspedas as spd
import pytplot as ptt
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
from dateutil import parser
import datetime
import h5py as hf

def get_shear(b_vec_1, b_vec_2, angle_unit="radians"):
    r"""
    Get the shear angle between two magnetic field lines.

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        Input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Input magnetic field vector.
    angle_unit : str, optional
        Preferred unit of angle returned by the code. Ddefault is "radians".

    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
    angle: float
        Angle between the two vectors in radians by default
    """
    unit_vec_1 = b_vec_1/np.linalg.norm(b_vec_1)
    unit_vec_2 = b_vec_2/np.linalg.norm(b_vec_2)
    angle = np.arccos(np.dot(unit_vec_1, unit_vec_2))

    #dp = b_vec_1[0] * b_vec_2[0] + b_vec_1[1] * b_vec_2[1] + b_vec_1[2] * b_vec_2[2]

    #mag1 = np.linalg.norm( b_vec_1)
    #mag2 = np.linalg.norm( b_vec_2)
    #angle = np.arccos(dp/(mag1*mag2))

    if (angle_unit == "radians"):
        return angle
    elif (angle_unit == "degrees"):
        return angle * 180/np.pi
    else:
        raise KeyError("angle_unit must be radians or degrees")


def get_rxben(b_vec_1, b_vec_2):
    r"""
    Get rxben between two magnetic field lines.

    As of now I am not sure what this is, but it has the following mathematical expression:

    .. math:: rexben = 0.5 ( |\vec{B_1}| + \vec{B_2} ) (1 - \hat{B_1} \cdot \hat{B_2})

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        Input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Input magnetic field vector.
    angle_unit : str, optional
        Preferred unit of angle returned by the code. Default is "radians".

    Returns
    -------
    rxben : Find out
    """
    # TODO: Update the documentation of this function

    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)

    return 0.5 * (mag_vec_1 + mag_vec_2) * (1 + b1_b2_dotp)


def get_vcs(b_vec_1, b_vec_2, n_1, n_2):
    r"""
    Get vec code.

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        Input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Input magnetic field vector.
    n_1 : float
        Density of first component
    n_2 : float
        Density of second component

    Returns
    -------
    Some scalar

    """
    va_p1 = 21.812  # conv. nT, m_P/cm ^ 3 product to km/s cassak-shay

    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)

    # bisector = mag_vec_1*b_vec_1 + mag_vec_2*b_vec_2 u_bisect = bisector /
    # np.linalg.norm(bisector)

    rx_mag_1 = mag_vec_1 * (1 + b1_b2_dotp)/2
    rx_mag_2 = mag_vec_2 * (1 + b1_b2_dotp)/2

    return va_p1 * np.sqrt(rx_mag_1 * rx_mag_2 * (rx_mag_1 + rx_mag_2)/(rx_mag_1 * n_1 +
                                                                        rx_mag_2 * n_2))


def get_bis(b_vec_1, b_vec_2, angle_unit="radians"):
    r"""
    Get the shear angle between two magnetic field lines.

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        First input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Second input magnetic field vector.

    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
     Returns bisect stuff
    """
    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)

    angle = np.arccos(b1_b2_dotp)

    if (angle_unit == "radians"):
        return angle
    elif (angle_unit == "degrees"):
        return angle * 180/np.pi
    else:
        raise KeyError("angle_unit must be radians or degrees")


def get_ca(b_vec, angle_unit="radians"):
    r"""
    Get ca.

    Parameters
    ----------
    b_vec : array of shape 1x3
        Input magnetic field vector.
    angle_unit : str, optional
        Preferred unit of angle returned by the code. Default is "radians".

    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
    angle : float Returns arctan of y- and z-component of the magnetic field vector.

    """
    angle = np.arctan(b_vec[1]/b_vec[2])

    if (angle_unit == "radians"):
        return angle
    elif (angle_unit == "degrees"):
        return angle * 180/np.pi
    else:
        raise KeyError("angle_unit must be radians or degrees")


#def rx_model_batch(
#    probe=None,
#    omni_level='hro',
#    maximum_shear=True,
#    mms_probe='mms1',
#    movie=None,
#    trange=None
#):
for i in range(1):
    r"""
    RX model from the IDL code.

    Parameters
    ----------
    probe : TYPE, optional
        DESCRIPTION. The default is probe.
    maximum_shear : TYPE, optional
        DESCRIPTION. The default is maximum_shear.
    movie : TYPE, optional.
        The default is movie.
    mmsprobe : TYPE, optional
        DESCRIPTION. The default is mmsprobe.
    times : TYPE, optional
        DESCRIPTION. The default is times.

    Returns
    -------
    None.
    """
    probe = None
    omni_level = 'hro'
    maximum_shear = True
    mms_probe = 'mms1'
    movie = None
    trange = ['2016-12-24 15:08', '2016-12-24 15:12']

    # Get time range as a datetime object
    trange_unix = [parser.parse(xx) for xx in trange]
    # For MMS add 5 minutes buffer on either side of the time range and give out 'trange_mms' in a
    # format that can be read by PySpedas
    trange_mms_dt = [trange_unix[0] - datetime.timedelta(minutes=5), 
                     trange_unix[1] + datetime.timedelta(minutes=5)]
    trange_mms = [xx.strftime('%Y-%m-%d %H-%M-%S') for xx in trange_mms_dt]

    # Download the OMNI data (default level of 'hro_1min') for the specified timerange.
    if ('omni_vars' in locals()):
        pass
    else:
        omni_vars = spd.omni.data(trange=trange, level=omni_level)

    omni_time = ptt.get_data('BX_GSE')[0]
    # omni_time_unix = np.vectorize(datetime.utcfromtimestamp)(omni_time[:]) # converting omni_time
    # from unixtime to utc datetime object array in python

    omni_bx_gse = ptt.get_data('BX_GSE')[1]
    omni_by_gse = ptt.get_data('BY_GSE')[1]
    omni_bz_gse = ptt.get_data('BZ_GSE')[1]
    omni_by_gsm = ptt.get_data('BY_GSM')[1]
    omni_bz_gsm = ptt.get_data('BZ_GSM')[1]
    omni_np = ptt.get_data('proton_density')[1]
    omni_vx = ptt.get_data('Vx')[1]
    omni_vy = ptt.get_data('Vy')[1]
    omni_vz = ptt.get_data('Vz')[1]
    omni_sym_h = ptt.get_data('SYM_H')[1]

    if (probe is None):
        pass
    elif ~hasattr(probe, '__len__'):
        themis_stt_vars = spd.themis.state(probe=probe, trange=trange)
        themis_time = ptt.get_data(f'th{probe}_pos_gse')[0]
        themis_sc_pos = ptt.get_data(f'th{probe}_pos_gse')[1]
    elif hasattr(probe, '__len__'):
        themis_vars = []
        sc_pos = []
        for xx in probe:
            themis_vars.append(spd.themis.state(probe=xx, trange=trange))
            sc_pos.append(ptt.get_data(f'th{xx}_pos_gse')[1])

    if (mms_probe is not None):
        if('mms_mec_vars' in locals()):
            pass
        else:
            mms_mec_vars = spd.mms.mec(trange=trange_mms, data_rate='srvy', probe='1')
            mms_time = ptt.get_data('mms1_mec_r_gse')[0]
            mms_sc_pos = ptt.get_data('mms1_mec_r_gse')[1:3]
    else:
        pass

    if (movie is None):
        time_imf = np.nanmedian(omni_time)
        b_imf_x = np.nanmedian(omni_bx_gse)
        b_imf_y = np.nanmedian(omni_by_gsm)
        b_imf_z = np.nanmedian(omni_bz_gsm)

        if (b_imf_z > 15 or b_imf_z < -18):
            warnings.warn(
            f"The given parameters produced the z-component of IMF field (b_imf_z) {b_imf_z} nT,"
            f"which is out of range in which model is valid (-18 nT < b_imf_z < 15 nT)"
            )

        np_imf = np.nanmedian(omni_np)
        vx_imf = np.nanmedian(omni_vx)
        vy_imf = np.nanmedian(omni_vy)
        vz_imf = np.nanmedian(omni_vz)
        sym_h_imf = np.nanmedian(omni_sym_h)

        if ~(np.isfinite(np_imf)):
            np_imf = 5
        if ~(np.isfinite(vx_imf)):
            np_imf = 500
        if ~(np.isfinite(sym_h_imf)):
            np_imf = -1

        m_p = 1.672e-27  # Mass of proton in SI unit

        rho = np_imf * m_p * 1.15
        p_dyn = 1.6726e-6 * 1.15 * np_imf * (vx_imf**2 + vy_imf**2 + vz_imf**2) #  Solar wind ram
                                                                                #  pressure
                                                                                #  in nPa,
                                                                                #  including
        # roughly 4% Helium++ contribution

        if (p_dyn > 8.5 or p_dyn < 0.5):
            warnings.warn(
            f"The given parameters produced a dynamic pressure of {p_dyn} nPa which is out of"
            f" range in which model is valid (0.5 nPa < p_dyn < 8.5 nPa)",
            )
        param = [p_dyn, sym_h_imf, b_imf_y, b_imf_z, 0, 0, 0, 0, 0, 0]

        # Compute the dipole tilt angle
        ps = geopack.recalc(time_imf)

        n_arr = 150

        bx = np.full((n_arr, n_arr), np.nan)
        by = np.full((n_arr, n_arr), np.nan)
        bz = np.full((n_arr, n_arr), np.nan)

        bx_t96 = np.full((n_arr, n_arr), np.nan)
        by_t96 = np.full((n_arr, n_arr), np.nan)
        bz_t96 = np.full((n_arr, n_arr), np.nan)

        bx_igrf = np.full((n_arr, n_arr), np.nan)
        by_igrf = np.full((n_arr, n_arr), np.nan)
        bz_igrf = np.full((n_arr, n_arr), np.nan)

        b_msx = np.full((n_arr, n_arr), np.nan)
        b_msy = np.full((n_arr, n_arr), np.nan)
        b_msz = np.full((n_arr, n_arr), np.nan)

        x_shu = np.full((n_arr, n_arr), np.nan)
        y_shu = np.full((n_arr, n_arr), np.nan)
        z_shu = np.full((n_arr, n_arr), np.nan)

        rho_sh = np.full((n_arr, n_arr), np.nan)

        #rp = np.full((n_arr, n_arr), np.nan)

        # r = np.full((n_arr, n_arr), np.nan)
        # zp = np.full((n_arr, n_arr), np.nan)
        # x0 = np.full((n_arr, n_arr), np.nan)

        shear = np.full((n_arr, n_arr), np.nan)
        rx_en = np.full((n_arr, n_arr), np.nan)
        va_cs = np.full((n_arr, n_arr), np.nan)
        bisec = np.full((n_arr, n_arr), np.nan)
        y_coord = np.full((n_arr, n_arr), np.nan)
        z_coord = np.full((n_arr, n_arr), np.nan)
        b_sh_ca = np.full((n_arr, n_arr), np.nan)
        b_sh_mag = np.full((n_arr, n_arr), np.nan)
        n_sh = np.full((n_arr, n_arr), np.nan)

        d_theta = np.pi/100

        # Shue et al.,1998, equation 10
        #if (b_imf_z >=0):
        #    ro = (11.44 + 0.013 * b_imf_z) * (p_dyn)**(-1.0/6.6)
        #    #ro = (10.22 + 1.29 * np.tanh(0.184 * (b_imf_z + 8.14))) * (p_dyn)**(-1.0/6.6)
        #else:
        #    ro = (11.44 + 0.14 * b_imf_z) * (p_dyn)**(-1.0/6.6)
        ro = (10.22 + 1.29 * np.tanh(0.184 * (b_imf_z + 8.14))) * (p_dyn)**(-1.0/6.6)

        # Shue et al.,1998, equation 11
        #alpha = (0.58 - 0.010 * b_imf_z) * (1 + 0.010 * p_dyn)
        alpha = (0.58 - 0.007 * b_imf_z) * (1 + 0.024 * np.log(p_dyn))
        rmp = ro * (2/(1 + np.cos(0.0))) ** alpha  # Stand off position of the magnetopause

        A = 2
        len_y = 80
        len_z = 80
        count = 0
        for j in range(0, len_y):
            y0 = 40 - j
            for k in range(0, len_z):
                z0 = 40 - k
                rp = np.sqrt(y0**2 + z0**2)  # Projection of r into yz-plane

                for index in range(0, 100):

                    theta = index * d_theta
                    r = ro * (2/(1 + np.cos(theta))) ** alpha
                    zp = r * np.sin(theta)  # not really in z direction, but a
                                                        # distance in yz plane
                    x0 = r * np.cos(theta)

                    if x0 == 0:
                        signx = 1.0
                    else:
                        signx = np.sign(x0)

                    if y0 == 0:
                        signy = 1.0
                    else:
                        signy = np.sign(y0)

                    if z0 == 0:
                        signz = 1.0
                    else:
                        signz = np.sign(z0)

                    if (rp <= zp):
                        #print(index, rp, zp)
                        #print(f'Value of theta = {theta}')

                        y_coord[j, k] = y0
                        z_coord[j, k] = z0
                        #print( j, k, y0, z0)
                        x_shu[j, k] = (r - 0.5) * np.cos(theta)
                        phi = np.arctan2(z0, y0)
    
                        if (abs(y0) == 0 or abs(z0) == 0):
                            if(abs(y0) == 0):
                                y_shu[j, k] = 0
                                z_shu[j, k] = (r - 0.5) * np.sin(theta)
                            elif (abs(z0) == 0):
                                z_shu[j, k] = 0
                                y_shu[j, k] = (r - 0.5) * np.sin(theta)
                        else:
                            z_shu[j, k] = np.sqrt((rp - 0.5)**2/(1 + np.tan(phi)**(-2)))
                            y_shu[j, k] = z_shu[j, k]/np.tan(phi)
    
                        rho_sh[j, k] = rho * (1.509 * np.exp(x_shu[j, k]/rmp) + .1285)
    
                        y_shu[j, k] = abs(y_shu[j, k])*signy
                        z_shu[j, k] = abs(z_shu[j, k])*signz
    
                        # Cooling JGR 2001 Model, equation 9 to 12
                        ll = 3 * rmp/2 - x0  # the distance from the focus to the
                                                   # magnetopause surface
                        b_msx[j, k] = - A * (- b_imf_x * (1 - rmp / (2 * ll)) + b_imf_y * (y0 / ll)
                                      + b_imf_z * (z0 / ll))
                        b_msy[j, k] = A * (- b_imf_x * (y0 / (2 * ll)) + b_imf_y * (2 - y0**2/(
                                      ll * rmp)) - b_imf_z * (y0 * z0 / (ll * rmp)))
                        b_msz[j, k] = A * (- b_imf_x * (z0 / (2 * ll)) - b_imf_y * (y0 * z0 / (ll
                                      * rmp)) + b_imf_z * (2 - z0**2 / (ll * rmp)))
    
                        # TODO: Implement Geopack T96!!!
                        # Compute the external magnetic field from the T95 model for a given
                        # position and time, in GSM coordinate
                        bx_t96[j, k], by_t96[j, k], bz_t96[j, k] = t96(param, ps, x_shu[j, k],
                                                                       y_shu[j, k], z_shu[j, k])
    
                        # Compute the internal magnetic field from the IGRF model for a given
                        # position in GSM coordinate
                        bx_igrf[j, k], by_igrf[j, k], bz_igrf[j, k] = geopack.igrf_gsm(x_shu[j, k],
                                                                                       y_shu[j, k],
                                                                                       z_shu[j, k])
    
                        bx[j, k] = bx_t96[j, k] + bx_igrf[j, k]
                        by[j, k] = by_t96[j, k] + by_igrf[j, k]
                        bz[j, k] = bz_t96[j, k] + bz_igrf[j, k]
    
                        if (np.sqrt(y_shu[j, k]**2 + z_shu[j, k]**2) > 31):
                            shear[j, k] = np.nan
                        else:
                            shear[j, k] = get_shear([bx[j, k], by[j, k], bz[j, k]], [b_msx[j, k],
                                                   b_msy[j, k], b_msz[j, k]], angle_unit="degrees")
    
                        rx_en[j, k] = get_rxben([bx[j, k], by[j, k], bz[j, k]], [b_msx[j, k],
                                                 b_msy[j, k], b_msz[j, k]])
                        va_cs[j, k] = get_vcs([bx[j, k], by[j, k], bz[j, k]], [b_msx[j, k],
                                              b_msy[j, k], b_msz[j, k]], 0.1, rho_sh[j, k])
                        bisec[j, k] = get_bis([bx[j, k], by[j, k], bz[j, k]],
                                              [b_msx[j, k], b_msy[j, k], b_msz[j, k]])
                        b_sh_ca[j, k] = get_ca([b_msx[j, k], b_msy[j, k], b_msz[j, k]])
                        b_sh_mag[j, k] = np.linalg.norm([b_msx[j, k], b_msy[j, k], b_msz[j, k]])
                        n_sh[j, k] = rho_sh[j, k]
                        break

    # Save the outputs to a file
    data_file = hf.File('../data/all_data_rx_model_0.5re_20211005_v02.h5', 'w')

    data_file.create_dataset('bx', data=bx)
    data_file.create_dataset('by', data=by)
    data_file.create_dataset('bz', data=bz)

    data_file.create_dataset('bx_t96', data=bx_t96)
    data_file.create_dataset('by_t96', data=by_t96)
    data_file.create_dataset('bz_t96', data=bz_t96)
                                                  
    data_file.create_dataset('bx_igrf', data=bx_igrf)
    data_file.create_dataset('by_igrf', data=by_igrf)
    data_file.create_dataset('bz_igrf', data=bz_igrf)
                                                  
    data_file.create_dataset('b_msx', data=b_msx)
    data_file.create_dataset('b_msy', data=b_msy)
    data_file.create_dataset('b_msz', data=b_msz)
                                                      
    data_file.create_dataset('x_shu', data=x_shu)
    data_file.create_dataset('y_shu', data=y_shu)
    data_file.create_dataset('z_shu', data=z_shu)
                                                      
    data_file.create_dataset('rho_sh', data=rho_sh)
    #data_file.create_dataset('rp', data=rp)
                                                  
    #data_file.create_dataset('r', data=r)
    #data_file.create_dataset('zp', data=zp)
    #data_file.create_dataset('x0', data=x0)
                                                      
    data_file.create_dataset('shear', data=shear)
    data_file.create_dataset('rx_en', data=rx_en)
    data_file.create_dataset('va_cs', data=va_cs)
    data_file.create_dataset('bisec', data=bisec)
    data_file.create_dataset('y_coord', data=y_coord)
    data_file.create_dataset('z_coord', data=z_coord)
    data_file.create_dataset('b_sh_ca', data=b_sh_ca)
    data_file.create_dataset('b_sh_mag', data=b_sh_mag)
    data_file.create_dataset('n_sh', data=n_sh)

    data_file.close()

    # TODO: Code for plotting rx_en contours/energy density