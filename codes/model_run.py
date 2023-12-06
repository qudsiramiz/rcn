import numpy as np


def model_run(*args):
    """
    Returns the value of the magnetic field at a given point in the model grid using three different
    models
    """
    j = args[0][0]
    k = args[0][1]
    dr = args[0][2]
    y_max = args[0][3]
    z_max = args[0][4]
    model_type = args[0][5]


    y0 =y_max - int(j * dr)
    z0 = z_max - int(k * dr)
    rp = np.sqrt(y0**2 + z0**2)  # Projection of r into yz-plane

    for index in range(0, 100):

        theta = index * d_theta
        r = ro * (2/(1 + np.cos(theta))) ** alpha
        zp = r * np.sin(theta)  # not really in z direction, but a distance in yz plane
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
            # print(index, rp, zp)
            # print(f'Value of theta = {theta}')

            y_coord = y0
            z_coord = z0
            x_shu = (r -m_p) * np.cos(theta)
            phi = np.arctan2(z0, y0)
            # print( j, k, theta, x_shu[j,k])

            if (abs(y0) == 0 or abs(z0) == 0):
                if(abs(y0) == 0):
                    y_shu = 0
                    z_shu = (r -m_p) * np.sin(theta)
                elif (abs(z0) == 0):
                    z_shu = 0
                    y_shu = (r -m_p) * np.sin(theta)
            else:
                z_shu = np.sqrt((rp - 1.0)**2/(1 + np.tan(phi)**(-2)))
                y_shu = z_shu/np.tan(phi)

            rho_sh = rho * (1.509 * np.exp(x_shu/rmp) + .1285)
            n_sh = rho_sh/m_proton

            y_shu = abs(y_shu)*signy
            z_shu = abs(z_shu)*signz

            # Cooling JGR 2001 Model, equation 9 to 12
            # the distance from the focus to the magnetopause surface
            ll = 3 * rmp/2 - x0
            b_msx = - A * (- b_imf_x * (1 - rmp / (2 * ll)) + b_imf_y * (y0 / ll)
                           + b_imf_z * (z0 / ll))
            b_msy = A * (- b_imf_x * (y0 / (2 * ll)) + b_imf_y * (2 - y0**2/(
                           ll * rmp)) - b_imf_z * (y0 * z0 / (ll * rmp)))
            b_msz = A * (- b_imf_x * (z0 / (2 * ll)) - b_imf_y * (y0 * z0 / (ll
                         * rmp)) + b_imf_z * (2 - z0**2 / (ll * rmp)))
            try:
                if model_type == 't96':
                    bx_ext, by_ext, bz_ext = gp.t96.t96(param, ps, x_shu, y_shu, z_shu)
                elif model_type == 't01':
                    bx_ext, by_ext, bz_ext = gp.t01.t01(param, ps, x_shu, y_shu, z_shu)
            except:
                    print(f'Skipped for {x_shu, y_shu, z_shu}')
                    pass

            bx_igrf, by_igrf, bz_igrf = gp.igrf_gsm(x_shu, y_shu, z_shu)

            bx = bx_ext + bx_igrf
            by = by_ext + by_igrf
            bz = bz_ext + bz_igrf

            #if (np.sqrt(y_shu**2 + z_shu**2) > 31):
            #    shear = np.nan
            #    rx_en = np.nan
            #    va_cs = np.nan
            #    bisec_msp = np.nan
            #    bisec_msh = np.nan
            #else:
            shear = get_shear([bx, by, bz], [b_msx, b_msy, b_msz], angle_unit="degrees")

            rx_en = get_rxben([bx, by, bz], [b_msx, b_msy, b_msz])
            va_cs = get_vcs([bx, by, bz], [b_msx, b_msy, b_msz], n_sh, 0.1)
            bisec_msp, bisec_msh = get_bis([bx, by, bz], [b_msx, b_msy, b_msz])
            break

    return j, k, bx, by, bz, shear, rx_en, va_cs, bisec_msp, bisec_msh
