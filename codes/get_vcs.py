def get_vcs(b_vec_1, b_vec_2, n_1, n_2):
    r"""
    Get vcs code.

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
    vcs : float
        The exhaust velocity in km/s
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

    vcs = va_p1 * np.sqrt(rx_mag_1 * rx_mag_2 * (rx_mag_1 + rx_mag_2)/(rx_mag_1 * n_2 +
                                                                       rx_mag_2 * n_1))

    # vcs = va_p1 * np.sqrt(mag_vec_1 * mag_vec_2 * (mag_vec_1 + mag_vec_2)/(mag_vec_1 * n_2 +
    #                                                                         mag_vec_2 * n_1))

    return vcs
