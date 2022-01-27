def get_rxben(b_vec_1, b_vec_2):
    r"""
    Get the reconnection energy between two magnetic field lines.

    It has the following mathematical expression:

    .. math:: rexben = 0.5 (|\vec{B_1}| + \vec{B_2}) (1 - \hat{B_1} \cdot \hat{B_2})

    Parameters
    ----------
    b_vec_1 : array of shape 1x3
        Input magnetic field vector.
    b_vec_2 : array of shape 1x3
        Input magnetic field vector.

    Returns
    -------
    rxben : float
        Reconnection field energy density in nPa
    """

    alpha = - 14.87 * np.pi / 180  # radians (From Hesse2013)
    b_vec_1 = np.array(b_vec_1)
    b_vec_2 = np.array(b_vec_2)
    mag_vec_1 = np.linalg.norm(b_vec_1)
    mag_vec_2 = np.linalg.norm(b_vec_2)

    unit_vec_1 = b_vec_1/mag_vec_1
    unit_vec_2 = b_vec_2/mag_vec_2

    # The bisector vector of thwo input vectors
    unit_vec_bisec = (unit_vec_1 + unit_vec_2) / np.linalg.norm(unit_vec_1 + unit_vec_2)

    # Cross product of the two input vectors with the bisector vector to get the reconnection
    # component of the magnetic field.
    rx_b_1 = np.cross(b_vec_1, unit_vec_bisec)
    rx_b_2 = np.cross(b_vec_2, unit_vec_bisec)

    rx_b_mag_1 = np.linalg.norm(rx_b_1)
    rx_b_mag_2 = np.linalg.norm(rx_b_2)

    # The guide field of the reconnection
    gd_b_1 = b_vec_1 - rx_b_1
    gd_b_2 = b_vec_2 - rx_b_2

    gd_b = gd_b_1 + gd_b_2
    gd_b_mag = np.linalg.norm(gd_b)

    # The reconnection energy
    b_u_prime = rx_b_mag_1 * np.cos(alpha) + gd_b_mag * np.sin(alpha)
    b_d_prime = rx_b_mag_2 * np.cos(alpha) + gd_b_mag * np.sin(alpha)

    # Reconnection energy
    #rx_en = b_u_prime ** 2 * b_d_prime ** 2
    rx_en = rx_b_mag_1 ** 2 * rx_b_mag_2 ** 2
    #unit_vec_1 = b_vec_1/mag_vec_1
    #unit_vec_2 = b_vec_2/mag_vec_2
#
    #u_bisect = (unit_vec_1 + unit_vec_2) / 2
    #rx_bmag1 = np.dot(u_bisect, b_vec_1)
    #rx_bmag2 = np.dot(u_bisect, b_vec_2)
    #b1_b2_dotp = np.dot(unit_vec_1, - unit_vec_2)
#
    ##rx_en = 0.5 * (mag_vec_1 * mag_vec_2) * (1 + b1_b2_dotp)
    #rx_en = (rx_bmag1**2 * rx_bmag2**2)   # nPa
    ##rx_en = (rx_bmag1**2 + rx_bmag2**2) * 1.03  # MJ/RE^3

    # Angle between the two vectors
    #angle_b1_b2 = get_shear(b_vec_1, b_vec_2, angle_unit="radians")
#
    #angle_bisect = angle_b1_b2/2.
#
    ## Reconnecting component of the magnetic field
    #rx_bmag1 = mag_vec_1 * np.sin(angle_bisect)
    #rx_bmag2 = mag_vec_2 * np.sin(angle_bisect)
#
    ## Non-reconnecting/guide field component of the magnetic field
    #gd_bmag1 = mag_vec_1 * np.cos(angle_bisect)
    #gd_bmag2 = mag_vec_2 * np.cos(angle_bisect)
#
    ##TODO: Check if this is correct
    ##gd_mag = np.sqrt(gd_bmag1**2 + gd_bmag2**2)  # Magnitude of guide field?
    #gd_mag = np.sqrt(gd_bmag1**2 + gd_bmag2**2)
    #alpha = 14.87 * np.pi/180.
#
    #b_u_prime = rx_bmag1 * np.cos(alpha)# + gd_mag * np.sin(alpha)
    #b_d_prime = rx_bmag2 * np.cos(alpha)# + gd_mag * np.sin(alpha)
#
    #rx_en = b_u_prime**2 * b_d_prime**2

    #angle_bisect = np.arccos(np.dot(b_vec_1, b_vec_2)/(np.linalg.norm(b_vec_1)*np.linalg.norm(b_vec_2)))

    return rx_en
