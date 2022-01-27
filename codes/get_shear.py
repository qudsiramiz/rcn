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
        Preferred unit of angle returned by the code. Default is "radians".

    Raises
    ------
    KeyError If the key is not input_angle is not set to "radians" or "degrees" then the code raises
        a key error.

    Returns
    -------
    angle: float
        Angle between the two vectors, in radians by default
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
