import numpy as np
import tabulate

# Define a function to calculate the angle bisector of two vectors
def angle_bisector(v1, v2):
    # Normalize the vectors
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    # Calculate the angle bisector
    v3 = v1 + v2
    v3 = v3 / np.linalg.norm(v3)
    return v3


# Define two vectors generated from random numbers, with 3 components each and with a seed value of
# 0
for i in range(3):
    np.random.seed(i)
    v1 = np.random.rand(3) * 10
    v2 = np.random.rand(3) * 30

    # Calculate the angle bisector of the two vectors
    vb = angle_bisector(v1, v2)

    # Compute the cross product of v1 and v3
    v1_cross = np.cross(v1, vb)

    # Compute the cross product of v2 and v3
    v2_cross = np.cross(v2, vb)

    # Get the magnitude of the cross products
    v1_cross_mag = np.linalg.norm(v1_cross)
    v2_cross_mag = np.linalg.norm(v2_cross)

    # Compute the dot product of v1 and v3
    v1_dot = np.dot(v1, vb)
    v2_dot = np.dot(v2, vb)

    # Find the angle between v1_cross and v2_cross
    v1_v2_cross_angle = np.arccos(np.round(np.dot(v1_cross, v2_cross) / (v1_cross_mag * v2_cross_mag), 3))

    # Print all the values in the terminal in a fancy table with 2 decimal places
    print(tabulate.tabulate([['v1', f'{v1[0]:.2f}', f'{v1[1]:.2f}', f'{v1[2]:.2f}'],
                            ['v2', f'{v2[0]:.2f}', f'{v2[1]:.2f}', f'{v2[2]:.2f}'],
                            ['vb', f'{vb[0]:.2f}', f'{vb[1]:.2f}', f'{vb[2]:.2f}'],
                            ['v1_cross', f'{v1_cross[0]:.2f}', f'{v1_cross[1]:.2f}', f'{v1_cross[2]:.2f}'],
                            ['v2_cross', f'{v2_cross[0]:.2f}', f'{v2_cross[1]:.2f}', f'{v2_cross[2]:.2f}'],
                            ['v1_cross_mag', f'{v1_cross_mag:.2f}', '', ''],
                            ['v2_cross_mag', f'{v2_cross_mag:.2f}', '', ''],
                            ['v1_dot', f'{v1_dot:.2f}', '', ''],
                            ['v2_dot', f'{v2_dot:.2f}', '', ''],
                            ['v1_v2_cross_angle', f'{v1_v2_cross_angle:.2f}', '', '']],
                            headers=['Vector', 'x', 'y', 'z'], tablefmt='fancy_grid'))
