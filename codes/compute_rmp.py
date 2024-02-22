import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Activate the latex text rendering
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

# Equation for pdyn
# p_dyn = 1.6726e-6 * 1.15 * np_imf * v_imf**2
# where p_dyn is in nPa, np_imf is in cm^-3, v_imf is in km/s


def get_magnetopause(sw_params):
    """
    NOTE: This function is only valid under the following conditions:
        1. p_dyn is in the range of 0.1 to 10 nPa
    """
    # Shue et al.,1998, equation 9
    ro = (10.22 + 1.29 * np.tanh(0.184 * (sw_params["b_imf"][2] + 8.14))) * (
        sw_params["p_dyn"]
    ) ** (-1.0 / 6.6)

    # Shue et al.,1998, equation 11
    alpha = (0.58 - 0.007 * sw_params["b_imf"][2]) * (
        1 + 0.024 * np.log(sw_params["p_dyn"])
    )
    rmp = ro * (2 / (1 + np.cos(0.0))) ** alpha  # Stand off position of the magnetopause

    return rmp


# Create a grid of points with bz in range of -5 to 5 and p_dyn in range of 0.1 to 10
# and compute the magnetopause for each point
# Plot the results
bz_grid = np.linspace(-10, -9, 100)
p_dyn_grid = np.linspace(0.1, 40, 100)
bz_grid, p_dyn_grid = np.meshgrid(bz_grid, p_dyn_grid)

rmp_grid = np.zeros(bz_grid.shape)
for i in range(bz_grid.shape[0]):
    for j in range(bz_grid.shape[1]):
        sw_params = {
            "b_imf": np.array([0.0, 0.0, bz_grid[i, j]]),
            "p_dyn": p_dyn_grid[i, j],
        }
        rmp_grid[i, j] = get_magnetopause(sw_params)

# Find the grid location where Rmp is closest to 6.1 RE
rmp_diff = np.abs(rmp_grid - 6.1)
rmp_diff_min = np.nanmin(rmp_diff)
rmp_diff_min_loc = np.where(rmp_diff == rmp_diff_min)
# print(f"rmp_diff_min = {rmp_diff_min:.2f} RE")
# print(f"rmp_diff_min_loc = {rmp_diff_min_loc}")
print(f"rmp = {rmp_grid[rmp_diff_min_loc]}")
print(f"bz = {bz_grid[rmp_diff_min_loc]}")
print(f"p_dyn = {p_dyn_grid[rmp_diff_min_loc]}")


"""
# if np_imf = 30 cm^-3, find the velocity of the solar wind at which the magnetopause is at 6.1 RE
# p_dyn = 1.6726e-6 * 1.15 * np_imf * v_imf**2
v_imf = np.sqrt(50 / (1.6726e-6 * 1.15 * 30))


print(f"v_imf = {v_imf:.2f} km/s")

# Plot the results
fig = plt.figure()
ax = fig.gca(projection="3d")
ax.plot_surface(bz_grid, p_dyn_grid, rmp_grid)
ax.set_xlabel("$B_{\rm z}$[nT]")
ax.set_ylabel("$P_{\rm dyn}$[nPa]")
ax.set_zlabel("$R_{\rm mp}$ [$R_\oplus$]")
# denote the Rmp = 6.1 RE line
ax.plot(
    bz_grid[0, :],
    p_dyn_grid[0, :],
    rmp_grid[0, :] * 0 + 6.1,
    color="red",
    linestyle="--",
    linewidth=2,
)

# Display the minimum and maximum values of Rmp on the plot
ax.text2D(
    0.05,
    0.95,
    f"$R_{{mp}}$ min = {np.min(rmp_grid):.2f} $R_\oplus$\n$R_{{mp}}$ max = {np.max(rmp_grid):.2f} $R_\oplus$",
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment="top",
)

# Set the viewing angle
ax.view_init(90, 0)

# Color the points where Rmp < 6.1 RE
rmp_grid[rmp_grid >= 6.1] = np.nan
ax.plot_surface(bz_grid, p_dyn_grid, rmp_grid, color="red", alpha=1)


plt.show()
"""
