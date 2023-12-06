import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Activate the latex text rendering
plt.rc("text", usetex=True)
plt.rc("font", family="serif")


def get_magnetopause(sw_params):
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
bz_grid = np.linspace(-10, 10, 100)
p_dyn_grid = np.linspace(0.1, 20, 100)
bz_grid, p_dyn_grid = np.meshgrid(bz_grid, p_dyn_grid)

rmp_grid = np.zeros(bz_grid.shape)
for i in range(bz_grid.shape[0]):
    for j in range(bz_grid.shape[1]):
        sw_params = {
            "b_imf": np.array([0.0, 0.0, bz_grid[i, j]]),
            "p_dyn": p_dyn_grid[i, j],
        }
        rmp_grid[i, j] = get_magnetopause(sw_params)

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

plt.show()
