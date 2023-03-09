import matplotlib.pyplot as plt
import numpy as np

# define data
values = np.array([12, 10, 8, 6, 4])
colors = np.array(['r', 'b', 'g', 'y', 'm', 'c'])
labels = np.array(['A', 'B', 'C', 'D', 'E', 'F'])
radii = np.array([1, 4, 6, 8, 10, 12])

# Set latex font
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# set up figure
fig, ax = plt.subplots(figsize=(6, 6))
ax.axis('equal')

# draw circles
for i in range(len(radii)):
    circle = plt.Circle((0, radii[i]), radii[i], color=colors[i], alpha=1, fill=False, linewidth=2)
    ax.add_artist(circle)


# add labels
for i in range(len(radii)):
    plt.annotate(labels[i], xy=(0, 2*radii[i] - 0.5), ha='center', va='center')

# add legend
plt.legend(labels, loc='lower right')

# Set the axes limits
ax.set_xlim(-12, 12)
ax.set_ylim(0, 24)

# show plot
plt.savefig('../figures/nested_circles.pdf')
plt.close("all")
