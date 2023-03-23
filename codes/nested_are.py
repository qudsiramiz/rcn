'''
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
'''
'''
import matplotlib.pyplot as plt
import numpy as np

# Generate some sample data
x = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
y = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# Set the size of the marker and the offset position
marker_size = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
offset_position = marker_size/2

# Plot the scatter plot with 'o' marker
# plt.scatter(x, y, s=marker_size**2, facecolors='none', edgecolors='b', linewidths=2)

# Adjust the position of the circles
for i in range(len(x)):
    plt.scatter(x[i], y[i], s=marker_size[i]**2, facecolors='none', edgecolors='b', offset_position=offset_position[i])

# Show the plot
plt.savefig('../figures/nested_circles.pdf')
plt.close("all")
'''


import matplotlib.pyplot as plt

x = [1, 2, 3, 4, 5]
y = [2, 3, 1, 4, 2]

fig, ax = plt.subplots()
ax.scatter(x, y, marker='o')

for i in range(len(x)):
    ax.annotate(f"({x[i]}, {y[i]})", (x[i], y[i]), xytext=(0, -10), textcoords="offset points", ha='center')
    
plt.show()
