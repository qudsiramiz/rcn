wimport matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(8, 8))

# Create a circle with radius 1
circle = plt.Circle((0, 0), 1, color='k', alpha=0.5, fill=False)

# Add circle to plot
ax.add_artist(circle)

# Set x and y limits to create a square plot
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)

# Create a rage of x and y values within the circle's radius
x = np.linspace(-1, 1, 1000)
y = np.sqrt(1 - x**2)

# Shade the lower h-af the circle
ax.fill_between(x, y, -y, where=x>=0, color='k', alpha=0.5)

plt.savefig('../figures/circle.png', dpi=300, bbox_inches='tight')
plt.close()
