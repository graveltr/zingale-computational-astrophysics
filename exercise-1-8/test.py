import matplotlib.pyplot as plt
import numpy as np

# Generate sample data
x = np.linspace(0, 2 * np.pi, 100)  # 100 points from 0 to 2*pi
y = np.sin(x)

# Create a plot
plt.plot(x, y, label="Sine Wave")

# Add labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Simple Sine Wave Plot')
plt.legend()

# Show the plot
plt.show()
