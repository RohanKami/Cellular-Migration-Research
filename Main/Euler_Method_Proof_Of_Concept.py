import matplotlib.pyplot as plt
import numpy as np

# Function to calculate y_numerical and y_analytical
def calculate_values(deltaT):
    steps = int(1 / deltaT)
    t_values = np.arange(0, 1, deltaT)
    y_numerical = np.zeros(steps)
    y_analytical = np.zeros(steps)
    
    y = 0
    for i in range(steps):
        t = i * deltaT
        y += deltaT * -t**2
        y_numerical[i] = y
        y_analytical[i] = -1/3 * t**3
    
    return t_values, y_numerical, y_analytical

# Calculate values for different deltaT
deltaT_values = [0.1, 0.01, 0.001]
results = [calculate_values(deltaT) for deltaT in deltaT_values]

# Create the plot
plt.figure(figsize=(3.5, 2.5), dpi=300)

# Plot analytical solution once
t_values_analytical, _, y_analytical = results[-1]  # Get the last set for t_values and y_analytical
plt.plot(t_values_analytical, y_analytical, label='Analytical Solution', linestyle='--', color='black')

# Plot results for each deltaT
colors = ['blue', 'green', 'red']
for deltaT, (t_values, y_numerical, _), color in zip(deltaT_values, results, colors):
    plt.plot(t_values, y_numerical, label=f'Numerical (Δt={deltaT})', color=color)

# Add labels with LaTeX formatting, title, legend, and grid
plt.xlabel('Time, $t$ (h)', fontsize=10)
plt.ylabel('Value', fontsize=10)
plt.title('Comparison of Numerical Integration and Analytical Solution', fontsize=12)
plt.legend(loc='best', fontsize=8)
plt.grid(True)

# Improve the style of the grid and axes
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=8)

# Show the plot
plt.tight_layout()
plt.show()
