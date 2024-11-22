import numpy as np
import matplotlib.pyplot as plt

# Parameters for the Gaussian distribution
mean = 0          # Mean (μ)
std_dev = 1       # Standard deviation (σ)

# Generate x values (range of the distribution)
x = np.linspace(mean - 4*std_dev, mean + 4*std_dev, 1000)

# Gaussian formula
y = (1 / (std_dev * np.sqrt(2 * np.pi))) * \
    np.exp(-0.5 * ((x - mean) / std_dev)**2)

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(x, y, label=f'Gaussian (μ={mean}, σ={std_dev})', color='blue')

# Adding labels, title, and legend
plt.title("Gaussian Distribution", fontsize=16)
plt.xlabel("x", fontsize=14)
plt.ylabel("Probability Density", fontsize=14)
plt.legend(fontsize=12)
plt.grid(alpha=0.3)

# Show the plot
plt.show()
