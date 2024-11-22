import ctypes
import numpy as np
import matplotlib.pyplot as plt

# Load the shared library
gaussian_lib = ctypes.CDLL('./gaussian.so')

# Define the function prototype: double gaussian(double x, double mean, double std_dev)
gaussian_lib.gaussian.argtypes = [
    ctypes.c_float, ctypes.c_float, ctypes.c_float]
gaussian_lib.gaussian.restype = ctypes.c_float

# Parameters for the Gaussian distribution
mean = 0.0         # Mean (μ)
std_dev = 1.0      # Standard deviation (σ)

# Generate x values (range of the distribution)
x_values = np.linspace(mean - 4 * std_dev, mean + 4 * std_dev, 1000)

# Call the C Gaussian function for each x value
y_values = np.array([gaussian_lib.gaussian(x, mean, std_dev)
                    for x in x_values])

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values, label=f'Gaussian (μ={
         mean}, σ={std_dev})', color='blue')

# Adding labels, title, and legend
plt.title("Gaussian Distribution (C Implementation)", fontsize=16)
plt.xlabel("x", fontsize=14)
plt.ylabel("Probability Density", fontsize=14)
plt.legend(fontsize=12)
plt.grid(alpha=0.3)

# Show the plot
plt.show()
