import matplotlib.pyplot as plt
import numpy as np

# Read the data from the file
with open('junk/incident.txt', 'r') as file:
    data0 = file.readlines()

with open('junk/incident_strip.txt', 'r') as file:
    data1 = file.readlines()

with open('junk/Gauss.txt', 'r') as file:
    Guass = file.readlines()

# Convert the data to a list of floats
data0 = [float(line.strip()) / 3 for line in data0]
data1 = [float(line.strip()) / 3 for line in data1]
Guass = [float(line.strip()) for line in Guass]

# Calculate the difference between data0 and data1
data_diff = [d0 - d1 for d0, d1 in zip(data0, data1)]

# Plot the time domain data
plt.figure()
plt.plot(data0, label='antenna', linestyle=':', color='blue')
plt.plot(data1, linestyle='--', label='Strip only', color='green')
plt.plot(data_diff, label='Net reflected', color='red')
plt.plot(Guass, label='Gauss', color='black')
plt.title('Incident Data Plot')
plt.xlabel('Index')
plt.ylabel('Value')
plt.grid(True)
plt.legend()

# Perform Fourier transform on data_diff
fft_data_diff = np.fft.fft(data_diff)
fft_gauss = np.fft.fft(Guass)
# fft_gauss = np.fft.fft(data1)
s = fft_data_diff / fft_gauss

# Calculate the frequency axis
time_interval = 0.441e-12  # seconds
fft_freq = np.fft.fftfreq(len(data_diff), d=time_interval)

# Select only the positive frequencies
positive_freq_indices = np.where(fft_freq >= 0)
fft_freq = fft_freq[positive_freq_indices]
fft_data_diff = np.abs(fft_data_diff[positive_freq_indices])
fft_gauss = np.abs(fft_gauss[positive_freq_indices])
s = np.abs(s[positive_freq_indices])
s = 20 * np.log10(s)

# Plot the frequency domain data
plt.figure()
plt.plot(fft_freq, fft_data_diff, color='blue', label='Antenna')
plt.plot(fft_freq, fft_gauss, color='green', label='Strip')
plt.plot(fft_freq, s, color='red', label='Return Loss')
plt.title('Frequency Spectrum of Net Reflected Data')
plt.xlim([0, 2e10])
plt.ylim([-20, 5])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.grid(True)
plt.legend()

plt.show()
