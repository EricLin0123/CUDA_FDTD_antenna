import numpy as np
import matplotlib.pyplot as plt
import io
from tqdm import tqdm
import multiprocessing as mp

frames = 1600

# Constants
LIMX = 60
LIMY = 100

output_video = 'movie.mp4'

delX = 0.389e-3
delY = 0.4e-3
delT = 0.441e-12

x = np.arange(1, LIMX + 1) * delX
y = np.arange(1, LIMY + 1) * delY

# Create figure for 3D plotting
fig = plt.figure(figsize=(10, 7.5))
ax = fig.add_subplot(111, projection='3d')

# Initialize video frames list
frames = []

# Initialize the patch position
B = np.zeros((LIMX, LIMY))
heightP = 0.2

# Microstrip
B[19:26, 0:50] = heightP  # Note: Python uses 0-based indexing

# Patch antenna
B[14:47, 50:90] = heightP

# Initial plot
X, Y = np.meshgrid(y, x)
surf = ax.plot_surface(Y, X, B, cmap='viridis')
ax.set_xlim([0, LIMY * delY])
ax.set_ylim([0, LIMX * delX])
ax.set_zlim([-1, 1])
ax.set_xlabel('x-axis (m)')
ax.set_ylabel('y-axis (m)')
ax.set_zlabel('E_z (V/m)')
ax.set_title(
    'Here is the patch and microstrip line.\nIn a second a pulse will be launched')

# Capture initial frame multiple times
for _ in range(25):  # 25 frames of the initial plot
    # Convert plot to image
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100)
    with open('buffer/initial_frame.png', 'wb') as f:
        f.write(buf.getvalue())
    buf.close()


def process_frame(iter):
    fig = plt.figure(figsize=(10, 7.5))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([0, LIMY * delY])
    ax.set_ylim([0, LIMX * delX])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('x-axis (m)')
    ax.set_ylabel('y-axis (m)')
    ax.set_zlabel('E_z (V/m)')

    # Read data file
    filename = f'junk/junk.{iter}'
    try:
        A = np.loadtxt(filename)
    except:
        print(f"Could not read file {filename}")
        return

    # Reshape data into B matrix
    B = np.zeros((LIMX, LIMY))
    for a in range(LIMX):
        B[a, :] = A[a*LIMY:(a+1)*LIMY]

    # Create new surface plot
    surf = ax.plot_surface(Y, X, B, cmap='viridis')

    # Set labels and limits
    ax.set_title(f'Time = {
                 iter*delT:.2e} sec\nPlot of z-component of E-field under Patch Antenna and Microstrip')

    # Convert plot to image and save
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100)
    with open(f'buffer/{iter}.png', 'wb') as f:
        f.write(buf.getvalue())
    buf.close()
    plt.close(fig)


# Main animation loop
with mp.Pool(processes=16) as pool:
    list(tqdm(pool.imap(process_frame, range(1600)),
         total=1600, desc="Processing frames"))

print("Done!")
