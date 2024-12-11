# Author: 林萬荃 (Eric Lin)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

dx = 0.389
dy = 0.4
dz = 0.265

x_offset = (60 - 32) / 2 * dx
y_offset = 0
z_offset = 3 * dz
L1 = 32 * dx
L2 = 40 * dy
L3 = 50 * dy  # microstrip line length
L4 = 5 * dx
L5 = 6 * dx  # microstrip width

x_boundary = 60 * dx
y_boundary = 100 * dy
z_boundary = 3 * dz


antenna_vertices = [
    [0 + x_offset, L3 + y_offset, 0 + z_offset],
    [0 + x_offset, L3 + L2 + y_offset, 0 + z_offset],
    [L1 + x_offset, L2 + L3 + y_offset, 0 + z_offset],
    [L1 + x_offset, L3 + y_offset, 0 + z_offset],
    [L4 + L5 + x_offset, L3 + y_offset, 0 + z_offset],
    [L4 + L5 + x_offset, 0 + y_offset, 0 + z_offset],
    [L4 + x_offset, 0 + y_offset, 0 + z_offset],
    [L4 + x_offset, L3 + y_offset, 0 + z_offset]
]

center_x = (antenna_vertices[1][0] + antenna_vertices[2][0]) / 2
center_y = (antenna_vertices[1][1] + antenna_vertices[2][1]) / 2
center_z = (antenna_vertices[1][2] + antenna_vertices[2][2]) / 2
# Add the text at the center point
ax.text(center_x, center_y, center_z, f"{L1}mm", color='black')
# Calculate the center point between antenna_vertices[2] and antenna_vertices[3]
center_x = (antenna_vertices[2][0] + antenna_vertices[3][0]) / 2
center_y = (antenna_vertices[2][1] + antenna_vertices[3][1]) / 2
center_z = (antenna_vertices[2][2] + antenna_vertices[3][2]) / 2
# Add the text at the center point
ax.text(center_x, center_y, center_z, f"{L2}mm", color='black')
# Calculate the center point between antenna_vertices[5] and antenna_vertices[6]
center_x = (antenna_vertices[5][0] + antenna_vertices[6][0]) / 2
center_y = (antenna_vertices[5][1] + antenna_vertices[6][1]) / 2
center_z = (antenna_vertices[5][2] + antenna_vertices[6][2]) / 2
# Add the text at the center point
ax.text(center_x, center_y, center_z, f"{L5}mm", color='black')


# Add faces to the antenna
antenna_faces = [[antenna_vertices[j] for j in [0, 1, 2, 3, 4, 5, 6, 7]]]
copper_color = '#ffc600'
# Plot the antenna face
ax.add_collection3d(Poly3DCollection(
    antenna_faces, facecolors=copper_color, linewidths=1, edgecolors='#939393', alpha=1))


vertices = [[0, 0, 0], [x_boundary, 0, 0], [x_boundary, y_boundary, 0], [0, y_boundary, 0],
            [0, 0, z_boundary], [x_boundary, 0, z_boundary], [x_boundary, y_boundary, z_boundary], [0, y_boundary, z_boundary]]

# Define the six faces of the cube
faces = [[vertices[j] for j in [0, 1, 2, 3]],  # Bottom face
         [vertices[j] for j in [4, 5, 6, 7]],  # Top face
         [vertices[j] for j in [0, 1, 5, 4]],  # Front face
         [vertices[j] for j in [2, 3, 7, 6]],  # Back face
         [vertices[j] for j in [0, 3, 7, 4]],  # Left face
         [vertices[j] for j in [1, 2, 6, 5]]]  # Right face

# Plot the faces of the cube
ax.add_collection3d(Poly3DCollection(faces, facecolors='#ffffff',
                    linewidths=1, edgecolors='black', alpha=.25))

ax.add_collection3d(Poly3DCollection(
    [faces[0]], facecolors=copper_color, linewidths=1, edgecolors='#939393', alpha=0))


# Set plot limits and labels

ax.set_box_aspect([1, 1, 1])  # Aspect ratio is 1:1:1

# Set the same scale for all axes
max_range = max(x_boundary, y_boundary, z_boundary)
ax.set_xlim([0, max_range])
ax.set_ylim([0, max_range])
ax.set_zlim([0, max_range])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
