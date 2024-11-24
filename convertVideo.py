import imageio

# Define parameters
fps = 24  # Frames per second
output_filename = 'output_video.mp4'

# Generate the list of file names
image_files = [f"buffer/{i}.png" for i in range(600)]

# Create the MP4 video
with imageio.get_writer(output_filename, fps=fps, format='mp4') as writer:
    for file in image_files:
        image = imageio.imread(file)
        writer.append_data(image)

print(f"Video saved as {output_filename}")
