make clean

echo 'compiling...'
make main

echo 'running...'
cd junk
../fdtd
cd ..

echo 'plotting PNG ...'
python plotCPU.py

echo 'Using GPU H.264 hardware to convert to video'
/usr/bin/ffmpeg -hwaccel cuda -framerate 60 -pattern_type glob -i 'buffer/*.png' -c:v h264_nvenc -preset fast output.mp4
