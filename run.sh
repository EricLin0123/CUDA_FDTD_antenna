make clean

echo 'compiling...'
make main

echo 'running...'
cd junk
../fdtd
cd ..

echo 'plotting...'
python plotCPU.py

echo 'convert to video'
python convertVideo.py
