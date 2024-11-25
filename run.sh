make clean

echo 'compiling...'
make fdtd

echo 'running...'
cd junk
../fdtd
cd ..

echo 'plotting...'
python plotCPU.py

# echo 'convert to video'
# python convertVideo.py
