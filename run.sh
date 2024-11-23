echo 'compiling...'
make fdtd
echo 'running...'
cd junk
../fdtd
cd ..
echo 'plotting...'
python plot.py
