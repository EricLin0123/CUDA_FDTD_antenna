gauss: gauss.c
	gcc -shared -o gaussian.so -fPIC gauss.c -lm

fdtd: fdtd.cpp
	gcc -o fdtd fdtd.cpp -lm

main: main.cu
	nvcc -o fdtd main.cu -lcuda -lcudart

PHONY: clean
clean:
	rm -f gaussian.so fdtd main