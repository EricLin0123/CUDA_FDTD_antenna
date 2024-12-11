main: main.cu
	nvcc -o fdtd main.cu -lcuda -lcudart
	
fdtd: fdtd.cpp
	gcc -o fdtd fdtd.cpp -lm -lstdc++

PHONY: clean
clean:
	rm -f fdtd fdtd