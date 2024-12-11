main: main.cu
	nvcc -o fdtd main.cu -lcuda -lcudart
	
fdtd: fdtd.cpp
	gcc -o fdtd_seq fdtd.cpp -lm


PHONY: clean
clean:
	rm -f fdtd fdtd_seq