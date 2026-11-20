NVCC ?= nvcc
CXX ?= g++
NVCCFLAGS ?= -O2 -std=c++14
CXXFLAGS ?= -O2 -std=c++14

.PHONY: all main cpu clean
all: main
main: build/fdtd build/fdtd-strip

build:
	mkdir -p $@

build/fdtd: main.cu output.h | build
	$(NVCC) $(NVCCFLAGS) -o $@ main.cu

build/fdtd-strip: main.cu output.h | build
	$(NVCC) $(NVCCFLAGS) -DSTRIP -o $@ main.cu

# Original CPU implementation, retained as a reference (legacy text output).
cpu: build/fdtd-cpu
build/fdtd-cpu: fdtd.cpp | build
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -rf build
