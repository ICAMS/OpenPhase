.PHONY: all build clean

all: build

build:
	mpicxx -g --shared -std=c++17 -fPIC -Wl,-soname=libmpi_wrapper.so mpi_wrapper.cpp -o ../lib/libmpi_wrapper.so -I../include -lfftw3_mpi -lfftw3_omp -lfftw3

clean:
	rm -rf ../lib/libmpi_wrapper.so

