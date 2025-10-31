.PHONY: all build clean

all: build

build:
	g++ -O3 -Wall -shared -fopenmp -std=c++17 -fPIC $$(python3 -m pybind11 --includes) OpenPhase.cpp -o OpenPhase$$(python3-config --extension-suffix) -I../external -I../external/WinBase64 -I../include -L../lib -lOpenPhase ../external/miniz.o ../external/WinBase64/libbase64.o -Wl,-rpath='./' -Wl,-z,origin

clean:
	rm OpenPhase.cpython-39-x86_64-linux-gnu.so




