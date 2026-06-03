.PHONY: all build clean

DEPTH = ..
include $(DEPTH)/Makefile.defs

PYTHON ?= python3
MODULE = OpenPhase$(shell $(PYTHON)-config --extension-suffix)
PYBIND11_INCLUDES = $(shell $(PYTHON) -m pybind11 --includes)
LOCAL_RUNPATH = -Wl,-rpath='$$ORIGIN/../lib' -Wl,-z,origin
CONDA_RUNPATH = -Wl,-rpath='$(CONDA_LIB_DIR)'

all: build

build:
	$(CXX) $(CXXFLAGS) -shared -fPIC $(PYBIND11_INCLUDES) OpenPhase.cpp -o $(MODULE) $(INCLUDES) $(EXTERNAL_INCLUDES) -L$(DEPTH)/lib -lOpenPhase $(EXTERNAL_LIBS) $(LOCAL_RUNPATH) $(CONDA_RUNPATH)

clean:
	rm -f $(MODULE)
