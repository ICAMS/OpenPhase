#   This file is part of the OpenPhase (R) software library.
#
#   Copyright (c) 2009-2020 Ruhr-Universitaet Bochum,
#                 Universitaetsstrasse 150, D-44801 Bochum, Germany
#             AND 2018-2020 OpenPhase Solutions GmbH,
#                 Wasserstrasse 494, D-44795 Bochum, Germany.
#
#    All rights reserved.
#
#
#    DEVELOPMENT VERSION, DO NOT PUBLISH OR DISTRIBUTE.
#
#
#   OpenPhase (R) is a joint development of Interdisciplinary Centre for
#   Advanced Materials Simulation (ICAMS), Ruhr University Bochum
#   and OpenPhase Solutions GmbH.

# Standard makefile for OpenPhase executable

.PHONY : all clean cleanall

DEPTH = ../..

SRC := $(wildcard *.cpp)
PROGS := $(basename $(SRC))

include $(DEPTH)/Makefile.defs

INCLUDES = -I$(DEPTH)/include ${EXTERNAL_INCLUDES}
LDFLAGS = -L$(DEPTH)/lib ${EXTERNAL_LIBS}
LIBSRC = $(DEPTH)/lib/$(LIBNAME)

MAKE := make SETTINGS="$(SETTINGS)"


# Make target for compilation of the example
all: $(PROGS)

%: %.cpp $(SRC) $(LIBSRC) Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS) $(LIBSOP) $(STDLIBS) $(RUNPATH)

# Make target for removing of compiled executables

clean:
	rm -rf *.bak *.o .g* .d* $(PROGS)

# Make target for removing all traces of older compilations and simulation results

cleanall:
	rm -f *.bak *.o .g* .d* *.log $(PROGS)
	rm -rf VTK/
	rm -rf RawData/
	rm -rf TextData/
