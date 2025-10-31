#   This file is part of the OpenPhase (R) software library.
#  
#  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
#                Universitaetsstrasse 150, D-44801 Bochum, Germany
#            AND 2018-2025 OpenPhase Solutions GmbH,
#                Universitaetsstrasse 136, D-44799 Bochum, Germany.
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#     
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Standard makefile for OpenPhase executable

.PHONY : all clean cleanall

DEPTH = ../..

SRC := $(wildcard *.cpp)
PROGS := $(basename $(SRC))

include $(DEPTH)/Makefile.defs

INCLUDES += ${EXTERNAL_INCLUDES}
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
	rm -f *.bak *.o .g* .d* *.log *.sim $(PROGS)
	rm -rf VTK/
	rm -rf RawData/
	rm -rf TextData/
