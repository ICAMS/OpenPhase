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


# Root make file for OpenPhase project.

.PHONY : all clean cleanall clean_doc examples exercises benchmarks docs
SETTINGS = default
include Makefile.defs

PYTHON ?= python3


all: Makefile
	$(call print_header,Release)
#	$(call build_info_header)
#	@echo "Extracting libraries: Do not abort process!"
#	if [ ! -d "external/xyz" ]; then \
#		tar -xzf external/xyz.tar.gz -C external; fi
#	@echo "Libraries extracted."
#	@echo "---------------------------------"
	@echo "Compiling external libraries:"
	@echo "---------------------------------"
	make -C external --no-print-directory
	@echo "---------------------------------"
ifneq ($(findstring mpi-parallel, $(SETTINGS)),)
	@echo "Compiling MPI wrapper:"
	@echo "---------------------------------"
	make -C mpi_wrapper --no-print-directory
	@echo "---------------------------------"
endif
	@echo "Compiling main library:"
	@echo "---------------------------------"

	$(MAKE) --no-print-directory -C src

ifeq ($(findstring mpi-parallel, $(SETTINGS)),)
	@$(PYTHON) -c "import pybind11" 2>/dev/null \
	&& $(MAKE) --no-print-directory -C pythonbindings \
	|| (echo 'pybind11 not found. Skipping pythonbindings compilation. Install pybind11 if you need pythonbindings and re-run make' && exit 0)
endif
	$(call print_ending,$(COMPILE_CURRENT_TIME),$(COMPILE_START_TIME))

examples:
	make -C examples

benchmarks:
	make -C benchmarks

maintained:
	$(call print_header,Release)
#	@echo "Extracting libraries: Do not abort process!"
#	if [ ! -d "external/xyz" ]; then \
#		tar -xzf external/xyz.tar.gz -C external; fi
#	@echo "Libraries extracted."
#	@echo "---------------------------------"
#	@echo "Compile external libraries."
#	@echo "---------------------------------\n"
	make -C external --no-print-directory
#	@echo "---------------------------------\n"
	@echo "Compile main library."
	$(MAKE) --no-print-directory -C src
	make maintained --no-print-directory -C benchmarks
	make maintained --no-print-directory -C examples
	$(call print_ending,$(COMPILE_CURRENT_TIME),$(COMPILE_START_TIME))
docs:
	doxygen

clean:
	rm -rf hdf5 HighFive tinyxml2
	rm -rf external/miniz.o external/WinBase64/libbase64.o
	make --no-print-directory -C mpi_wrapper clean
	make --no-print-directory -C src clean
	make --no-print-directory -C benchmarks clean
	make --no-print-directory -C examples clean

cleanall:
	rm -rf hdf5 HighFive tinyxml2
	rm -rf external/miniz.o external/WinBase64/libbase64.o
	make --no-print-directory -C mpi_wrapper clean
	make --no-print-directory -C src clean
	make --no-print-directory -C benchmarks cleanall
	make --no-print-directory -C examples cleanall

clean_doc:
	-rm -r ./documentation/doxygen
