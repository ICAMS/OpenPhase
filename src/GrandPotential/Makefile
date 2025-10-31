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


# Depth relative to Openphase main folder
DEPTH = ../..

# Include makefile definitions
include $(DEPTH)/Makefile.defs

# Object SUBDIRS
OBJDIR := $(DEPTH)$(OBJDIR)
# Finds cpp files in src/ folder and remove path
SRC := $(wildcard *.cpp)
SRC := $(notdir $(SRC))
# Exlude entries from list above
SRC := $(filter-out $(EXCLUDE), $(SRC))

# List of object files to create
OBJS = $(SRC:%.cpp=$(OBJDIR)/%.o)
# List of directories in src/
SUBDIRS = $(sort $(dir $(wildcard */)))
# List of header files in include and its subfolders
HEADERS = $(wildcard $(DEPTH)/include/*.h) $(wildcard $(DEPTH)/include/**/*.h)
# Dependency files
DEPS = $(patsubst %.o,%.d,$(OBJS))

LIB = $(DEPTH)/lib/$(LIBNAME)

.PHONY : all clean subdirs $(SUBDIRS)

all: subdirs $(LIB)
subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

ifneq ($(findstring static, $(SETTINGS)),)
$(LIB): $(OBJS)
	ar $(AROPTIONS) $(LIB) $(OBJS)
else
$(LIB): $(OBJS) $(SUBDIRS)
	$(CXX) $(SOOPTIONS) -Wl,-soname=$(LIBNAME) -o $(LIB) $(OFILES) $(STDLIBS)
endif

$(OBJDIR)/%.o: %.cpp Makefile
	$(CXX) -MMD $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:

# Include List of dependency files in /obj folder. Will not complain
# if not present at first compilation.
-include $(DEPS)
