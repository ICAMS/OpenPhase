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


# Collective makefile for OpenPhase example

# List of directories in examples/

EXAMPLES = $(sort $(dir $(wildcard */)))

MAINTAINED=${filter-out\
           AdditiveNiAl/\
           Bainite_FID/\
           CPmulti/\
           CoSuperalloys/\
           EutecticII/\
           EutecticPtC/\
           Flow/\
           FlowMovingDendrites/\
           HeatDiffusion/\
           LatentHeat/\
           InterfaceStressPhaseTransformation/\
           MartensiteEXP/\
           MartensiteFeC_2D/\
           MartensiteTempering/\
           MultigrainGeneration/\
           MultigrainPlasticity/\
           OxideMulti/\
           Parallelism/\
           Pearlite/\
           Piezoelectric/\
           Recrystallization/\
           ResolvedTwin/\
           SharpPhaseField_ShrinkingSphere/\
           SolidSolidInteraction/\
           SolidificationAlCuMg/\
           SolidificationAlMn/\
           SolidificationFe/\
           SolidificationFeC/\
           SolidificationMgAl/\
           SolidificationMgAl.Old/\
           SolidificationNiAl/\
           SolidificationNiAlwithElasticity/\
           SuperalloysDD/\
           SuperalloysSSDD/\
           SuperalloysSimple/\
           ThermalGrooving.Alt/\
           , ${EXAMPLES}}

ifeq ($(MAKECMDGOALS),clean)
	FLAG = clean
endif
ifeq ($(MAKECMDGOALS),cleanall)
	FLAG = cleanall
endif
SETTINGS = default
.PHONY : all clean cleanall $(EXAMPLES)

all: $(EXAMPLES)
clean: $(EXAMPLES)
cleanall: $(EXAMPLES)
maintained: $(MAINTAINED)

$(EXAMPLES):
	-make SETTINGS="$(SETTINGS)" $(FLAG) -C $@
