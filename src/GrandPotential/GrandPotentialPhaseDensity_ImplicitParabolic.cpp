/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *     
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *   File created :   2021
 *   Main contributors :   Raphael Schiedung
 *
 */

#include "GrandPotential/GrandPotentialPhaseDensity_ImplicitParabolic.h"
#include "Settings.h"

namespace openphase
{
void GrandPotentialPhaseDensity_ImplicitParabolic::Initialize(Settings& locSettings)
{
    GrandPotentialPhaseDensity_Parabolic::Initialize(locSettings);
    Zeta.assign(Ncomp, 0.0);
}
void GrandPotentialPhaseDensity_ImplicitParabolic::ReadInput(std::stringstream& InputFile, int moduleLocation)
{
    GrandPotentialPhaseDensity_Parabolic::ReadInput(InputFile,moduleLocation);
    ConsoleOutput::Write("");
    ConsoleOutput::Write("Implicit parabolic grand potential density material properties for "+PhaseNames[PhaseIdx]);
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+std::to_string(PhaseIdx)+"_"+ElementNames[comp];
        Zeta[comp] = FileInterface::ReadParameterD(InputFile, moduleLocation, "Zeta_"+converter);
        if (Zeta[comp] <= 0) ConsoleOutput::WriteWarning("Zeta needs to positive!",thisclassname,"Read Input");
    }
    std::string converter = ""+std::to_string(PhaseIdx);
    precision     = FileInterface::ReadParameterD(InputFile, moduleLocation, "PREC_"+converter);
    MaxIterations = FileInterface::ReadParameterI(InputFile, moduleLocation, "MAXI_"+converter);
    ConsoleOutput::Write("");
}
}// namespace openphase

