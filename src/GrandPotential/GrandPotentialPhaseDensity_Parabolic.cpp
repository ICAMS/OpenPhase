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

#include "GrandPotential/GrandPotentialPhaseDensity_Parabolic.h"
#include "Settings.h"

namespace openphase
{
void GrandPotentialPhaseDensity_Parabolic::Initialize(Settings& locSettings)
{
    GrandPotentialPhaseDensity::Initialize(locSettings);

    EPS.assign(Ncomp, 0.0);
    C0 .assign(Ncomp, 0.0);
}
void GrandPotentialPhaseDensity_Parabolic::ReadInput (std::stringstream& InputFile, int moduleLocation)
{
    ConsoleOutput::Write("");
    ConsoleOutput::Write("Parabolic grad potential density material properties for "+PhaseNames[PhaseIdx]);
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+std::to_string(PhaseIdx)+"_"+ElementNames[comp];
        C0[comp] = FileInterface::ReadParameterD(InputFile, moduleLocation, "C0_"+converter);
        if (C0[comp] < 0) ConsoleOutput::WriteWarning("C0 needs to positive!",thisclassname,"Read Input");
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+std::to_string(PhaseIdx)+"_"+ElementNames[comp];
        EPS[comp] = FileInterface::ReadParameterD(InputFile, moduleLocation, "EPS_"+converter, false, 0.0);
        if (EPS[comp]==0.0)
        {
            double BulkModulus = FileInterface::ReadParameterD(InputFile, moduleLocation, "K0_"+converter);
            EPS[comp] = BulkModulus/C0[comp]/C0[comp];
        }
        if (EPS[comp] <= 0) ConsoleOutput::WriteWarning("EPS needs to positive!",thisclassname,"Read Input");
    }
    ConsoleOutput::Write("");
}
}// namespace openphase::GrandPotential
