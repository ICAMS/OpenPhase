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

#include "GrandPotential/GrandPotentialPhaseDensity_IdealGas.h"
#include "PhysicalConstants.h"
#include "Settings.h"

namespace openphase
{
void GrandPotentialPhaseDensity_IdealGas::Initialize(Settings& locSettings)
{
    GrandPotentialPhaseDensity::Initialize(locSettings);
    MolarMass.assign(Ncomp, 0.0);
}
void GrandPotentialPhaseDensity_IdealGas::ReadInput(std::stringstream& InputFile, int moduleLocation)
{
    ConsoleOutput::Write("");
    ConsoleOutput::Write("Ideal gas grand potential density properties of "+PhaseNames[PhaseIdx]);
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        MolarMass[comp] = FileInterface::ReadParameterD(InputFile, moduleLocation, "M_"+ElementNames[comp]);
        if (MolarMass[comp] <= 0) ConsoleOutput::WriteWarning("MOS needs to positive!",thisclassname,"Read Input");
    }
    GravityAcceleration = FileInterface::ReadParameterD(InputFile, moduleLocation, "g");
    ConsoleOutput::Write("");
}
}// namespace openphase
