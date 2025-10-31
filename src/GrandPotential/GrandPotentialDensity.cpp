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

#include "Settings.h"
#include "GrandPotential/GrandPotentialDensity.h"
#include "GrandPotential/GrandPotentialPhaseDensity_CompressibleFluid.h"
#include "GrandPotential/GrandPotentialPhaseDensity_IdealGas.h"
#include "GrandPotential/GrandPotentialPhaseDensity_ImplicitParabolic.h"
#include "GrandPotential/GrandPotentialPhaseDensity_Parabolic.h"
#include "VTK.h"

namespace openphase
{
GrandPotentialDensity::GrandPotentialDensity(Settings& locSettings, std::string filename)
{
    InitializeAndReadInput(locSettings, filename);
}
void GrandPotentialDensity::InitializeAndReadInput(Settings& locSettings, std::stringstream& inp_data)
{
    PhaseNames = locSettings.PhaseNames;
    Nphases    = locSettings.Nphases;

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert(thisclassname);

    int moduleLocation = FileInterface::FindModuleLocation(inp_data, thisclassname);

    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    {

        std::string converter = ""+std::to_string(PhaseIdx);
        std::string str_omega = FileInterface::ReadParameterK(inp_data, moduleLocation, "OMEGA_"+converter);
        if      (str_omega == "PARABOLIC" )         storage.push_back(new GrandPotentialPhaseDensity_Parabolic         (PhaseIdx));
        else if (str_omega == "IMPLICITPARABOLIC")  storage.push_back(new GrandPotentialPhaseDensity_ImplicitParabolic (PhaseIdx));
        else if (str_omega == "IDEALGAS")           storage.push_back(new GrandPotentialPhaseDensity_IdealGas          (PhaseIdx));
        else if (str_omega == "COMPRESSIBLEFLUID")  storage.push_back(new GrandPotentialPhaseDensity_CompressibleFluid (PhaseIdx));
        else
        {
            std::stringstream message;
            message << "Phase Model "+str_omega+" is not implemented!\n";
            ConsoleOutput::WriteExit(message.str(),thisclassname,"ReadInput");
            OP_Exit(EXIT_FAILURE);
        }
    }

    DoGravity = FileInterface::ReadParameterB(inp_data, moduleLocation, "Gravity", false, false);
    if (DoGravity)
    {
        GravityDirection = FileInterface::ReadParameterI(inp_data, moduleLocation, "GDir");
        for (auto omega : storage)
        {
            if (!omega->GravityImplemented)
            {
                std::stringstream message;
                message << "Gravity for " << omega->thisclassname << " is not yet implemented!\n";
                ConsoleOutput::WriteExit(message.str(),thisclassname,"ReadInput");
                OP_Exit(EXIT_FAILURE);
            }
        }
    }

    for (auto omega : storage) omega->Initialize(locSettings);
    for (auto omega : storage) omega->ReadInput(inp_data, moduleLocation);

    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}
void GrandPotentialDensity::InitializeAndReadInput(Settings& locSettings, std::string filename)
{
    std::fstream inp(filename.c_str(), std::ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + filename + "\" could not be opened", thisclassname, "InitializeAndReadInput");
        std::exit(EXIT_FAILURE);
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();
    ConsoleOutput::WriteStandard("Source", filename);
    InitializeAndReadInput(locSettings,inp_data);
}
}// namespace openphase::GrandPotential
