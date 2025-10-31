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

 *   File created :   2017
 *   Main contributors :   Oleg Shchyglo
 *
 */

#include "SymmetryVariants.h"
#include "Settings.h"

namespace openphase
{
using namespace std;

SymmetryVariants::SymmetryVariants(Settings& locSettings,
                                                const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void SymmetryVariants::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "SymmetryVariants";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    TransformationMatrices.resize(Nphases);
    for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        TransformationMatrices[pIndex].resize(locSettings.Nvariants[pIndex], dMatrix3x3::UnitTensor());
    }
    set = false;
    initialized = true;
    ConsoleOutput::WriteStandard("SymmetryVariants", "Initialized");
}

void SymmetryVariants::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteStandard("Source", InputFileName);
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File " + InputFileName + " could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ReadInput(data);
}

void SymmetryVariants::ReadInput(stringstream& inp)
{
    assert(initialized);

    ConsoleOutput::WriteLineInsert("Crystal symmetry transformation matrices");
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    // Reading variants for each phase
    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        for(size_t n = 0; n < TransformationMatrices[pIndex].size(); n++)
        {
            stringstream converter;
            converter << "V" << "_" << pIndex << "_" << n;

            TransformationMatrices[pIndex][n] = FileInterface::ReadParameterM3x3(inp, moduleLocation, converter.str(),false,dMatrix3x3::UnitTensor());
            set = true;
        }
    }
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

}// namespace openphase
