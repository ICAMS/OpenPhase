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

 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Marvin Tegeler;
 *
 */

#include "MovingFrame.h"
#include "Settings.h"

namespace openphase
{
using namespace std;

MovingFrame::MovingFrame()
{
    trigger_phase_idx = -1;
    trigger_position  = -1;
    moving_direction  = "NN";
}

void MovingFrame::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "MovingFrame";
    thisobjectname = thisclassname + ObjectNameSuffix;

    trigger_phase_idx = -1;
    trigger_position  = -1;
    moving_direction  = "NN";

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void MovingFrame::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("MovingFrame input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);
    inp.close();

    ConsoleOutput::WriteLine();
}

void MovingFrame::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    trigger_phase_idx = FileInterface::ReadParameterI(inp, moduleLocation, string("TriggerPhaseIndex"), false, -1);
    trigger_position  = FileInterface::ReadParameterI(inp, moduleLocation, string("TriggerPosition"), false, -1);
    moving_direction  = FileInterface::ReadParameterS(inp, moduleLocation, string("MovingDirection"), false, "No");

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

}// namespace openphase
