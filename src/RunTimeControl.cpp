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

 *   File created :   2018
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#include "RunTimeControl.h"
#include "Settings.h"

namespace openphase
{

using namespace std;

void RunTimeControl::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "RunTimeControl";
    thisobjectname = thisclassname + ObjectNameSuffix;

    MaxTimeStep           = 0;
    dt                    = 0.0;
    SimulationTime        = 0.0;
    RestartSwitch         = false;
    StartTimeStep         = 0;
    TimeStep              = 0;
    CheckpointInterval    = 1;
    VTKOutputInterval     = 1;
    ConsoleOutputInterval = 1;
    OpenMPThreads         = 1;

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;
    TextDir = locSettings.TextDir;

    CheckTime = mygettime();
    StopTrigger = false;

    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void RunTimeControl::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("RunTimeControl input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::string filetype = FileInterface::getFileExtension(InputFileName); 
	if (filetype == "opi")
	{
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
	}
	else
	if (filetype == "json")
	{
		ReadJSON(InputFileName);
	}
	else
	{
		std::cerr << "Filetype " << filetype << " not recognized. Filetype must be opi or json." << std::endl;
		OP_Exit(1);
	}
}

void RunTimeControl::ReadInput(stringstream& inp)
{
    int moduleLocation    = FileInterface::FindModuleLocation(inp, thisclassname);

    SimulationTitle       = FileInterface::ReadParameterF(inp, moduleLocation, {string("SimulationTitle"), string("SimTtl")});

    UnitsOfLength         = FileInterface::ReadParameterS(inp, moduleLocation, {string("UnitsOfLength"), string("LUnits")}, false, "m");
    UnitsOfMass           = FileInterface::ReadParameterS(inp, moduleLocation, {string("UnitsOfMass"), string("MUnits")}, false, "kg");
    UnitsOfTime           = FileInterface::ReadParameterS(inp, moduleLocation, {string("UnitsOfTime"), string("TUnits")}, false, "s");
    UnitsOfEnergy         = FileInterface::ReadParameterS(inp, moduleLocation, {string("UnitsOfEnergy"), string("EUnits")}, false, "J");

#ifdef _OPENMP
    OpenMPThreads         = FileInterface::ReadParameterI(inp, moduleLocation, {string("OpenMPThreads"), string("nOMP")});
    omp_set_num_threads(OpenMPThreads);
#endif

    MaxTimeStep           = FileInterface::ReadParameterI(inp, moduleLocation, {string("MaxTimeStep"), string("nSteps")});
    dt                    = FileInterface::ReadParameterD(inp, moduleLocation, string("dt"));

    RestartSwitch         = FileInterface::ReadParameterB(inp, moduleLocation, {string("RestartSwitch"), string("Restrt")}, false, false);
    if(RestartSwitch)
    {
        StartTimeStep     = FileInterface::ReadParameterI(inp, moduleLocation, {string("StartTimeStep"), string("RestartTimeStep"), string("tStart")});
    }
    else
    {
        StartTimeStep     = 0;
    }
    CheckpointInterval    = FileInterface::ReadParameterI(inp, moduleLocation, {string("CheckpointInterval"), string("tRstrt")});
    VTKOutputInterval     = FileInterface::ReadParameterI(inp, moduleLocation, {string("VTKOutputInterval"), string("FTime")});
    ConsoleOutputInterval = FileInterface::ReadParameterI(inp, moduleLocation, {string("ConsoleOutputInterval"), string("STime")});
    StopCheckInterval     = FileInterface::ReadParameterD(inp, moduleLocation, string("StopCheckInterval"), false, 0.0);

    LogModeScreen         = FileInterface::ReadParameterB(inp, moduleLocation, string("LogScreen"),  false, false);
    LogModeVTK            = FileInterface::ReadParameterB(inp, moduleLocation, string("LogVTK"),     false, false);
    LogScreenFactor       = FileInterface::ReadParameterI(inp, moduleLocation, string("LogScreenF"), false, 1);
    LogVTKFactor          = FileInterface::ReadParameterI(inp, moduleLocation, string("LogVTKF"),    false, 1);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLine();

    if (VTKOutputInterval < 1)
    {
        ConsoleOutput::WriteSimple(string("Illegal value for VTKOutputInterval => ") + to_string(VTKOutputInterval) + string(" ! VTKOutputInterval is set to 1"));
        VTKOutputInterval = 1;
    }
    if (ConsoleOutputInterval < 1)
    {
        ConsoleOutput::WriteSimple(string("Illegal value for ConsoleOutputInterval => ") + to_string(ConsoleOutputInterval) + string(" ! ConsoleOutputInterval is set to 1"));
        ConsoleOutputInterval = 1;
    }

    if (std::filesystem::create_directories(VTKDir))
    {
        ConsoleOutput::WriteStandard("Created directory", VTKDir);
    }

    if (std::filesystem::create_directories(RawDataDir))
    {
        ConsoleOutput::WriteStandard("Created directory", RawDataDir);
    }

    if (std::filesystem::create_directories(TextDir))
    {
        ConsoleOutput::WriteStandard("Created directory", TextDir);
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void RunTimeControl::ReadJSON(const string InputFileName)
{

	std::ifstream f(InputFileName);
	json data = json::parse(f);
	if (data.contains(thisclassname))
	{
		json RTC = data[thisclassname];

    	SimulationTitle       = FileInterface::ReadParameter<std::string>(RTC, {"SimulationTitle"}, "NoTitle");

		UnitsOfLength         = FileInterface::ReadParameter<std::string>(RTC, {"UnitsOfLength"}, "m");
		UnitsOfMass           = FileInterface::ReadParameter<std::string>(RTC, {"UnitsOfMass"}, "kg");
		UnitsOfTime           = FileInterface::ReadParameter<std::string>(RTC, {"UnitsOfTime"}, "s");
		UnitsOfEnergy         = FileInterface::ReadParameter<std::string>(RTC, {"UnitsOfTime"}, "J");

	#ifdef _OPENMP
		OpenMPThreads         = FileInterface::ReadParameter<int>(RTC, {"OpenMPThreads"}, 1);
		omp_set_num_threads(OpenMPThreads);
	#endif

		MaxTimeStep           = FileInterface::ReadParameter<int>(RTC, {"MaxTimeStep"});
		dt                    = FileInterface::ReadParameter<double>(RTC, {"dt"});

		RestartSwitch         = FileInterface::ReadParameter<bool>(RTC, {"RestartSwitch"}, false);
		if(RestartSwitch)
		{
		    StartTimeStep     = FileInterface::ReadParameter<int>(RTC, {"StartTimeStep"}, 0);
		}
		else
		{
		    StartTimeStep     = 0;
		}
		CheckpointInterval    = FileInterface::ReadParameter<int>(RTC, {"CheckpointInterval"}, 1);
		VTKOutputInterval     = FileInterface::ReadParameter<int>(RTC, {"VTKOutputInterval"}, 1);
		ConsoleOutputInterval = FileInterface::ReadParameter<int>(RTC, {"ConsoleOutputInterval"}, 1);
		StopCheckInterval     = FileInterface::ReadParameter<double>(RTC, {"StopCheckInterval"}, 0.0);

		LogModeScreen         = FileInterface::ReadParameter<bool>(RTC, {"LogScreen"}, false);
		LogModeVTK            = FileInterface::ReadParameter<bool>(RTC, {"LogVTK"}, false);
		LogScreenFactor       = FileInterface::ReadParameter<double>(RTC, {"LogScreenF"}, 1);
		LogVTKFactor          = FileInterface::ReadParameter<double>(RTC, {"LogVTKF"}, 1);

		ConsoleOutput::WriteLine();
		ConsoleOutput::WriteLine();

		if (VTKOutputInterval < 1)
		{
		    ConsoleOutput::WriteSimple(string("Illegal value for VTKOutputInterval => ") + to_string(VTKOutputInterval) + string(" ! VTKOutputInterval is set to 1"));
		    VTKOutputInterval = 1;
		}
		if (ConsoleOutputInterval < 1)
		{
		    ConsoleOutput::WriteSimple(string("Illegal value for ConsoleOutputInterval => ") + to_string(ConsoleOutputInterval) + string(" ! ConsoleOutputInterval is set to 1"));
		    ConsoleOutputInterval = 1;
		}

		if (std::filesystem::create_directories(VTKDir))
		{
		    ConsoleOutput::WriteStandard("Created directory", VTKDir);
		}

		if (std::filesystem::create_directories(RawDataDir))
		{
		    ConsoleOutput::WriteStandard("Created directory", RawDataDir);
		}

		if (std::filesystem::create_directories(TextDir))
		{
		    ConsoleOutput::WriteStandard("Created directory", TextDir);
		}
	}
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

RunTimeControl& RunTimeControl::operator= (const RunTimeControl& rhs)
{
    if (this != &rhs) // protect against self-assignment
    {
        thisclassname         = rhs.thisclassname;

        SimulationTitle       = rhs.SimulationTitle;

        UnitsOfLength         = rhs.UnitsOfLength;
        UnitsOfMass           = rhs.UnitsOfMass;
        UnitsOfTime           = rhs.UnitsOfTime;
        UnitsOfEnergy         = rhs.UnitsOfEnergy;

        dt                    = rhs.dt;
        TimeStep              = rhs.TimeStep;
        SimulationTime        = rhs.SimulationTime;
        MaxTimeStep           = rhs.MaxTimeStep;
        VTKOutputInterval     = rhs.VTKOutputInterval;
        ConsoleOutputInterval = rhs.ConsoleOutputInterval;
        StartTimeStep         = rhs.StartTimeStep;
        RestartSwitch         = rhs.RestartSwitch;
        CheckpointInterval    = rhs.CheckpointInterval;
        OpenMPThreads         = rhs.OpenMPThreads;

        VTKDir                = rhs.VTKDir;
        RawDataDir            = rhs.RawDataDir;
        TextDir               = rhs.TextDir;

#ifdef _OPENMP
        omp_set_num_threads(OpenMPThreads);
#endif
    }
    return *this;
}
void RunTimeControl::Write()
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        std::ofstream outp(RawDataDir + thisobjectname + ".dat");
        outp << TimeStep << std::flush;
    }
}
bool RunTimeControl::RestartPossible()
{
    std::ifstream inp(RawDataDir + thisobjectname + ".dat");
    if (inp)
    {
        inp >> StartTimeStep;
        RestartSwitch = true;
    }
    else
    {
        RestartSwitch = false;
    }
    return RestartSwitch;
}

} //namespace openphase
