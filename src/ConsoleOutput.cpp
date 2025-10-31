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

 *   File created :   2013
 *   Main contributors :   Philipp Engels; Raphael Schiedung
 *
 */

#include "ConsoleOutput.h"
#include "RunTimeControl.h"
#include "Includes.h"
#include "BuildInfo.h"

namespace openphase
{

using namespace std;

std::string ConsoleOutput::get_time()
{
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(30, '\0');
    size_t size =  std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    s.resize(size);
    return s;
}

void ConsoleOutput::WriteTimeStep(const RunTimeControl& RTC, const string Message,
        size_t ColumnWidth)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        WriteLine("=");
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Time step" << ": " << to_string(RTC.TimeStep) + "/" + to_string(RTC.MaxTimeStep) << "\n";
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Simulation time" << ": " << RTC.SimulationTime << "\n";
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Wall clock time" << ": " << get_time() << "\n";
        if (Message != "")
        {
            WriteLine("-");
            cout << Message;
        }
        WriteLine("=");
    }
}

void ConsoleOutput::WriteTimeStep(const int tStep, const int nSteps, const string Message,
        size_t ColumnWidth)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        WriteLine("_");
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "TimeStep" << " " << to_string(tStep) + "/" + to_string(nSteps) << "\n";
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Time" << " " << get_time() << "\n";
        if (Message != "")
        {
            cout << Message << "\n";
        }
        WriteLine("_");
    }
}

void ConsoleOutput::WriteTimeStep(const int tScreenWrite, const int tStep,
        const int nSteps, const string Message, size_t ColumnWidth)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        if (!(tStep%tScreenWrite))
        {
            WriteTimeStep(tStep, nSteps, Message, ColumnWidth);
        }
    }
}

void ConsoleOutput::WriteLine(const string LineType)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        cout << setfill(LineType[0]) << setw(LineLength) << "" << "\n";
    }
}

void ConsoleOutput::WriteLineInsert(const string Insert, const string LineType)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        cout << LineType[0] << LineType[0] << "< " << Insert << " >" << setfill(LineType[0]) <<
                setw(std::max<int>(0,LineLength-Insert.size()-6)) << "" << "\n";
    }
}

void ConsoleOutput::WriteBlankLine(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        cout << "\n";
    }
}

void ConsoleOutput::WriteSimple(const string Message)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        cout << Message << "\n";
    }
}

void ConsoleOutput::WriteCoordinate(const int x, const int y, const int z, const double dx)
{
    WriteStandard("Point",iVector3{x,y,x});
    WriteStandard("Coordinate",dVector3{x*dx,y*dx,x*dx});
}

void ConsoleOutput::WriteWithinMethod(const string Message, const string Instance,
                                      const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        string thisInstance = Instance;
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        cout << setfill(' ') << setw(10)  << left << thisInstance << "\n"
                                                  << Message      << "\n";
        WriteLine("-");
    }
}

void ConsoleOutput::WriteWarning(const string Message, const string Instance,
                        const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Warning)
    {
        string thisInstance = Instance;
        cerr << setfill('~') << setw(LineLength) << "" << "\n";
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        cerr << setfill(' ') << setw(10)  << left << "Warning: "  << thisInstance << "\n"
                                                  << "          " << Message      << "\n";
        cerr << setfill('~') << setw(LineLength) << "" << "\n";
    }
}

void ConsoleOutput::WriteExit(const string Message, const string Instance, const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        string thisInstance = Instance;
        cerr << "\n";
        cerr << setfill('*') << setw(LineLength) << "" << "\n";
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        cerr << setfill(' ') << setw(10) << left << "Calculation terminated!" << "\n"
                             << setw(10) << left << "Instance:" << thisInstance << "\n"
                             << setw(10) << left << "Reason: " << Message << "\n"
                             << setw(10) << left << "Time: " << get_time() << "\n";
        cerr << setfill('*') << setw(LineLength) << "" << "\n";
        cerr << "\n";
    }
}

void ConsoleOutput::WriteStartScreen(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        cout << "\n";
        WriteLine(">");
        cout << "  OpenPhase\n\n"
             << "  Copyright (c) Ruhr-Universitaet Bochum, Universitaetsstrasse 150, 44801 Bochum, Germany\n"
             << "            and OpenPhase Solutions GmbH, Universitaetsstrasse 136, 44799 Bochum, Germany.\n"
             << "  All rights reserved.\n";
        WriteLine("<");
        cout << "\n";
        cout << "Build time: " << BUILD_TIME << "\n";
        cout << "Git commit SHA: " << GIT_COMMIT_SHA << "\n";
        cout << endl;
    }
}

void ConsoleOutput::PressEnterToContinue(void)
{
    cout << "Press ENTER to continue... " << flush;
    cin.ignore( numeric_limits <streamsize> ::max(), '\n' );
}

void ConsoleOutput::StartProgressIndicator(std::string IndicatorTitle)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        ConsoleOutput::WriteSimple(IndicatorTitle+":");
        cout << "0------------------25------------------50------------------75---------------100%" << endl;
      //cout << "0__________________25__________________50__________________75_______________100%" << endl;
    }
}

void ConsoleOutput::AdvanceProgressIndicator(double start_value, double end_value, double current_value, int& pos)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        int current_value_normalized = 80.0*fabs((current_value - start_value)/(end_value - start_value));
        for(int n = pos; n < current_value_normalized; n++)
        {
            pos++;
            cout << "#";
        }
        cout << flush;
    }
}

void ConsoleOutput::EndProgressIndicator(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Normal)
    {
        cout << endl;
        ConsoleOutput::WriteSimple("Done!");
    }
}

}// namespace openphase
