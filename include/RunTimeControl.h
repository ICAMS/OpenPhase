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

#ifndef RUNTIMECONTROL_H
#define RUNTIMECONTROL_H

#include "Includes.h"
//#include "Macros.h"
//#include "Definitions.h"

#include <cstdio>
#include <fstream>
#include <string>

namespace openphase
{
class Settings;

class OP_EXPORTS RunTimeControl                                                 ///< Run time control module.
{
 public:
    std::string thisclassname;                                                  ///< Object's implementation class name
    std::string thisobjectname;                                                 ///< Object's name

    std::string UnitsOfLength;                                                  ///< Units of length
    std::string UnitsOfMass;                                                    ///< Units of mass
    std::string UnitsOfTime;                                                    ///< Units of time
    std::string UnitsOfEnergy;                                                  ///< Units of energy

    double dt;                                                                  ///< Initial time step (s)
    double SimulationTime;                                                      ///< Current simulation time: \Sum_{tStep} dt(tStep)
    int TimeStep;                                                               ///< Current time step index
    int& tStep = TimeStep;
    int MaxTimeStep;                                                            ///< Maximum number of time steps in simulation
    int& nSteps = MaxTimeStep;
    int VTKOutputInterval;                                                      ///< Visualization files writing frequency (in time steps)
    int& tFileWrite = VTKOutputInterval;
    int ConsoleOutputInterval;                                                  ///< Screen writing frequency (in time steps)
    int& tScreenWrite = ConsoleOutputInterval;
    int StartTimeStep;                                                          ///< Starting time step index. Used for setting the restart time step
    int& tStart = StartTimeStep;
    bool RestartSwitch;                                                         ///< Restart switch ("No" if restart is OFF, "Yes" if ON)
    bool& Restart = RestartSwitch;
    int CheckpointInterval;                                                     ///< Restart data writing frequency (in time steps)
    int& tRestartWrite = CheckpointInterval;
    int OpenMPThreads;                                                          ///< Number of OpenMP threads for parallel execution
    bool StopTrigger;                                                           ///< True if stop condition has been triggered
    myclock_t CheckTime;                                                        ///< Time when the program has last checked if a stop has been triggered
    double StopCheckInterval;                                                   ///< Interval between stop checks

    bool LogModeScreen;                                                         ///< Logarithmic writing to Screen
    bool LogModeVTK;                                                            ///< Logarithmic writing to VTK
    int LogScreenFactor;                                                        ///< Outputs to screen per decade = 10 x LogScreenFactor
    int LogVTKFactor;                                                           ///< Outputs to VTK    per decade = 10 x LogVTKFactor

    std::string SimulationTitle;                                                ///< Simulation title string
    std::string VTKDir;                                                         ///< Directory name for the VTK files
    std::string RawDataDir;                                                     ///< Directory name for the raw data files
    std::string TextDir;                                                        ///< Directory name for the text files

    RunTimeControl(){};                                                         ///< Default constructor
    RunTimeControl(Settings& locSettings,
                   const std::string InputFileName = DefaultInputFileName)      ///< Constructor
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "");  ///< Initializes run time control object
    void ReadInput(const std::string InputFileName);                            ///< Reads run time control parameters
    void ReadInput(std::stringstream& inp);                                     ///< Reads run time control parameters
    void ReadJSON(const std::string InputFileName);                                       ///< Reads run time control parameters
    void Write();                                                               ///< Write tStep to file
    bool RestartPossible();                                                     ///< Reads tStep from file returns true if successful
    bool WriteToScreen()                                                        ///< Returns true at regular console output intervals
    {
        if (LogModeScreen)
        {
            if (TimeStep > 1)
            {
                int tWrite = std::pow(10,std::floor(std::log10(TimeStep)))/LogScreenFactor;
                if (tWrite == 0) tWrite = 1;
                return !(TimeStep%tWrite);
            }
            return true;
        }
        return (!(TimeStep%ConsoleOutputInterval) or StopTrigger);
    }
    bool WriteVTK()                                                             ///< Returns true at regular VTK files writing intervals
    {
        if (LogModeVTK)
        {
            if (TimeStep > 1)
            {
                int tWrite = std::pow(10,std::floor(std::log10(TimeStep)))/LogVTKFactor;
                if (tWrite == 0) tWrite = 1;
                return !(TimeStep%tWrite);
            }
            return true;
        }
        return (!(TimeStep%VTKOutputInterval) or StopTrigger);
    }
    bool WriteRawData()                                                         ///< Returns true at regular raw data files writing intervals
    {
        return (!(TimeStep%CheckpointInterval) or StopTrigger);
    }
    void IncrementTimeStep()                                                    ///< Increments simultaneously the time step and simulation time
    {
        #ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
        {
        #endif
            StopTrigger = CheckStop();
            if (StopTrigger) MaxTimeStep = TimeStep + 1;
        #ifdef MPI_PARALLEL
        }
        OP_MPI_Bcast(&nSteps, 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
        int iStop = StopTrigger;
        OP_MPI_Bcast(&iStop, 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
        StopTrigger = iStop;
        #endif
        SimulationTime += dt;
        TimeStep++;
    }
    void SetNewTimeStep(double new_dt)                                          ///< Sets new time step
    {
        dt = new_dt;
    }
    bool CheckStop()
    {
        double time = mygettime();
        double diff = time - CheckTime;
        if (diff < StopCheckInterval)
        {
            return false;
        }
        bool stop = false;
        CheckTime = time;
        std::string filename;
        filename = TextDir + "stop";
        std::ifstream file(filename);
        if(file.is_open())
        {
            stop = true;
            file.close();
            std::remove(filename.c_str());
        }
        return stop;
    }
    RunTimeControl& operator= (const RunTimeControl& rhs);                      ///< Copy operator
};

}// namespace openphase

#endif
