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

#include "AdvectionHR.h"
#include "Globals.h"
#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "GrandPotential/GrandPotentialDensity.h"
#include "GrandPotential/GrandPotentialSolver.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Temperature.h"
#include "Tools/AnalysisSintering.h"
#include "Tools/TimeInfo.h"
#include "Velocities.h"
#include "DrivingForce.h"

#include <iostream>

namespace op = openphase;

#ifdef MPI_PARALLEL
   int MPI_RANK;
   int MPI_SIZE;
#endif

void ReadInitializationInuptParameters(op::Settings& locSettings,
        std::string InputFileName, double& RelativeDensity, double& MeanRadius,
        double& StdRadius, long& SX, long& SY, long& SZ, size_t& RandomSeed);

int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_SINGLE, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
    std::cout << "MPI_RANK: " << MPI_RANK << "/" << MPI_SIZE << "\n";
#endif

    std::string InputFileName = op::DefaultInputFileName;
    if (argc > 1) InputFileName = argv[1];

    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID  );  // domain error occurred

    op::Settings OPSettings(InputFileName);

    op::AdvectionHR              ADHR  (OPSettings, InputFileName);
    op::BoundaryConditions       BC    (OPSettings, InputFileName);
    op::DoubleObstacle           DO    (OPSettings, InputFileName);
    op::GrandPotentialDensity    omega (OPSettings, InputFileName);
    op::GrandPotentialSolver     GPS   (OPSettings, InputFileName);
    op::InteractionSolidSolid    ISS   (OPSettings, InputFileName);
    op::InteractionSolidFluid    ISF   (OPSettings, InputFileName);
    op::InterfaceProperties      IP    (OPSettings, InputFileName);
    op::Temperature              Temp  (OPSettings, InputFileName);
    op::PhaseField               Phase (OPSettings, InputFileName);
    op::RunTimeControl           RTC   (OPSettings, InputFileName);
    op::AnalysisSintering        AS    (OPSettings, InputFileName);
    op::DrivingForce             dG    (OPSettings, InputFileName);
    op::TimeInfo                 Timer (OPSettings, "Execution Time Statistics");
    op::Velocities               Vel   (OPSettings);

    omega.Set(Temp.Tx,GPS.ChemicalPotential);

    if (RTC.RestartPossible()) // Restart Simulation
    {
        Phase.Read(OPSettings, BC, -1);
        Temp.SetInitial(BC);//Temp.Read(OPSettings, BC, -1); TODO 
        if (RTC.tStart != 0) GPS.Read(BC, Phase, omega);
        else GPS.SetInitial(Phase, omega, BC);
    }
    else // Initialize Simulation
    {
        RTC.tStart = 0;
        double RelativeDensity, MeanRadius, StdRadius;
        long SizeX, SizeY, SizeZ;
        size_t RandomSeed;
        ReadInitializationInuptParameters(OPSettings, InputFileName,
                RelativeDensity,MeanRadius, StdRadius, SizeX, SizeY, SizeZ,
                RandomSeed);
        auto SolidPhaseIndex = [](double i, double j, double k)
        {
            return 1.0;
            //TODO Fix code below
            //int idx = 1;
            //double step = (XN-X0)/(OPSettings.Nphases-1);
            //for (size_t n = 2; n < OPSettings.Nphases; n++)
            //    if (i >= X0 + (n-1)*step) idx++;
            //return idx;
        };

        op::Initializations::RectangularGreenBody(Phase, BC, OPSettings, SolidPhaseIndex,
                0, SizeX, SizeY, SizeZ, MeanRadius, StdRadius,
                RandomSeed, RelativeDensity);
        Temp.SetInitial(BC);
        GPS.SetInitial(Phase, omega, BC);
    }

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteToScreen())
        {
            op::ConsoleOutput::WriteTimeStep(RTC);
            op::ConsoleOutput::Write("Temperature", Temp.Tx(0,0,0));
            AS.DoDiagnostics(OPSettings, Phase, DO, IP, omega, GPS, RTC);
            Phase.FieldsProperties.WriteTableVolumes(RTC.tStep,OPSettings);
            if (RTC.tStep != 0) Timer.PrintWallClockSummary();
        }

        if (RTC.WriteVTK())
        {
            //Vel  .WriteVTK (RTC.tStep, OPSettings);
            GPS  .WriteVTK (OPSettings, Phase, omega, RTC.tStep, 8);
            Phase.WriteVTK (OPSettings, RTC.tStep, false, 8);
            Temp.WriteVTK (OPSettings,RTC.tStep);
        }

        if (RTC.WriteRawData())
        {
            GPS  .Write(OPSettings, -1);
            Phase.Write(OPSettings, -1);
            RTC  .Write();
            //Temp .Write(OPSettings, -1);
        }

        {
            Timer.SetStart();
            Phase.ClearGrainsForcesAndAccelerations();
            auto MassDensity = [&GPS,&Phase,&omega] (int i ,int j, int k){return GPS.MassDensity(i,j,k,Phase,omega);};
            ISS.Calculate(Phase, BC, MassDensity, RTC.dt);                      Timer.SetTimeStamp("Calculate solid-solid interactions");
            ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);               Timer.SetTimeStamp("Calculate solid selocities");
            //Phase.Advect(ADHR, Vel, BC, OPSettings, RTC.dt, RTC.tStep);         Timer.SetTimeStamp("Advect phase-field");
            GPS.Advect(ADHR, Vel, Phase, BC, RTC.dt, RTC.tStep);                Timer.SetTimeStamp("Advection of Mole per volume");
            IP.Set(Phase,BC);                                                   Timer.SetTimeStamp("Calculate interface properties");
            DO.CalculatePhaseFieldIncrements(Phase, IP, dG);                    Timer.SetTimeStamp("DO.CalculatePhaseFieldIncrements");
            GPS.Solve(Phase,omega,BC,IP,Temp,RTC.dt);                           Timer.SetTimeStamp("Solve grand potential method");
            Phase.MergeIncrements(BC,RTC.dt);                                   Timer.SetTimeStamp("Phase merge increments");

            // TODO add condition for checking of time step. I does not need to be
            // checked for each model every time step! 
            double dt_max = GPS.MaximumTimeStep(Phase, Temp, omega);
            if (dt_max < RTC.dt)
            {
                std::stringstream message;
                message << "Time step may be too large \n";
                message << "Current time step [s]: " << RTC.dt << "\n";
                message << "Maximun time step [s]: " << dt_max << "\n";
                op::ConsoleOutput::WriteWarning(message.str(), "SolidPhaseSintering", "main");
            }
        }
    }
    return EXIT_SUCCESS;
}

void ReadInitializationInuptParameters(op::Settings& locSettings,
        std::string InputFileName, double& RelativeDensity, double& MeanRadius,
        double& StdRadius, long& SX, long& SY, long& SZ, size_t& RandomSeed)
{
    // Read initialization input parameters
    op::ConsoleOutput::WriteBlankLine();
    op::ConsoleOutput::WriteLineInsert("ProjectInput");
    op::ConsoleOutput::WriteStandard("Source", InputFileName);
    std::fstream inp(InputFileName);
    if (!inp)
    {
        std::cerr << "File \"" <<  InputFileName << "\" could not be opened\n";
        std::exit(EXIT_FAILURE);
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();
    int moduleLocation = op::FileInterface::FindModuleLocation(inp_data, "ProjectInput");
    RelativeDensity = op::FileInterface::ReadParameterD(inp_data, moduleLocation, "Rho");
    MeanRadius = op::FileInterface::ReadParameterD(inp_data, moduleLocation, "AvR" )/locSettings.Grid.dx;
    StdRadius  = op::FileInterface::ReadParameterD(inp_data, moduleLocation, "StdR")/locSettings.Grid.dx;
    SX = op::FileInterface::ReadParameterD(inp_data, moduleLocation, "SX")/locSettings.Grid.dx;
    SY = op::FileInterface::ReadParameterD(inp_data, moduleLocation, "SY")/locSettings.Grid.dx;
    SZ = op::FileInterface::ReadParameterD(inp_data, moduleLocation, "SZ")/locSettings.Grid.dx;
    std::string SRandomSeed = op::FileInterface::ReadParameterS(inp_data, moduleLocation, "RandomSeed");
    if (SRandomSeed == "random")
    {
        std::random_device rd;
        RandomSeed = rd();
        op::ConsoleOutput::Write("RandomSeed",RandomSeed);
    }
    else RandomSeed = std::stol(SRandomSeed);
}
