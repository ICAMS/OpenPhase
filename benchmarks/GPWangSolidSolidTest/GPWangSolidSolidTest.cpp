/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2023
 *
 *   Authors:    Raphael Schiedung
 *
 *   Copyright (c) 2009-2022 Interdisciplinary Centre for Advanced Materials
 *                 Simulation (ICAMS). Ruhr-Universitaet Bochum. Germany
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "AdvectionHR.h"
#include "BoundaryConditions.h"
#include "ConsoleOutput.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "GrandPotential/GrandPotentialDensity.h"
#include "GrandPotential/GrandPotentialSolver.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "Macros.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Temperature.h"
#include "Tools/TimeInfo.h"
#include "Velocities.h"

using namespace openphase;

int main(int argc, char **argv)
{
    std::string InputFileName = DefaultInputFileName;
    if (argc > 1) InputFileName = argv[1];

#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_FUNNELED, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);  // the result was too large
    //feenableexcept(FE_UNDERFLOW); // the result was too small

    Settings           OPSettings;
    OPSettings.ReadInput(InputFileName);

    AdvectionHR           ADHR  (OPSettings,InputFileName);
    BoundaryConditions    BC    (OPSettings,InputFileName);
    DoubleObstacle        DO    (OPSettings,InputFileName);
    DrivingForce          dG    (OPSettings,InputFileName);
    GrandPotentialDensity omega (OPSettings,InputFileName);
    GrandPotentialSolver  GPS   (OPSettings,InputFileName);
    InteractionSolidFluid ISF   (OPSettings,InputFileName);
    InteractionSolidSolid ISS   (OPSettings,InputFileName);
    InterfaceProperties   IP    (OPSettings,InputFileName);
    PhaseField            Phase (OPSettings,InputFileName);
    RunTimeControl        RTC   (OPSettings,InputFileName);
    Temperature           Temp  (OPSettings,InputFileName);
    TimeInfo              Timer (OPSettings, "Execution Time Statistics");
    Velocities            Vel   (OPSettings);

    omega.Set(Temp.Tx,GPS.ChemicalPotential);

    if(RTC.Restart)
    {
        Phase.Read (OPSettings, BC, RTC.tStart);
        Vel  .Read (OPSettings, BC, RTC.tStart);
    }
    else
    {
        // Set start time step
        RTC.tStart = 0;

        // Determine dimension of plate
        const double Lx = OPSettings.Grid.Ny/3.0;
        const double Ly = OPSettings.Grid.Ny/3.0;
        const double Lz = OPSettings.Grid.Ny/3.0;
        const double x0 = OPSettings.Grid.Nx/2.0 - Lx/2.0 - Phase.Grid.iWidth/2.0 + 0.75;
        const double y0 = OPSettings.Grid.Ny/2.0 - Lx/3.0;
        const double z0 = OPSettings.Grid.Nz/2.0;
        const double x1 = OPSettings.Grid.Nx/2.0 + Lx/2.0 + Phase.Grid.iWidth/2.0 - 0.75;
        const double y1 = OPSettings.Grid.Ny/2.0 + Lx/3.0;

        const size_t FluidPhaseIdx = 0;
        const size_t SolidPhaseIdx = 1;

        // Initialize phase fields
        Initializations::Single(Phase, FluidPhaseIdx, BC);
        Initializations::Rectangular(Phase, SolidPhaseIdx, Lx , Ly, Lz, x0, y0, z0, BC);
        Initializations::Rectangular(Phase, SolidPhaseIdx, Lx , Ly, Lz, x1, y1, z0, BC);
        //Initializations::Sphere(Phase, SolidPhaseIdx, Ly/2.0, x0, y0, z0, BC);
        //Initializations::Sphere(Phase, SolidPhaseIdx, Ly/2.0, x1, y0, z0, BC);

        Phase.Finalize(BC);
        GPS.SetInitial(Phase, omega, BC);
        Temp.SetInitial(BC);
    }

    ConsoleOutput::WriteSimple("Entering the Time Loop!!!");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteRawData())
        {
            Phase.Write(OPSettings, RTC.tStep);
            Vel  .Write(OPSettings, RTC.tStep);
        }

        if(RTC.WriteVTK())
        {
            Phase.WriteVTK(OPSettings, RTC.tStep, false, 8);
            GPS  .WriteVTK(OPSettings, Phase, omega, RTC.tStep, 8);
        }

        //  Output to screen and do some diagnostics
        if(RTC.WriteToScreen())
        {
            ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps);
            for (size_t idx = 0u; idx < Phase.FieldsProperties.size(); idx++)
            if (Phase.FieldsProperties[idx].State == AggregateStates::Solid)
            {

                Grain& grain = Phase.FieldsProperties[idx];
                double dV = Phase.Grid.CellVolume(false);
                double grain_mass = grain.Density*grain.Volume*dV;
                double cfl = grain.Vcm.abs()*RTC.dt/Phase.Grid.dx;
                ConsoleOutput::Write("Grain "+std::to_string(idx));
                ConsoleOutput::Write("-");
                ConsoleOutput::Write("Force",grain.Force);
                ConsoleOutput::Write("Mass",grain_mass);
                ConsoleOutput::Write("Center of Mass",grain.Rcm);
                ConsoleOutput::Write("Torque",grain.Torque);
                ConsoleOutput::Write("Moment of Inertia",grain.InertiaM);
                ConsoleOutput::Write("cfl",cfl);
                ConsoleOutput::Write(" ");
            }

            ConsoleOutput::Write("Mole of [mol]", GPS.TotalAmountOfComponent(0));
            ConsoleOutput::Write("Delta Mole of [%]", GPS.TotalAmountOfComponentChange(0)*100);

            Timer.PrintWallClockSummary();
        }

        {
            Timer.SetStart();
            dG.Clear();
            Phase.ClearGrainsForcesAndAccelerations();
            auto MassDensity = [&GPS,&Phase,&omega] (int i ,int j, int k){return GPS.MassDensity(i,j,k,Phase,omega);};
            ISS.Calculate(Phase, BC, MassDensity, RTC.dt);                      Timer.SetTimeStamp("Calculate solid-solid interactions");
            ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);               Timer.SetTimeStamp("Calculate solid selocities");
            //Phase.Advect(ADHR, Vel, BC, OPSettings, RTC.dt, RTC.tStep);       Timer.SetTimeStamp("Advect phase-field");
            GPS.Advect(ADHR, Vel, Phase, BC, RTC.dt, RTC.tStep);                Timer.SetTimeStamp("Advect chemical potential");
            IP.Set(Phase, BC);                                                  Timer.SetTimeStamp("Calculate interface properties");
            DO.CalculatePhaseFieldIncrements(Phase, IP, dG);                    Timer.SetTimeStamp("DO.CalculatePhaseFieldIncrements");
            GPS.Solve(Phase,omega,BC,IP,Temp,RTC.dt);                           Timer.SetTimeStamp("Solve grand potential method");
            Phase.MergeIncrements(BC,RTC.dt);                                   Timer.SetTimeStamp("Phase.MergeIncrements");
        }
    }
#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}

