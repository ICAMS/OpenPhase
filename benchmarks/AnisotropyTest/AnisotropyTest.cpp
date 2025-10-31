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

#include "BoundaryConditions.h"
#include "ConsoleOutput.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"
#include "UserDrivingForce.h"
#include "Velocities.h"

#ifdef MPI_PARALLEL
#include "mpi_wrapper.h"
#endif

using namespace openphase;

int main(int argc, char **argv)
{
    std::string InputFileName = DefaultInputFileName;

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

    BoundaryConditions    BC    (OPSettings,InputFileName);
    DrivingForce          dG    (OPSettings,InputFileName);
    DoubleObstacle        DO    (OPSettings,InputFileName);
    InterfaceProperties   IP    (OPSettings,InputFileName);
    PhaseField            Phase (OPSettings,InputFileName);
    UserDrivingForce      UDF   (OPSettings,InputFileName);
    RunTimeControl        RTC   (OPSettings,InputFileName);
    TimeInfo              Timer (OPSettings, "Execution Time Statistics");
    Velocities            Vel   (OPSettings);

    //bool Finialize = false;

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
        //const double Lx = OPSettings.Dimensions.Ny/3.0;
        //const double Ly = OPSettings.Dimensions.Ny/3.0;
        //const double Lz = OPSettings.Dimensions.Ny/3.0;

        const double x0 = OPSettings.Grid.Nx/2.0 + OPSettings.Grid.Nx/5.0;
        const double x1 = OPSettings.Grid.Nx/2.0 - OPSettings.Grid.Nx/5.0;
        const double y0 = OPSettings.Grid.Ny/2.0;
        const double z0 = OPSettings.Grid.Nz/2.0;
        const double R  = OPSettings.Grid.Nx/8.0;

        const size_t FluidPhaseIdx = 0;
        const size_t SolidPhaseIdx = 1;


        // Initialize phase fields
        Initializations::Single(Phase, FluidPhaseIdx, BC);
        Initializations::Sphere(Phase, SolidPhaseIdx, R, x0, y0, z0, BC);
        const size_t idxB = Initializations::Sphere(Phase, SolidPhaseIdx, R, x1, y0, z0, BC);
        //Initializations::Rectangular(Phase, SolidPhaseIdx, Lx , Ly, Lz, x0, y0, z0, BC);
        Phase.Finalize(BC);

        const double Q1 =  0 * Pi/180.0;
        const double Q2 =  0 * Pi/180.0;
        const double Q3 = 45 * Pi/180.0;
        EulerAngles locAngles({Q1, Q2, Q3}, XYZ);
        Phase.FieldsProperties[idxB].Orientation = locAngles.getQuaternion();
    }

    ConsoleOutput::WriteSimple("Entering the Time Loop!!!");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteRawData())
        {
            Phase.Write(OPSettings, RTC.tStep);
        }

        if(RTC.WriteVTK())
        {
            Phase.WriteVTK(OPSettings, RTC.tStep, false, 8);
        }

        //  Output to screen and do some diagnostics
        if(RTC.WriteToScreen())
        {
            ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps);
            dG.PrintDiagnostics();
            Timer.PrintWallClockSummary();
        }

        Timer.SetStart();
        dG.Clear();
        IP.Set(Phase, BC);
        Timer.SetTimeStamp("Calculate interface properties");
        DO.CalculatePhaseFieldIncrements(Phase, IP, dG);
        Timer.SetTimeStamp("Calculate curvature contribution");
        UDF.SetDrivingForce(Phase, dG);
        Timer.SetTimeStamp("Calculate user driving force");
        dG.Average(Phase, BC);
        Timer.SetTimeStamp("Average driving force");
        if (RTC.WriteVTK()) dG.WriteVTK(OPSettings, Phase, RTC.tStep);
        dG.MergePhaseFieldIncrements(Phase,IP);
        Timer.SetTimeStamp("Merge driving force increments");
        Phase.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Normalize Phase-field Increments");
        Phase.MergeIncrements(BC,RTC.dt);
        Timer.SetTimeStamp("Merge phase-field increments");
    }
#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif

    return 0;
}

