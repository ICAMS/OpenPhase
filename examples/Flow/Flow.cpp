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
 *   Main contributors :   Oleg Shchyglo; Marvin Tegeler; Raphael Schiedung
 *
 */

#include "AdvectionHR.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Temperature.h"
#include "TextOutput.h"
#include "Tools/TimeInfo.h"
#include "Velocities.h"

using namespace openphase;

int main(int argc, char *argv[])
{

#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_FUNNELED, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    std::string InputFile = DefaultInputFileName;
    if (argc > 1) InputFile = argv[1];

    Settings                            OPSettings(InputFile);
    RunTimeControl                      RTC(OPSettings,InputFile);
    PhaseField                          Phi(OPSettings);
    BoundaryConditions                  BC(OPSettings,InputFile);
    Temperature                         Tx(OPSettings,InputFile);
    Composition                         Cx(OPSettings,InputFile);
    InterfaceProperties                 IP(OPSettings,InputFile);
    DoubleObstacle                      DO(OPSettings);
    DrivingForce                        dG(OPSettings,InputFile);
    EquilibriumPartitionDiffusionBinary DF(OPSettings,InputFile);
    FlowSolverLBM                       FL(OPSettings, RTC.dt, InputFile);
    Velocities                          Vel(OPSettings);
    AdvectionHR                         ADHR (OPSettings);
    TimeInfo                            Timer;

    int x1 = OPSettings.Grid.Nx/2;
    int y1 = OPSettings.Grid.Ny/2;
    int z1 = OPSettings.Grid.Nz/2;

    if(RTC.Restart)
    {
        std::cout << "Restart data being read!\n";
        Phi.Read(OPSettings, BC, RTC.tStart);
        Cx .Read(OPSettings, BC, RTC.tStart);
        Tx .Read(OPSettings, BC, RTC.tStart);
        FL .Read(OPSettings, BC, RTC.tStart);
        RTC.tStart += 1;
        std::cout << "Done\n";
    }
    else
    {
        std::cout << "Set initial state\n";
        // Single liquid phase starting configuration
        size_t idx1 = Initializations::Single(Phi, 0, BC);
        // Adding the nucleus in the middle of the simulation domain
        size_t idx2 = Initializations::Sphere(Phi, 1, 5.0, x1, y1, z1, BC);
         // Orienting the nucleus as desired
        //EulerAngles locAngle({Pi/2.0, 0.0*Pi/8.0, Pi/8.0}, XYZ);
        EulerAngles locAngle({5.0*Pi/8.0, 3.0*Pi/8.0, Pi/8.0}, XYZ);
        Phi.FieldsProperties[idx2].Orientation = locAngle.getQuaternion();
        Phi.FieldsProperties[idx1].Mobile  = true;
        Phi.FieldsProperties[idx2].Mobile  = true;

        // Liquid channel surrounded by alpha phase
        /*std::vector<size_t> ids = Initializations::TwoWalls(Phi, 0, 1, 10,BC,OPSettings);*/

        FL.SetUniformVelocity(BC, dVector3::ZeroVector());
        Cx.SetInitialMoleFractions(Phi);
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        std::cout << "Done\n";

        std::cout << "Write initial state to file\n";
        //  Output to file
        {
            Cx  .WriteVTK(OPSettings, 0);
            Tx  .WriteVTK(OPSettings, 0);
            Vel .WriteVTK(OPSettings, 0);
            FL  .WriteVTK(OPSettings, Phi, 0);
        }

        // Write restart output
        {
            Phi.Write(OPSettings, 0);
            Cx .Write(OPSettings, 0);
            Tx .Write(OPSettings, 0);
            FL .Write(OPSettings, 0);
        }
        RTC.tStart = 1;
        std::cout << "Done\n";
    }

    std::cout << "Entering the Time Loop!!!\n";
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();

        dG.Clear();
        Timer.SetTimeStamp("Clear Driving Force");

        IP.Set(Phi, Tx, BC);
        Timer.SetTimeStamp("Interface properties");

        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        Timer.SetTimeStamp("DO.CalculatePhaseFieldIncrements");

        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("DF.GetDrivingForce");

        dG.Average(Phi, BC);
        Timer.SetTimeStamp("dG.Average");

        dG.MergePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("dG.MergePhaseFieldIncrements");

        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Phi.NormalizeIncrements");

        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
        Timer.SetTimeStamp("DF.Solve");

        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
        Timer.SetTimeStamp("Tx.Set");

        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Phi.MergeIncrements");

        FL.Solve(Phi, Cx, Vel, BC);
        Timer.SetTimeStamp("FL.Solve");

        //InteractionSolidFluid::CalculateSolidVelocities(Phi, Vel, BC, OPSettings, RTC.dt);
        //ADHR.AdvectPhaseField(Phi, Vel, BC, OPSettings.dx, RTC.dt, RTC.tStep);
        //Timer.SetTimeStamp("ADHR.AdvectPhaseField Phi");
        
        Cx.Advect(ADHR,Vel,Phi,BC,RTC.dt,RTC.tStep);
        //DF.CalculatePhaseConcentrations(Phi, Cx, Tx);
        Timer.SetTimeStamp("ADHR.AdvectField Cx");

        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx .WriteVTK(OPSettings, RTC.tStep);
            Tx .WriteVTK(OPSettings, RTC.tStep);
            Vel.WriteVTK(OPSettings, RTC.tStep);
            FL .WriteVTK(OPSettings, Phi, RTC.tStep);
            //Cx.WriteStatistics(tStep, OPSettings.dt);
        }

        if(RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
            Cx .Write(OPSettings, RTC.tStep);
            Tx .Write(OPSettings, RTC.tStep);
            FL .Write(OPSettings, RTC.tStep);
        }

        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message;
            message  = ConsoleOutput::GetStandard("Interface energy density", I_En);
            message += ConsoleOutput::GetStandard("Interface energy",         DO.Energy(Phi, IP));
            message += ConsoleOutput::GetStandard("Center of Mass position",  Phi.FieldsProperties[1].Rcm);
            message += ConsoleOutput::GetStandard("Center of Mass velocity",  Phi.FieldsProperties[1].Vcm);
            message += ConsoleOutput::GetStandard("Acceleration",             Phi.FieldsProperties[1].Acm);
            message += ConsoleOutput::GetStandard("Force",                    Phi.FieldsProperties[1].Force);
            message += ConsoleOutput::GetStandard("Torque",                   Phi.FieldsProperties[1].Torque);
            message += ConsoleOutput::GetStandard("Moment of inertia",        Phi.FieldsProperties[1].InertiaM);
            //message += Info::GetStandard("Orientation",             Phi.FieldsProperties[1].Orientation);
            ConsoleOutput::WriteTimeStep(RTC, message);

            //  Statistics
            Phi.PrintPointStatistics(x1,y1,z1);
            Cx .PrintPointStatistics(x1,y1,z1);
            Tx .PrintPointStatistics(x1,y1,z1);
            Vel.PrintPointStatistics(x1,y1,z1);
            Vel.PrintPointStatistics(x1,y1, 0);
            dG .PrintDiagnostics();
            Phi.PrintPFVolumes();
            Timer.PrintWallClockSummary();
        }
    } //end time loop

#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif

    return 0;
}
