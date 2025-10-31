/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2017
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
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "ConsoleOutput.h"
#include "Initializations.h"
#include "LocalLBM.h"
#include "Macros.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"

#include <complex>

using namespace openphase;

int main(int argc, char **argv)
{
    std::string InputFileName = DefaultInputFileName;
    double ErrorGap   = DBL_MAX;
    double ErrorForce = DBL_MAX;

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

    AdvectionHR           ADHR  (OPSettings);
    BoundaryConditions    BC    (OPSettings);
    InteractionSolidFluid ISF   (OPSettings);
    PhaseField            Phase (OPSettings);
    RunTimeControl        RTC   (OPSettings);
    Velocities            Vel   (OPSettings);

    LocalLBM              LB    (OPSettings,RTC.dt);

    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("Lattice Boltzmann Capillary Bridge Benchmark");
    ConsoleOutput::WriteStandard("Source", DefaultInputFileName);

    std::fstream inpF(InputFileName, std::ios::in);
    if (!inpF)
    {
        std::stringstream message;
        message << "File \"" << InputFileName << "\" could not be opened";
        throw std::runtime_error(message.str());
    };
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    int moduleLocation = FileInterface::FindModuleLocation(inp, "LBCapillaryBridge");

          int    tAdvect          = FileInterface::ReadParameterI(inp, moduleLocation, std::string("iTAdvect"));
    const bool   AdvectSolids     = FileInterface::ReadParameterB(inp, moduleLocation, std::string("bAdvect"));
    const bool   EnforceVolume    = FileInterface::ReadParameterB(inp, moduleLocation, std::string("bEnforceV"));
    const bool   EnforceDiameter  = FileInterface::ReadParameterB(inp, moduleLocation, std::string("bEnforceD"));
    const bool   EnforceCurvature = FileInterface::ReadParameterB(inp, moduleLocation, std::string("bEnforceC"));

    // Define short hand for system size
    const size_t TotalNx = OPSettings.Grid.TotalNx;
    const size_t Nx = OPSettings.Grid.Nx;
    const size_t Ny = OPSettings.Grid.Ny;
    const size_t Nz = OPSettings.Grid.Nz;

    // Determine the number of time-steps the system needs to reach the equilibrium
    const double dx = OPSettings.Grid.dx;
    const double dt = RTC.dt;

    // Determine the liquid radius in lattice units
    const double lbR0 = TotalNx/3;
    const double R0   = TotalNx/3*dx;

    // Correct extensive quantities in 2D
    const bool ThreeDim = Nx > 1 and Ny > 1 and Nz > 1;

    // Define phase indices
    size_t FluidPhaseIdx = 0;
    size_t WallPhaseIdx  = 1;

    // Define phase-field indices
    //size_t FluidIdx = 0;
    size_t WallIdx  = 1;

    if(RTC.Restart)
    {
        LB   .Read (OPSettings, BC, RTC.tStart);
        Phase.Read (OPSettings, BC, RTC.tStart);
        Vel  .Read (OPSettings, BC, RTC.tStart);
    }
    else
    {
        // Set start time step
        RTC.tStart = 0;

        // Determine dimension of plate
        const double Lx = 4.00*lbR0 + 1.0*Phase.Grid.iWidth;
        const double Ly = 1.00*lbR0 + 1.0*Phase.Grid.iWidth;
        const double Lz = 4.00*lbR0 + 1.0*Phase.Grid.iWidth; // Thickness of the walls
        const double j0 = 1.50*lbR0;

        // Initialize phase fields
        Initializations::Single(Phase, FluidPhaseIdx, BC);
        WallIdx  = Initializations::Rectangular(Phase, WallPhaseIdx, Lx , Ly, Lz, 0, j0, 0, BC);

        // Initialize liquid bridge of Radius R
        STORAGE_LOOP_BEGIN(i,j,k,LB.DensityWetting,0)
        {
            if (j >= 0 and j <= j0)
            {
                const dVector3 pos  ({double(i+OPSettings.Grid.OffsetX),double(j+OPSettings.Grid.OffsetY),double(k+OPSettings.Grid.OffsetZ)});
                const double r = std::sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
                LB.DensityWetting(i,j,k,{0}) = LB.DensityProfile(r-lbR0);
            }
            else LB.DensityWetting(i,j,k,{0}) = LB.VaporDensity[0];

            LB.MomentumDensity(i,j,k,{0}) = {0.,0.,0.};
        }
        STORAGE_LOOP_END
        LB.FinalizeInitialiation(Phase, Vel, BC);

        ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);
    }

    // Open log files
    std::string FileName = OPSettings.TextDir + "TimeLog.csv";
    std::fstream log(FileName, std::ios::trunc | std::ios::out);
    log << std::scientific << std::setprecision(16);

    double volumeW = Phase.FieldsProperties[WallIdx].Volume*dx*dx*dx;
    if (ThreeDim) volumeW *= 4; else volumeW *= 2*lbR0;
    const double SolidMass = Phase.FieldsProperties[WallIdx].Density*volumeW;
    double tc = std::sqrt(SolidMass/LB.SurfaceTension[0]);
    assert(tc > 0.0);

    // Start time loop
    std::vector<double> WallForces;
    ConsoleOutput::WriteSimple("Entering the Time Loop!!!");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        // NOTE: In the planar case
        volumeW = Phase.FieldsProperties[WallIdx].Volume*dx*dx*dx;
        if (ThreeDim) volumeW *= 4; else volumeW *= 2*lbR0;

        // Set density of solid wall to ensure mass conservation
        Phase.FieldsProperties[WallIdx].Density = (volumeW!=0)? SolidMass/volumeW : 0.0;

        if (RTC.WriteRawData())
        {
            LB   .Write(OPSettings, RTC.tStep);
            Phase.Write(OPSettings, RTC.tStep);
            Vel  .Write(OPSettings, RTC.tStep);
        }

        if(RTC.WriteVTK())
        {
            LB   .WriteVTK(OPSettings, Phase, RTC.tStep);
            Phase.WriteVTK(OPSettings, RTC.tStep);
        }

        // Calculate mid radius of the capillary bridge
        const Grain Wall  = Phase.FieldsProperties[WallIdx];
        WallForces.push_back(Wall.Force[1]);

        //  Output to screen and do some diagnostics
        if(RTC.WriteToScreen())
        {
            const double t            = RTC.tStep*dt;
            const double t0           = tAdvect*dt;
            const double Time         = t-t0;
            const double TimeRescaled = Time/tc;

            // Calculate interface tensions (via density profile)
            const double sigma_lv = LB.CalculateSigmaX(0, 0);
            const double sigma_sv = LB.CalculateSigmaY(Wall.Rcm[1], Ny-1, 0);
            const double sigma_sl = LB.CalculateSigmaY(0, Wall.Rcm[1], 0);
            const double kappa    = LB.AverageCurvature(0, 0.25*lbR0, 0, Nz);
            const double theta    = LB.CalculateContactAngle(sigma_sl, sigma_sv, sigma_lv);

            // TODO Calculate pressures
            // const double pL       = LB.Pressure(   0,    0,    0);
            // const double pV       = LB.Pressure(Nx-1, Ny-1, Nz-1);
            // const double dp       = pL - pV;

            // TODO Calculate interface tensions (via curvature and pressure)
            //const double sigma_laplace = (std::abs(kappa) > DBL_EPSILON) ? dp/kappa : 0;

            // Calculate volumes
            std::array<double,2> lbVolumes = LB.CalculateFluidVolumes();
            const double volumeL  = lbVolumes[0]*dx*dx*dx;
            const double volumeV  = lbVolumes[1]*dx*dx*dx;

            // Calculate masses
            const double massF = LB.CalculateFluidMass()[0];

            // Calculate densities
#ifdef MPI_PARALLEL
            double locDensityVapor  = (MPI_RANK==MPI_SIZE-1) ? LB.Density(Phase, Phase.Grid.Nx, Ny, (Nz == 1)?0:1) : 0;
            double locDensityLiquid = (MPI_RANK==0) ? LB.Density(Phase,0,0,0) : 0;
            double DensityVapor  = 0;
            double DensityLiquid = 0;
            OP_MPI_Allreduce(&locDensityVapor , &DensityVapor , 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(&locDensityLiquid, &DensityLiquid, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#else
            double DensityVapor  = LB.Density(Phase, Nx, Ny, (Nz == 1)?0:1);
            double DensityLiquid = LB.Density(Phase,0,0,0);
#endif

            // Sort forces in ascending order to achieve better numerical accuracy
            std::sort(WallForces.begin(), WallForces.end(), [](double i, double j) {return std::abs(i) < std::abs(j);});
            const size_t N = WallForces.size() > 0 ? WallForces.size() : 1;

            // Calculate average Force on wall
            long double AverageWallForce = 0.0;
            for (auto WallForce: WallForces) AverageWallForce += WallForce;
            AverageWallForce /= N;

            // Calculate standard deviation of Force on wall
            long double stdDevWallForce = 0.0;
            for (auto WallForce: WallForces)
            {
                stdDevWallForce += std::pow(WallForce - AverageWallForce,2);
            }
            stdDevWallForce /= N;
            stdDevWallForce = std::sqrt(stdDevWallForce);
            WallForces.clear();

            const double RMid       = LB.CalculateRadius(OPSettings, 0, 0);
            const double WallForce  = Wall.Force[1];
            const double F_ex       = ThreeDim ? -0.25*Pi*RMid*sigma_lv : -sigma_lv*dx;
            if (TimeRescaled < 0)
            if (std::abs(F_ex) > 100*DBL_EPSILON)
            {
                 ErrorForce = std::abs(AverageWallForce-F_ex)/std::abs(F_ex);
            }

            //Calculate Wall velocity 
            const double WallVelocity = Phase.FieldsProperties[WallIdx].Vcm[1];
            const double WallVelocityExpected = (Time > 0.0) ? 2.0*LB.SurfaceTension[0]*dx/(SolidMass/lbR0)*Time : 0.0;

            // Calculate wall distance between the walls
            //const double WallSpeed         = Wall.Vcm[1];
            const double WallGap           = (2.0*Wall.Rcm[1]*dx - R0);
            const double WallGapExpected   = (Time > 0.0) ? LB.SurfaceTension[0]*dx/(SolidMass/lbR0)*Time*Time : WallGap;
            const double WallGapRescaled   = (2.0*Wall.Rcm[1]*dx - R0)/R0;
            double WallGapRescaledEx = 2.0;
            if (TimeRescaled>0)
            {
                if (ThreeDim) WallGapRescaledEx = 2.0 -Pi*TimeRescaled*TimeRescaled;
                else          WallGapRescaledEx = 2.0 - 2*TimeRescaled*TimeRescaled;
            }

            if (TimeRescaled>0)
            {
                ErrorGap  = std::abs(WallGapRescaled-WallGapRescaledEx);
                ErrorGap /= WallGapRescaledEx;
            }

            std::array<std::stringstream,2> line;
            line[1] << std::scientific << std::setprecision(16);
            ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Time [s]", Time);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Time Rescaled", TimeRescaled);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Volume of Vapor[m^3]", volumeV);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Volume of Liquid [m^3]", volumeL);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Volume of Wall [m^3]", volumeW);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Mass of Fluid [kg]", massF);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Mass of Wall [kg]", SolidMass);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Desity of Vapor [kg/m^3]", DensityVapor);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Desity of Liquid [kg/m^3]", DensityLiquid);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Desity of Wall [kg/m^3]", Phase.FieldsProperties[WallIdx].Density);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Surface Tension (integral) [kg/s^2]", sigma_lv);
            //Info::WriteWigthLog("Surface tension (dp/kappa) [kg/s^2]", sigma_laplace);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Expected Surface Tension [kg/s^2]", LB.SurfaceTension[0]);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Contact Angle [deg]", theta);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Wetting Parameter [%]", LB.Wetting[FluidPhaseIdx][WallPhaseIdx]);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Capillary Diameter [m]", 2*RMid);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Capillary Diameter Rescaled", 2*RMid/R0);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Curvature [1/m]", kappa);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Curvature Rescaled", kappa*R0);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Force on Wall [kg*m/s^2]", WallForce);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Average Force on Wall [kg*m/s^2]", AverageWallForce);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Expected Force on Wall [kg*m/s^2]", F_ex);
            ConsoleOutput::Write("-");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Relative Force Error", ErrorForce);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Wall Velocity [m/s]", WallVelocity);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Expected Wall Velocity [m/s]", WallVelocityExpected);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Distance Between Walls [m]", WallGap);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Expected Distance Between Walls [m]", WallGapExpected);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Rescaled Distance Between Walls", WallGapRescaled);
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Expected Rescaled Distance between Walls", WallGapRescaledEx);
            ConsoleOutput::Write("-");
            ConsoleOutput::WriteWithLog(line, RTC.tStep, "Relative Distance Error", ErrorGap);
            ConsoleOutput::Write("");
            ConsoleOutput::WriteLineToLogfile(log, line, RTC.tStep);

            if (RTC.tStep < tAdvect and sigma_lv > DBL_EPSILON) tc = std::sqrt(SolidMass/sigma_lv);
            if (WallGapRescaled < 1.25) break;
        }

        Phase.ClearGrainsForcesAndAccelerations();
        LB.Solve(Phase, Vel, BC);

        // Let walls move
        if (AdvectSolids and RTC.tStep > tAdvect)
        {
            // Allow only uni-axial motion
            Phase.FieldsProperties[WallIdx].Force.setX(0);
            Phase.FieldsProperties[WallIdx].Force.setZ(0);

            // Calculate wall velocities
            ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);

            // Move solid walls
            Phase.Advect(ADHR, Vel, BC, LB, RTC.dt, RTC.tStep);
        }
        else if (RTC.tStep > R0*R0/dx/dx)
        {

            if (EnforceDiameter)
            {
                const double DMid = 2*LB.CalculateRadius(OPSettings, 0, 0);
                LB.EnforceLiquidDiameter(BC, DMid, 2*R0);
            }

            if (EnforceVolume)
            {
                std::array<double,2> volume = LB.CalculateFluidVolumes();
                const double lbVCurrent = volume[0];
                const double lbVGoal    = ThreeDim ? lbR0*lbR0*lbR0 : lbR0*lbR0;
                LB.EnforceVolume(BC, lbVCurrent, lbVGoal);
            }

            if (EnforceCurvature)
            {
                const double kappa  = LB.AverageCurvature(0, 0.25*lbR0, 0, Nz);
                LB.EnforceCurvature(LB.Wetting[FluidPhaseIdx][WallPhaseIdx], kappa, ThreeDim ? 1.0/R0 : 0.0);
            }
        }
    }// End of time loop

    std::ofstream resultsSim ("Results.sim");
    if (resultsSim.is_open())
    {
        resultsSim << 2 << "\n";
        resultsSim << "ErrorWallForce "    << ErrorForce << " " << ErrorForce << "\n";
        resultsSim << "ErrorWallDistance " << ErrorGap   << " " << ErrorGap;
        resultsSim.flush();
        resultsSim.close();
    }

#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif

    if (ErrorGap < 0.4493)
    {
        ConsoleOutput::Write("BENCHMARK SUCCESSFUL COMPLETED\n");
        return EXIT_SUCCESS;
    }
    else
    {
        ConsoleOutput::Write("BENCHMARK FAILED\n");
        return EXIT_FAILURE;
    }
}

