/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2019
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
#include "DoubleObstacle.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "ConsoleOutput.h"
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"

using namespace openphase;

class LocalLBM : public FlowSolverLBM                                           ///< Modification of FLowSolverLBM which are only used in this example
{
 public:
    LocalLBM(Settings& locSettings, double in_dt): FlowSolverLBM(locSettings, in_dt){};

    void InitializeSingle();                                                    ///<  Initializes a single density
    void InitializeSphere(const double Radius, const int i0, const int j0,
            const int k0, const double* rho = nullptr);  ///<  Initializes a single sphere
};

int main(int argc, char *argv[])
{

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
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small


    Settings OPSettings;
    OPSettings.ReadInput();

    AdvectionHR           ADHR  (OPSettings);
    BoundaryConditions    BC    (OPSettings);
    DoubleObstacle        DO    (OPSettings);
    InteractionSolidFluid ISF   (OPSettings);
    InteractionSolidSolid ISS   (OPSettings);
    InterfaceProperties   IP    (OPSettings);
    PhaseField            Phase (OPSettings);
    RunTimeControl        RTC   (OPSettings);
    Velocities            Vel   (OPSettings);

    LocalLBM              LB    (OPSettings, RTC.dt);

    // Set initial conditions
    if(RTC.Restart)
    {
        Phase.Read(OPSettings, BC, RTC.tStart);
        LB   .Read(OPSettings, BC, RTC.tStart);
        Phase.SetBoundaryConditions(BC);
        LB   .SetBoundaryConditions(BC);
    }
    else
    {
        // Read benchmark specific input parameters
        ConsoleOutput::WriteBlankLine();
        ConsoleOutput::WriteLineInsert("LBGravity");
        ConsoleOutput::WriteStandard("Source", DefaultInputFileName);
        std::fstream inpF(DefaultInputFileName, std::ios::in);
        std::stringstream inp;
        inp << inpF.rdbuf();
        inpF.close();
        int moduleLocation = FileInterface::FindModuleLocation(inp, "LBGravity");
        const int RadiusLiquid = FileInterface::ReadParameterI(inp, moduleLocation, "iRadiusLiquid"); // Initial radius of solids
        const int RadiusSolid  = FileInterface::ReadParameterI(inp, moduleLocation, "iRadiusSolid");  // Initial radius of solids

        const int Nx = OPSettings.Grid.Nx;
        const int Ny = OPSettings.Grid.Ny;
        const int Nz = OPSettings.Grid.Nz;

        // Initialize Fluid phase (lattice Boltzmann representation)
        LB.InitializeSingle();

        // Initialize wall phase
        const size_t wall = Initializations::Single(Phase, 1, BC);
        Phase.FieldsProperties[wall].Mobile = false;

        // Initialize Fluid phase (Phase-Field representation)
        Initializations::Rectangular(Phase, 0, Nx-16, Ny-16, Nz, Nx/2, Ny/2, Nz/2, BC);

        // Initialize liquid sphere
        LB.InitializeSphere(RadiusLiquid, Nx/4, Ny/2, Nz/2, &LB.LiquidDensity[0]);

        // Initialize solid sphere
        Initializations::Sphere(Phase, 1, RadiusSolid, 2*Nx/4, Ny/2, Nz/2, BC);
        // Wet solid sphere
        LB.InitializeSphere(RadiusSolid, 2*Nx/4, Ny/2, Nz/2, &LB.Wetting[0][1]);
        LB.FinalizeInitialiation(Phase, Vel, BC);

        // Initialize solid sphere
        Initializations::Sphere(Phase, 2, RadiusSolid, 3*Nx/4, Ny/2, Nz/2, BC);
        // Wet solid sphere
        LB.InitializeSphere(RadiusSolid, 3*Nx/4, Ny/2, Nz/2, &LB.Wetting[0][2]);
        LB.FinalizeInitialiation(Phase, Vel, BC);

        ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);
    }

    // Start time loop
    std::vector<double> Forces;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteVTK())
        {
            Phase.WriteVTK(OPSettings, RTC.tStep);
            LB   .WriteVTK(OPSettings, Phase, RTC.tStep);
            LB .lbWriteVTK(OPSettings, Phase, RTC.tStep);
        }

        if (RTC.WriteRawData())
        {
            Phase.Write(OPSettings, RTC.tStep);
            LB   .Write(OPSettings, RTC.tStep);
        }

        if(RTC.WriteToScreen())
        {
            // Print out statistics to screen
            ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps);
            double VolumeSolids = 0;
            int    grainNr      = 0;
            dVector3 MomentumSolids = {0.0,0.0,0.0};
            for (size_t idx = 0; idx < Phase.FieldsProperties.size();idx++)
            if (Phase.FieldsProperties[idx].State == AggregateStates::Solid and
                Phase.FieldsProperties[idx].Mobile)
            {
                VolumeSolids += Phase.FieldsProperties[idx].Volume;

                const double   dx       = OPSettings.Grid.dx;
                const double   Mass     = Phase.FieldsProperties[idx].Volume*
                                          Phase.FieldsProperties[idx].Density*dx*dx*dx;
                const dVector3 Momentum = Phase.FieldsProperties[idx].Vcm*Mass;
                MomentumSolids += Momentum;

                ConsoleOutput::WriteSimple("Grain " + std::to_string(grainNr));
                ConsoleOutput::WriteLine("-");
                ConsoleOutput::Write("Position",     Phase.FieldsProperties[idx].Rcm);
                ConsoleOutput::Write("Acceleration", Phase.FieldsProperties[idx].Acm);
                ConsoleOutput::Write("Velocity",     Phase.FieldsProperties[idx].Vcm);
                ConsoleOutput::Write("Force",        Phase.FieldsProperties[idx].Force);
                ConsoleOutput::Write("Torque",       Phase.FieldsProperties[idx].Torque);
                ConsoleOutput::WriteBlankLine();
            }

            const dVector3 FluidMomentum = LB.CalculateFluidMomentum(Phase)[0];
            const double   FluidMass     = LB.CalculateFluidMass()[0];

            ConsoleOutput::Write("Solid volume ",  VolumeSolids);
            ConsoleOutput::Write("Solid momentum", MomentumSolids);
            ConsoleOutput::WriteBlankLine();
            ConsoleOutput::Write("Fluid mass",     FluidMass);
            ConsoleOutput::Write("Fluid momentum", FluidMomentum);
            ConsoleOutput::WriteBlankLine();
        }

        Phase.ClearGrainsForcesAndAccelerations();
        LB.Solve(Phase, Vel, BC);
        ISS.Calculate(Phase, BC, [] (int i ,int j, int k){return 0.0;}, RTC.dt);
        ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);
        Phase.Advect(ADHR, Vel, BC, LB, RTC.dt, RTC.tStep);
    }

#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif

    return 0;
}

void LocalLBM::InitializeSingle()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle       (i,j,k) = 0;
        DensityWetting (i,j,k,{0}) = VaporDensity[0];
        MomentumDensity(i,j,k,{0}) = {0.,0.,0.};
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LocalLBM::InitializeSphere(const double Radius, const int i0, const int j0,
        const int k0, const double* rho)
{
    const double loc_rho = (rho == nullptr) ? LiquidDensity[0] : *rho;
    const double D = 3;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        const int ii = i - i0;
        const int jj = j - j0;
        const int kk = k - k0;
        const double rr = sqrt(ii*ii + jj*jj + kk*kk);

        DensityWetting(i,j,k,{0}) = std::max((0.5*(loc_rho+VaporDensity[0]/dRho) -0.5*(loc_rho-VaporDensity[0]/dRho)*tanh(2.*(rr-Radius)/(D)))*dRho, DensityWetting(i,j,k,{0}));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
