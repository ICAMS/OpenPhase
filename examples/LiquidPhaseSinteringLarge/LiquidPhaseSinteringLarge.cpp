/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2018
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
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "Macros.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"
#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <random>
#include <sstream>

using namespace openphase;

class LocalLBM : public FlowSolverLBM                                           ///< Modification of FLowSolverLBM which are only used in this example
{
 public:
    LocalLBM(){};
    LocalLBM(Settings& locSettings, double in_dt):
        FlowSolverLBM(locSettings, in_dt){};

    void CalculateLiquidVolume(double& VolumeLiquid, double& VolumeVapor);

    void InitializeSingle();                                                    ///<  Initializes a single density
    void InitializeFinalize(PhaseField& Phase, Velocities& Vel,
            BoundaryConditions& BC);
};

int main(int argc, char *argv[])
{

#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_SINGLE, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
#endif

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);   // the result was too large
    //feenableexcept(FE_UNDERFLOW);  // the result was too small

    std::string InputFileName = "ProjectInput.opi";

    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("Sinter Liquid State Large");
    ConsoleOutput::WriteStandard("Source", InputFileName);

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

    int moduleLocation = FileInterface::FindModuleLocation(inp, "SinterLiquidStateLarge");
    const int    tStartAdvection = FileInterface::ReadParameterI(inp, moduleLocation, std::string("iTStartAdvection"));        // time to start advection
    //const int    tStartMelting   = FileInterface::ReadParameterI(inp, moduleLocation, std::string("iTStartMelting"));

    const size_t CoFluidPhaseIdx = 0;
    const size_t CoPhaseIdx = 1;
    const size_t WCPhaseIdx = 2;

    Settings OPSettings (DefaultInputFileName);

    AdvectionHR           ADHR  (OPSettings);
    BoundaryConditions    BC    (OPSettings);
    DoubleObstacle        DO    (OPSettings);
    InteractionSolidFluid ISF   (OPSettings);
    InteractionSolidSolid ISS   (OPSettings);
    InterfaceProperties   IP    (OPSettings);
    PhaseField            Phase (OPSettings);
    RunTimeControl        RTC   (OPSettings);
    Velocities            Vel   (OPSettings);

    LocalLBM LB (OPSettings, RTC.dt);

    // Set initial conditions
    if(RTC.Restart)
    {
        LB   .Read(OPSettings, BC, RTC.tStart);
        Phase.Read(OPSettings, BC, RTC.tStart);
        Vel  .Read(OPSettings, BC, RTC.tStart);
    }
    else
    {
        // Set start time step
        RTC.tStart = 0;

        const int Nx = OPSettings.Grid.TotalNx;
        const int Ny = OPSettings.Grid.Ny;
        const int Nz = OPSettings.Grid.Nz;
        const double AvRadius        = FileInterface::ReadParameterD(inp, moduleLocation, std::string("dAvRadius"));
        const double StdRadius       = FileInterface::ReadParameterD(inp, moduleLocation, std::string("StdRadius"));
        const double MinDistance     = FileInterface::ReadParameterD(inp, moduleLocation, std::string("dMinDistance"));

        std::string SRandomSeed = FileInterface::ReadParameterS(inp, moduleLocation, "RandomSeed");
        size_t RandomSeed;
        if (SRandomSeed == "random")
        {
            std::random_device rd;
            RandomSeed = rd();
            ConsoleOutput::Write("RandomSeed",RandomSeed);
        }
        else RandomSeed = std::stol(SRandomSeed);

        // Initialize vapor phase
        LB.InitializeSingle();
        Initializations::Single(Phase, CoFluidPhaseIdx, BC);

        int x0 = (Nx > 10) ? 10    : 0;
        int x1 = (Nx > 10) ? Nx-10 : Nx;
        int y0 = (Ny > 10) ? 10    : 0;
        int y1 = (Ny > 10) ? Ny-10 : Ny;
        int z0 = (Nz > 10) ? 10    : 0;
        int z1 = (Nz > 10) ? Nz-10 : Nz;

        std::mt19937_64 generator(RandomSeed);
        std::uniform_real_distribution<double> phase_gen(0,1);
        auto PhaseIndex = [&generator, &phase_gen](int i, int j, int k){return std::round(phase_gen(generator)+1);};
        std::vector<size_t> SpheresIdx = Initializations::FillRectangularWithSpheres(
                Phase, BC, OPSettings, PhaseIndex, AvRadius, StdRadius,
                x0, x1, y0, y1, z0, z1, MinDistance, RandomSeed);

        // Finalise initialization
        LB.InitializeFinalize(Phase, Vel, BC);
        LB.SetBoundaryConditions(BC);
    }
    std::cout << "Finished Initialization\n";

    // Calculate Mass of Co
    double MassCoRef = 0.0;
    const double  dx = OPSettings.Grid.dx;
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if (Phase.FieldsProperties[idx].Phase == CoPhaseIdx)
    {
        const double GrainVolume = Phase.FieldsProperties[idx].Volume*dx*dx*dx;
        MassCoRef += GrainVolume*Phase.FieldsProperties[idx].Density;
    }
    MassCoRef += LB.CalculateFluidMass()[0];

    std::ofstream Log("TextData/Log_Observables.csv", std::ofstream::app);
    Log << "tStep,"
           "VolumeWc,"
           "VolumeCo,"
           "VolumeLiquid,"
           "VolumeVapor,"
           "VolumeSolids,"
           "MassWC,"
           "MassSolidCo,"
           "MassFluidCo,"
           "MassCo,"
           "VolumeFractionVapor,"
           "VolumeFractionLiquid,"
           "VolumeFractionLiquidSolid,"
           "KineticEnergySolids,"
           "VelocityLimitApplied,"
           "TotalVelocityLimitApplied" << std::endl;
    Log << std::scientific << std::setprecision(16) << std::endl;

    // Start time loop
    const double dt = RTC.dt;
    std::vector<double> Forces;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        //  Output raw data to file
        if (RTC.WriteRawData())
        {
            Phase.Write(OPSettings, RTC.tStep);
            LB   .Write(OPSettings, RTC.tStep);
            Vel  .Write(OPSettings, RTC.tStep);
        }

        //  Output to VTK file
        if (RTC.WriteVTK())
        {
            Phase.WriteVTK(OPSettings, RTC.tStep);
            Vel  .WriteVTK(OPSettings, RTC.tStep);
            LB   .WriteVTK(OPSettings, Phase, RTC.tStep);
        }

        //  Output to screen
        if(RTC.WriteToScreen())
        {
            Phase.FieldsProperties.WriteTable(RTC.tStep, OPSettings);

            ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps);
            double VolumeSolids = 0;

            // Print out statistics
            double KinEnergySolids = 0.0;
            dVector3 MomentumSolids = {0.0,0.0,0.0};
            for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
            {
                Grain& grain = Phase.FieldsProperties[idx];
                if (grain.State == AggregateStates::Solid and grain.Mobile)
                {
                    Phase.FieldsProperties.WriteTable(idx,OPSettings,RTC.tStep);

                    VolumeSolids += grain.Volume;
                    const double   Mass     = grain.Volume*grain.Density*dx*dx*dx;
                    const dVector3 Momentum = grain.Vcm*Mass;
                    MomentumSolids  += Momentum;
                    KinEnergySolids += grain.Vcm*grain.Vcm*Mass*0.5;
                }
            }

            //const dVector3 FluidMomentum   = LB.CalculateFluidMomentum(Phase)[0];
            const double   FluidMass = LB.CalculateFluidMass()[0];

            const size_t Volume = OPSettings.Grid.Nx*OPSettings.Grid.Ny*OPSettings.Grid.Nz;
            //const size_t VolumeObstacles = LB.CountObstacleNodes();
            double VolumeLiquid, VolumeVapor;
            LB.CalculateLiquidVolume(VolumeLiquid, VolumeVapor);
            const double VolumeFluid = VolumeLiquid + VolumeVapor;
            const double VolumeFractionLiquid      = double(VolumeLiquid)/double(Volume);
            const double VolumeFractionVapor       = double(VolumeVapor )/double(Volume);
            const double VolumeFractionLiquidSolid = (VolumeFluid > 0.0) ? double(VolumeLiquid)/double(VolumeFluid) : 0.0;

            double VolumeWC = 0.0;
            double VolumeCo = 0.0;
            double MassWC   = 0.0;
            double MassCo   = 0.0;
            for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
            if (Phase.FieldsProperties[idx].Phase == CoPhaseIdx)
            {
                const double GrainVolume = Phase.FieldsProperties[idx].Volume*dx*dx*dx;
                VolumeCo += GrainVolume;
                MassCo   += GrainVolume*Phase.FieldsProperties[idx].Density;
            }
            else if (Phase.FieldsProperties[idx].Phase == WCPhaseIdx)
            {
                const double GrainVolume = Phase.FieldsProperties[idx].Volume*dx*dx*dx;
                VolumeWC += GrainVolume;
                MassWC   += GrainVolume*Phase.FieldsProperties[idx].Density;
            }

            ConsoleOutput::Write("Volume of solid WC ",             VolumeWC);
            ConsoleOutput::Write("Volume of solid Co ",             VolumeCo);
            ConsoleOutput::Write("Volume of liquid Co",             VolumeLiquid);
            ConsoleOutput::Write("Volume of Vaporous Co",           VolumeVapor);
            ConsoleOutput::Write("Total volume of solids",          VolumeSolids);
            //ConsoleOutput::Write("Volume of obstacle nodes ",       VolumeObstacles);
            ConsoleOutput::WriteBlankLine();
            ConsoleOutput::Write("Mass of solid WC ",               MassWC);
            ConsoleOutput::Write("Mass of solid Co ",               MassCo);
            ConsoleOutput::Write("Total mass of fluid Co",          FluidMass);
            ConsoleOutput::Write("Total mass of Co",                MassCo+FluidMass);
            ConsoleOutput::WriteBlankLine();
            ConsoleOutput::Write("Solid momentum",                  MomentumSolids);
            ConsoleOutput::Write("Solid kinetic energy",            KinEnergySolids);
            ConsoleOutput::WriteBlankLine();
            ConsoleOutput::Write("Vapor volume fraction",           VolumeFractionVapor);
            ConsoleOutput::Write("Liquid volume fraction",          VolumeFractionLiquid);
            //ConsoleOutput::Write("Solid obstacle volume fraction",  VolumeFractionObstacles);
            ConsoleOutput::Write("Liquid to solid volume fraction", VolumeFractionLiquidSolid);
            ConsoleOutput::WriteBlankLine();

            Log  << RTC.tStep                 << ","
                 << VolumeWC                  << ","
                 << VolumeCo                  << ","
                 << VolumeLiquid              << ","
                 << VolumeVapor               << ","
                 << VolumeSolids              << ","
                 << MassWC                    << ","
                 << MassCo                    << ","
                 << FluidMass                 << ","
                 << MassCo+FluidMass          << ","
                 << VolumeFractionVapor       << ","
                 << VolumeFractionLiquid      << ","
                 << VolumeFractionLiquidSolid << ","
                 << KinEnergySolids           << "\n";
        }

        Phase.ClearGrainsForcesAndAccelerations();
        LB.Solve(Phase, Vel, BC);

        // Calculate Solid Motion
        if (RTC.tStep > tStartAdvection)
        {
            ISS.Calculate(Phase, BC, [](int,int,int){return 0.;}, RTC.dt);
            ISF.CalculateSolidVelocities(Phase, Vel, BC, dt);
            Phase.Advect(ADHR, Vel, BC, LB, RTC.dt, RTC.tStep);
        }

        //TODO Melting is broken need update see examples Wetting
        //if (RTC.tStep >= tStartMelting)
        //{
        //    LB.DetectObstacles(Phase);
        //    Phase.Finalize(BC);

        //    // Calculate Solid Melting
        //    IP.Set(Phase);
        //    DO.CalculatePhaseFieldIncrements(Phase, IP, dG);
        //    Phase.MergeIncrements(BC, dt, true);
        //    LB.SetBoundaryConditions(BC);
        //    LB.CalculateDensityAndMomentum();

        //    // Calculate Mass of Co
        //    double MassCoNew = 0.0;
        //    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
        //    if (Phase.FieldsProperties[idx].Phase ==  CoPhaseIdx)
        //    {
        //        const double GrainVolume = Phase.FieldsProperties[idx].Volume*dx*dx*dx;
        //        MassCoNew   += GrainVolume*Phase.FieldsProperties[idx].Density;
        //    }
        //    MassCoNew += LB.CalculateFluidMass()[0];

        //    // Ensure mass conservation
        //    const double DeltaMassCo = MassCoRef - MassCoNew;
        //    //std::cout << DeltaMassCo << "\n";
        //    size_t FluidNodes = LB.CountFluidNodes();

        //    const double lbDeltaRho = DeltaMassCo/(LB.dRho*FluidNodes*dx*dx*dx);
        //    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        //    {
        //        if (not LB.Obstacle(i,j,k))
        //        LB.lbPopulations(i,j,k)({0})(0,0,0) += lbDeltaRho;
        //    }
        //    OMP_PARALLEL_STORAGE_LOOP_END
        //    LB.CalculateDensityAndMomentum();
        //}

    }// End of time loop
    return 0;
}

void LocalLBM::CalculateLiquidVolume(double& VolumeLiquid, double& VolumeVapor)
{
    VolumeLiquid   = 0;
    VolumeVapor    = 0;

    double New_lbLiquidDensity = 0;
    double New_lbVaporDensity = 1;

    STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0)
    if (not Obstacle(i,j,k))
    {
        const double DistLiquid = std::abs(DensityWetting(i,j,k,{0})/dRho - LiquidDensity[0]/dRho);
        const double DistVapor  = std::abs(DensityWetting(i,j,k,{0})/dRho - VaporDensity[0]/dRho);

        if (DistLiquid < DistVapor)
        {
            VolumeLiquid += 1.0;
            if (DensityWetting(i,j,k,{0}) > New_lbLiquidDensity)
                New_lbLiquidDensity = DensityWetting(i,j,k,{0});
        }
        else
        {
            VolumeVapor += 1.0;
            if (DensityWetting(i,j,k,{0}) < New_lbVaporDensity)
                New_lbVaporDensity = DensityWetting(i,j,k,{0});
        }
    }
    STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    double tmpLiquid = VolumeLiquid;
    double tmpVapor  = VolumeVapor;
    OP_MPI_Allreduce(&tmpLiquid, &(VolumeLiquid), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(&tmpVapor,  &(VolumeVapor),  1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif
    VolumeLiquid *= Grid.CellVolume(true);
    VolumeVapor  *= Grid.CellVolume(true);
}

void LocalLBM::InitializeSingle()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle(i,j,k) = 0;
        DensityWetting (i,j,k,{0}) = VaporDensity[0];
        MomentumDensity(i,j,k,{0}) = {0.,0.,0.};
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LocalLBM::InitializeFinalize(PhaseField& Phase, Velocities& Vel,
        BoundaryConditions& BC)
{
    // Initialize lattice Boltzmann populations for liquid phase
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        lbPopulations(i,j,k,{0}) = FlowSolverLBM::EquilibriumDistribution(DensityWetting(i,j,k,{0})/dRho, lbWeights, MomentumDensity(i,j,k,{0})/dm);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    DetectObstacles(Phase);
    SetObstacleNodes(Phase,Vel);
    FixPopulations();

    // Enforce mas and momentum conservation
    EnforceMassConservation();
}
