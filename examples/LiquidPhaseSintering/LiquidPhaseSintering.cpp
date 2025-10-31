/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2016
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
#include "AdvectionHR.h"
#include "DoubleObstacle.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Velocities.h"
#include <iomanip>

bool SortAbs (double i, double j) {return std::abs(i) < std::abs(j);}

using namespace openphase;

class LocalLBM : public FlowSolverLBM                                           ///< Modification of FLowSolverLBM which are only used in this example
{
 public:
    LocalLBM(Settings& locSettings, double in_dt): FlowSolverLBM(locSettings, in_dt){};

    double rho_l = 1.8141;                                                      ///<  Initial liquid density
    double rho_v = 0.2;                                                         ///<  Vapor gas density
    double FluidMass;

    double CalculateLiquidVolume(void);
    void EnforceLiquidConcentration(const double& CCurrent,
            const double& CGoal, const int& InitialRadius);

    void InitializeSingle();                                                    ///<  Initializes a single density
    void InitializeSphere(const double Radius,
            const int i0, const int j0, const int k0);                          ///<  Initializes a single sphere
};

void WriteDistanceStatisticToFile(std::ofstream& file, const int tStep,
        PhaseField& Phase, const unsigned int N,
        double& OldDistanceAverage,
        const char sep = ',');

void DetermineSolidParticlePositions(
        std::vector<double>& i0, std::vector<double>& j0, std::vector<double>& k0,
        const PhaseField& Phase, const int NSolids, const double RadiusSolid,
        const double RadiusLiquid);

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

    std::string InputFileName = DefaultInputFileName;
    if (argc > 1) InputFileName = argv[1];

    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("Liquid Phase Sintering");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::fstream inp(InputFileName, std::ios::in);
    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << InputFileName << "\" could not be opened";
        throw std::runtime_error(message.str());
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();


    int moduleLocation = FileInterface::FindModuleLocation(inp_data, "SinterLiquidState");
    //const double SolidDensity        = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("dSolidDensity"));//4.00;   //Lattice units
    const int    NSolids             = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("iNSolids"));//3;     // Number of solid particles
    const int    tInit               = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("iTInit"));//2000;  // time to start advection
    const int    tAdvect             = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("iTAdvect"));//2000;  // time to start advection
    const double LiquidConcentration = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("dLiquidConcentration")); //0.20;  // Concentration of liquid
    const int    InitialRadius       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("iInitialRadius"));//10    // Initial radius of solids

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
        Phase.Read(OPSettings, BC, RTC.tStart);
        LB   .Read(OPSettings, BC, RTC.tStart);
    }
    else
    {
        RTC.tStart = 0;

        const int iDistance = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("iDistance"));//10;    // Initial distance between particles

        // Radius of liquid and solid spheres
        const double RadiusSolid  = InitialRadius;
        const double RadiusLiquid = (1.0+3.0*LiquidConcentration)*RadiusSolid;

        // Radius of liquid and solid spheres
        ConsoleOutput::WriteLine("-");
        ConsoleOutput::WriteStandard("Liquid radius", RadiusLiquid);
        ConsoleOutput::WriteStandard("Solid radius", RadiusSolid);

        // Set position of particles
        std::vector<double> i0;
        std::vector<double> j0;
        std::vector<double> k0;

        i0.resize(NSolids,0);
        j0.resize(NSolids,0);
        k0.resize(NSolids,0);

        // Initialize vapor phase
        LB.InitializeSingle();
        Initializations::Single(Phase, 0, BC);

        // If the distance is to small, a single liquid bubble is initialed
        DetermineSolidParticlePositions(i0, j0, k0, Phase, NSolids,
                RadiusSolid, RadiusSolid + iDistance/2.0);

        // Initialize liquid container
        LB.InitializeSphere(RadiusLiquid,
                OPSettings.Grid.Nx/2.0, OPSettings.Grid.Ny/2.0, OPSettings.Grid.Nz/2.0);

        for(int n = 0; n < NSolids; n++)
        {
            // Initialize solid core
            Initializations::Sphere(Phase, 1, RadiusSolid, i0[n], j0[n], k0[n], BC);
        }

        LB.FinalizeInitialiation(Phase, Vel, BC);
        ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);
    }

    // Declare and initialize statistic files
    std::ofstream LogDistance(DefaultTextDir + "stat_distance.csv", std::ofstream::out|std::ofstream::app);
    std::ofstream Log(DefaultTextDir + "stat.csv", std::ofstream::out|std::ofstream::app);
    std::ofstream* LogSolid;
    LogSolid = new std::ofstream [NSolids];
    for(int n = 0; n < NSolids; n++)
    {
        LogSolid[n].open(DefaultTextDir + "stat_solid" + std::to_string(n) + ".csv", std::ofstream::out|std::ofstream::app);
        LogSolid[n] << std::scientific << std::endl;
        LogSolid[n] << std::setprecision(16) << std::endl;
    }

    // Write header
    Log <<"tStep,LiquidConcentration,LiquidVolume,SolidVolume,FluidMass,"
        <<"SolidInterfceEnergy,ObstacleNodes,AvForce,StdDevForce" << std::endl;
    Log << std::scientific << std::endl;
    Log << std::setprecision(16) << std::endl;

    // Start time loop
    double OldDistance = 0;
    std::vector<double> Forces;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteVTK())
        {
            Phase.WriteVTK(OPSettings, RTC.tStep);
            Vel  .WriteVTK(OPSettings, RTC.tStep);
            LB   .WriteVTK(OPSettings, Phase, RTC.tStep);
            LB .lbWriteVTK(OPSettings, Phase, RTC.tStep);
        }

        if (RTC.WriteRawData())
        {
            Phase.Write(OPSettings, RTC.tStep);
            LB   .Write(OPSettings, RTC.tStep);
        }

        Forces.push_back(Phase.FieldsProperties[1].Force.abs());

        if(RTC.WriteToScreen())
        {
            ConsoleOutput::WriteTimeStep(RTC);

            double VolumeSolids = 0;
            int    grainNr      = 0;

            // Print out statistics
            dVector3 MomentumSolids = {0.0,0.0,0.0};
            for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
            {
                if (Phase.FieldsProperties[idx].State == AggregateStates::Solid)
                {
                    VolumeSolids += Phase.FieldsProperties[idx].Volume;

                    const double   dx       = OPSettings.Grid.dx;
                    const double   Mass     = Phase.FieldsProperties[idx].Volume*
                                              Phase.FieldsProperties[idx].Density*dx*dx*dx;
                    const dVector3 Momentum = Phase.FieldsProperties[idx].Vcm*Mass;
                    MomentumSolids += Momentum;

                    ConsoleOutput::WriteSimple("Grain " + std::to_string(grainNr));
                    ConsoleOutput::WriteLine("-");
                    ConsoleOutput::Write("Density",  Phase.FieldsProperties[idx].Density);
                    ConsoleOutput::Write("Position", Phase.FieldsProperties[idx].Rcm);
                    ConsoleOutput::Write("Velocity", Phase.FieldsProperties[idx].Vcm);
                    ConsoleOutput::Write("Force",    Phase.FieldsProperties[idx].Force);
                    ConsoleOutput::Write("Torque",   Phase.FieldsProperties[idx].Torque);
                    ConsoleOutput::WriteBlankLine();

                    Phase.FieldsProperties[idx].WriteTable(LogSolid[grainNr], RTC.tStep);
                    grainNr++;
                }
            }

            const double SolidIntEnergy = DO.Energy(Phase,IP);
            const double VolumeLiquid   = LB.CalculateLiquidVolume();
            const double TempCLiquid    = VolumeLiquid/(VolumeLiquid + VolumeSolids);
            const double FluidMass      = LB.CalculateFluidMass()[0];

            const dVector3 FluidMomentum = LB.CalculateFluidMomentum(Phase)[0];

            const size_t ObstacleNodes = LB.CountObstacleNodes();

            // Sort forces in ascending order to achieve better numerical accuracy
            std::sort(Forces.begin(),Forces.end(), SortAbs);
            const size_t N = Forces.size() > 0 ? Forces.size() : 1;

            // Calculate average Force on wall
            long double AverageForce = 0.0;
            for (auto Force: Forces) AverageForce += Force;
            AverageForce /= N;

            // Calculate standard deviation of Force on wall
            long double stdDevForce = 0.0;
            for (auto Force: Forces)
            {
                stdDevForce += std::pow(Force - AverageForce,2);
            }
            stdDevForce /= N;
            stdDevForce = std::sqrt(stdDevForce);
            Forces.clear();

            ConsoleOutput::Write("Solid volume ",          VolumeSolids);
            ConsoleOutput::Write("Solid obstacle nodes ",  ObstacleNodes);
            ConsoleOutput::Write("Solid momentum",         MomentumSolids);
            ConsoleOutput::Write("Solid interface energy", SolidIntEnergy);
            ConsoleOutput::WriteBlankLine();
            ConsoleOutput::Write("Liquid volume",          VolumeLiquid);
            ConsoleOutput::Write("Liquid density",         LB.rho_l * LB.dRho);
            ConsoleOutput::Write("Liquid concentration",   TempCLiquid);
            ConsoleOutput::WriteBlankLine();
            ConsoleOutput::Write("Fluid mass",             FluidMass);
            ConsoleOutput::Write("Fluid momentum",         FluidMomentum);
            ConsoleOutput::WriteBlankLine();
            ConsoleOutput::Write("Average force",          AverageForce);
            ConsoleOutput::Write("stdDev of force",        stdDevForce);
            ConsoleOutput::WriteBlankLine();

            Log << RTC.tStep      << ","
                << TempCLiquid    << ","
                << VolumeLiquid   << ","
                << VolumeSolids   << ","
                << FluidMass      << ","
                << SolidIntEnergy << ","
                << ObstacleNodes  << ","
                << AverageForce   << ","
                << stdDevForce    << ","
                << std::endl;

            WriteDistanceStatisticToFile(LogDistance, RTC.tStep, Phase, NSolids, OldDistance);
        }

        // Enforce Liquid concentration (NOTE Example specific function)
        if (RTC.tStep < tInit)
        {
            double VolumeSolids = 0;
            for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
            {
                if (Phase.FieldsProperties[idx].State == AggregateStates::Solid)
                {
                    VolumeSolids += Phase.FieldsProperties[idx].Volume;
                }
            }

            const double VolumeLiquid = LB.CalculateLiquidVolume();
            const double TempCLiquid  = VolumeLiquid/(VolumeLiquid + VolumeSolids);

            LB.EnforceLiquidConcentration(TempCLiquid, LiquidConcentration, InitialRadius);
        }

        IP.Set(Phase,BC);
        Phase.ClearGrainsForcesAndAccelerations();
        LB.Solve(Phase, Vel, BC);
        if (RTC.tStep > tAdvect)
        {
            ISS.Calculate(Phase, BC, [](int,int,int){return 0.;}, RTC.dt);
            ISF.CalculateSolidVelocities(Phase, Vel, BC, RTC.dt);
            Phase.Advect(ADHR, Vel, BC, LB, RTC.dt, RTC.tStep);
        }
    }// End of time loop

#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif

    return 0;
}

void DetermineSolidParticlePositions(
        std::vector<double>& i0, std::vector<double>& j0, std::vector<double>& k0,
        const PhaseField& Phase, const int NSolids, const double RadiusSolid,
        const double RadiusLiquid)
{
    const double offset = std::max(RadiusLiquid, RadiusSolid+Phase.Grid.iWidth/2.0);
    if (NSolids == 1)
    {
        i0[0] = Phase.Grid.Nx/2.0;
        j0[0] = Phase.Grid.Ny/2.0;
        k0[0] = Phase.Grid.Nz/2.0;
    }
    else if (NSolids == 2)
    {
        i0[0] = Phase.Grid.Nx/2.0 - offset;
        j0[0] = Phase.Grid.Ny/2.0;
        k0[0] = Phase.Grid.Nz/2.0;

        i0[1] = Phase.Grid.Nx/2.0 + offset;
        j0[1] = Phase.Grid.Ny/2.0;
        k0[1] = Phase.Grid.Nz/2.0;
    }
    else if (NSolids == 3)
    {
        const double h  = std::sqrt(3.0) * offset;
        const double O  = tan(Pi/180.0*30) * offset;       // distance from base line to center of equilateral triangle

        i0[0] = Phase.Grid.Nx/2.0 - offset;
        j0[0] = Phase.Grid.Ny/2.0 - O;
        k0[0] = Phase.Grid.Nz/2.0;

        i0[1] = Phase.Grid.Nx/2.0 + offset;
        j0[1] = Phase.Grid.Ny/2.0 - O;
        k0[1] = Phase.Grid.Nz/2.0;

        i0[2] = Phase.Grid.Nx/2.0;
        j0[2] = Phase.Grid.Ny/2.0 - O + h;
        k0[2] = Phase.Grid.Nz/2.0;
    }
    else if (NSolids == 4)
    {
        const double O1 = tan(Pi/180.0*30) * offset;       // distance from base line to center of equilateral triangle
        const double h1 = std::sqrt(3.0)/2.0 * 2.0*offset; // height of equilateral triangle
        const double h2 = std::sqrt(6.0)/3.0 * 2.0*offset; // median of tetrahedron

        i0[0] = Phase.Grid.Nx/2.0 - offset;
        j0[0] = Phase.Grid.Ny/2.0 - O1;
        k0[0] = Phase.Grid.Nz/2.0 - h2/3.0;

        i0[1] = Phase.Grid.Nx/2.0 + offset;
        j0[1] = Phase.Grid.Ny/2.0 - O1;
        k0[1] = Phase.Grid.Nz/2.0 - h2/3.0;

        i0[2] = Phase.Grid.Nx/2.0;
        j0[2] = Phase.Grid.Ny/2.0 - O1 + h1;
        k0[2] = Phase.Grid.Nz/2.0 - h2/3.0;

        i0[3] = Phase.Grid.Nx/2.0;
        j0[3] = Phase.Grid.Ny/2.0;
        k0[3] = Phase.Grid.Nz/2.0 + h2*2.0/3.0;
    }
    else
    {
        //TODO Random distribution
    }
}

void WriteDistanceStatisticToFile(std::ofstream& file, const int tStep,
        PhaseField& Phase, const unsigned int N, double& OldDistanceAverage,
        const char sep)
{
    file << std::scientific << std::endl;
    file << std::setprecision(16) << std::endl;
    if (tStep == 0)
    {
        file << "tStep" << sep;
        for(size_t n = 0; n < Phase.FieldsProperties.size(); ++n)
            for(size_t m = n + 1; m < Phase.FieldsProperties.size(); ++m)
            {
                Grain grainA = Phase.FieldsProperties[n];
                Grain grainB = Phase.FieldsProperties[m];

                if (grainA.State == AggregateStates::Solid and grainB.State == AggregateStates::Solid)
                {
                    file << "R" << n << "_" << m << sep;
                }
            }
        file << "RAverage" << sep << "dRAverage" << std::endl;
    }
    else
    {
        // Calculate distances
        double DistanceAverage = 0.0;
        file  << tStep << sep;
        for(size_t n = 0; n < Phase.FieldsProperties.size(); ++n)
            for(size_t m = n + 1; m < Phase.FieldsProperties.size(); ++m)
            {
                Grain grainA = Phase.FieldsProperties[n];
                Grain grainB = Phase.FieldsProperties[m];

                if (grainA.State == AggregateStates::Solid and grainB.State == AggregateStates::Solid)
                {
                    const double buffer = (grainA.Rcm - grainB.Rcm).abs();
                    DistanceAverage += buffer;
                    file << buffer << sep;
                }
            }
        DistanceAverage /= (N*(N+1.0)/2.0-N);
        file << DistanceAverage << sep;

        if (tStep > 1) file << 0 << std::endl;
        else file << DistanceAverage - OldDistanceAverage << std::endl;

        OldDistanceAverage = DistanceAverage;
    }
}

void LocalLBM::EnforceLiquidConcentration(const double& CCurrent, const double& CGoal, const int& InitialRadius)
{
    if (std::abs(CCurrent-CGoal) > 0.0001)
    {
        const double DeltaDensity = -0.1/InitialRadius*(CGoal-CCurrent)*rho_l;

        // Fix density
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (not Obstacle(i,j,k))
        {
            const double DistLiquid = abs(DensityWetting(i,j,k,{0}) - rho_l);
            const double DistVapor  = abs(DensityWetting(i,j,k,{0}) - rho_v);

            if (DistLiquid < DistVapor)
            {
                lbPopulations(i,j,k,{0}) = lbPopulations(i,j,k,{0}) -
                    FlowSolverLBM::EquilibriumDistribution(DensityWetting(i,j,k,{0})/dRho + DeltaDensity, lbWeights,{0.,0.,0.}) +
                    FlowSolverLBM::EquilibriumDistribution(DensityWetting(i,j,k,{0})/dRho, lbWeights, {0.,0.,0.});
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // Update entire fluid mass in order to ensure mass conservation
        FluidMass = CalculateFluidMass()[0];
    }
}

double LocalLBM::CalculateLiquidVolume(void)
{
    int VolumeLiquid   = 0;
    //int VolumeVapor    = 0;
    //int VolumeObstacle = 0;

    //double MassLiquid = 0;
    //double MassVapor  = 0;

    double New_rho_l = 0;
    double New_rho_v = 1;

    STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0)
    if (not Obstacle(i,j,k))
    {
        const double DistLiquid = abs(DensityWetting(i,j,k,{0})/dRho - rho_l);
        const double DistVapor  = abs(DensityWetting(i,j,k,{0})/dRho - rho_v);

        if (DistLiquid < DistVapor)
        {
            ++VolumeLiquid;
            if (DensityWetting(i,j,k,{0}) > New_rho_l)
                New_rho_l = DensityWetting(i,j,k,{0})/dRho;
            //MassLiquid += DensityWetting(i,j,k,{0})/dRho;
        }
        else
        {
            //++VolumeVapor;
            if (DensityWetting(i,j,k,{0}) < New_rho_v)
                New_rho_v = DensityWetting(i,j,k,{0})/dRho;
            //MassVapor += DensityWetting(i,j,k,{0})/dRho;
        }
    }
    //else VolumeObstacle++;
    STORAGE_LOOP_END

    //Update liquid and vapor density
    //rho_l = MassLiquid/VolumeLiquid;
    //rho_v = MassVapor/VolumeVapor;

    rho_l = New_rho_l;
    rho_v = New_rho_v;


    return VolumeLiquid;
}

void LocalLBM::InitializeSingle()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle(i,j,k) = 0;
        DensityWetting(i,j,k,{0}) = rho_v*dRho;
        MomentumDensity(i,j,k,{0}) = {0.,0.,0.};
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LocalLBM::InitializeSphere(const double Radius,
        const int i0, const int j0, const int k0)
{
    const double D = 3;
    const int Bcells = DensityWetting.Bcells();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,Bcells,)
    {
        const int ii = i - i0;
        const int jj = j - j0;
        const int kk = k - k0;
        const int rr = sqrt(ii*ii + jj*jj + kk*kk);

        DensityWetting(i,j,k,{0}) = std::max((0.5*(rho_l+rho_v) -0.5*(rho_l-rho_v)*tanh(2.*(rr-Radius)/(D)))*dRho, DensityWetting(i,j,k,{0}));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
