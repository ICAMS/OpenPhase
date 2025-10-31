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

 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Dmitri Medvedev; Amol Subhedar;
 *                         Marvin Tegeler; Raphael Schiedung; Reza Namdar
 *
 */

#include "BoundaryConditions.h"
#include "Composition.h"
#include "FluidDynamics/BenziGas.h"
#include "FluidDynamics/D3Q27.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/VanDerWaalsGas.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"
#include "Velocities.h"

namespace openphase
{
using namespace std;

constexpr double lbcs2 = 1.0/3.0;                                               ///< Speed of sound squared [lattice units]

inline double psi(double rho, double rho0)                                      ///< Density potential function
{
    assert(rho > 0.0);
    return 1.0 - exp(-rho/rho0);
}

inline double Phi(
        const double lbDensity,
        const double ReducedPressure,
        const double GasParameter)
{
    const double potential = GasParameter*ReducedPressure - lbDensity/3.0;
    assert(potential <= 0.0);
    return std::sqrt(-potential);
}

 void CalculateDistancePeriodic(dVector3 A, dVector3 B,
                dVector3 &dist, double Nx, double Ny, double Nz)
{
    dist[0] = A[0]-B[0];
    if (std::abs(dist[0]) > std::abs(dist[0]+Nx)) dist[0] += Nx;
    else if (std::abs(dist[0]) > std::abs(dist[0]-Nx)) dist[0] -= Nx;
    dist[1] = A[1]-B[1];
    if (std::abs(dist[1]) > std::abs(dist[1]+Ny)) dist[1] += Ny;
    else if (std::abs(dist[1]) > std::abs(dist[1]-Ny)) dist[1] -= Ny;
    dist[2] = A[2]-B[2];
    if (std::abs(dist[2]) > std::abs(dist[2]+Nz)) dist[2] += Nz;
    else if (std::abs(dist[2]) > std::abs(dist[2]-Nz)) dist[2] -= Nz;
};

FlowSolverLBM::FlowSolverLBM(Settings& locSettings, double in_dt, const std::string InputFileName)
{
    
    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname,"ReadInput()");
        std::exit(EXIT_FAILURE);
    };
    stringstream inp;
    inp << inpF.rdbuf();

    Initialize(locSettings,inp,in_dt);
    //ReadInput(InputFileName); // NOTE will be calles in Initialize
}

void FlowSolverLBM::Initialize(Settings& locSettings, std::stringstream& inp, double in_dt)
{
    thisclassname = "FlowSolverLBM";

    Grid = locSettings.Grid;

    dt      = in_dt;
    Nphases = locSettings.Nphases;
    Ncomp   = (locSettings.Ncomp > 0) ? locSettings.Ncomp - 1 : 0;

    Pth    = 1.0e05;
    PthOld = 1.0e05;
    Pth0   = 1.0e05;
    Poutlet = 1e-6;

    cs2 = lbcs2*Grid.dx*Grid.dx/dt/dt;

    ObstaclesChanged = false;

    Bcells = Grid.Bcells;
    ElementNames = locSettings.ElementNames;
    ReadInput(inp); // NOTE N_Fluid_Comp needs to be read

    Obstacle.Allocate               (Grid, Bcells);
    ObstacleAppeared.Allocate       (Grid, Bcells);
    ObstacleChangedDensity.Allocate (Grid, Bcells);
    ObstacleVanished.Allocate       (Grid, Bcells);
    DensityWetting.Allocate         (Grid, {N_Fluid_Comp}, Bcells);
    ForceDensity.Allocate           (Grid, {N_Fluid_Comp}, Bcells);
    MomentumDensity.Allocate        (Grid, {N_Fluid_Comp}, Bcells);
    lbPopulations.Allocate          (Grid, {N_Fluid_Comp}, Bcells);
    lbPopulationsTMP.Allocate       (Grid, {N_Fluid_Comp}, Bcells);
    nut.Allocate                    (Grid, {N_Fluid_Comp}, Bcells);
    HydroPressure.Allocate          (Grid, {N_Fluid_Comp}, Bcells);
    DivVel.Allocate                 (Grid, {N_Fluid_Comp}, Bcells);
    GradRho.Allocate                (Grid, {N_Fluid_Comp}, Bcells);
    
    switch(Grid.Active())
    {
        case 1:
        {
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                lbWeights[ii][jj][kk] = LBStencil1D[ii][jj][kk];
            }
            break;
        }
        case 2:
        {
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                lbWeights[ii][jj][kk] = LBStencil2D[ii][jj][kk];
            }
            break;
        }
        case 3:
        {
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                lbWeights[ii][jj][kk] = LBStencil3D[ii][jj][kk];
            }
            break;
        }
        default:
        {

        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle               (i,j,k) = false;
        ObstacleAppeared       (i,j,k) = false;
        ObstacleChangedDensity (i,j,k) = false;
        ObstacleVanished       (i,j,k) = false;

        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            DensityWetting   (i,j,k,{n}) = 0.0;
            HydroPressure    (i,j,k,{n}) = 0.0;
            DivVel           (i,j,k,{n}) = 0.0;
            GradRho          (i,j,k,{n}).set_to_zero();
            ForceDensity     (i,j,k,{n}).set_to_zero();
            MomentumDensity  (i,j,k,{n}).set_to_zero();
            lbPopulations    (i,j,k,{n}).set_to_zero();
            lbPopulationsTMP (i,j,k,{n}).set_to_zero();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    InitialTotalMass = -1.0;

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    initialized = true;
    ConsoleOutput::Write(thisclassname, "Initialized");
}

void FlowSolverLBM::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("FlowSolverLBM input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname,"ReadInput()");
        std::exit(EXIT_FAILURE);
    };
    stringstream inp;
    inp << inpF.rdbuf();

    ReadInput(inp);

    inpF.close();
}

void FlowSolverLBM::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    N_Fluid_Comp = FileInterface::ReadParameterD(inp, moduleLocation, std::string("N_FLUID_COMP"), false, 1);
    nu.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream NUstr;
        NUstr << "NU[" << n << "]";
        nu[n] = FileInterface::ReadParameterD(inp, moduleLocation, NUstr.str());
    }

    Do_BounceBack            = FileInterface::ReadParameterB(inp, moduleLocation, std::string("BOUNCEBACK"));
    Do_Benzi                 = FileInterface::ReadParameterB(inp, moduleLocation, std::string("BENZI"), false, false);
    Do_BounceBackElastic     = FileInterface::ReadParameterB(inp, moduleLocation, std::string("BBELASTIC"), false, false);
    Do_Drag                  = FileInterface::ReadParameterB(inp, moduleLocation, std::string("DRAG"));
    Do_FixPopulations        = FileInterface::ReadParameterB(inp, moduleLocation, std::string("bFixPop"),  false, true);
    Do_Gravity               = FileInterface::ReadParameterB(inp, moduleLocation, std::string("GRAVITY"), false, false);
    Do_Kupershtokh           = FileInterface::ReadParameterB(inp, moduleLocation, std::string("KUPERSHTOKH"), false, false);
    Do_ThermalComp           = FileInterface::ReadParameterB(inp, moduleLocation, std::string("THERMALCOMP"), false, false);
    Do_DI          			 = FileInterface::ReadParameterB(inp, moduleLocation, std::string("Diffuse_Interface"), false, false);
    FluidRedistributionRange = FileInterface::ReadParameterI(inp, moduleLocation, std::string("FluidRedistributionRange"), false, 1);
    ParaKuper                = FileInterface::ReadParameterD(inp, moduleLocation, std::string("ParaKuper"), false, -0.0152);
    h_star                   = FileInterface::ReadParameterD(inp, moduleLocation, std::string("H_STAR"), Do_Drag, 0.0);
    Do_TwoPhase = (Do_Benzi or Do_Kupershtokh) ? true : false;
    Do_GuoForcing            = FileInterface::ReadParameterB(inp, moduleLocation, std::string("GUO_FORCING"), false,  Do_TwoPhase);
    Do_EDForcing             = FileInterface::ReadParameterB(inp, moduleLocation, std::string("ED_FORCING"),  false, !Do_TwoPhase);
    GradRho_Upwind    	     = FileInterface::ReadParameterB(inp, moduleLocation, std::string("GradRho_Upwind"), false, false); 
	GradRho_Central          = FileInterface::ReadParameterB(inp, moduleLocation, std::string("GradRho_Central"), false, false); 
	GradRho_VanLeer          = FileInterface::ReadParameterB(inp, moduleLocation, std::string("GradRho_VanLeer"), false, true);

    if (Bcells < FluidRedistributionRange)
    {
        std::stringstream message;
            message << Bcells << " boundary cells are too few, at least "
                    << FluidRedistributionRange << " boundary cells are needed!";

        ConsoleOutput::WriteExit(message.str(), thisclassname, "ReadInput");
        std::exit(EXIT_FAILURE);
    }

    // Determine density discretization
    // NOTE that dRho is the equilibrium fluid density [Kg/m^2] if no liquid-vapor phase separation is considered!
    dRho = FileInterface::ReadParameterD(inp, moduleLocation, std::string("dRho"));
    dM   = dRho*Grid.CellVolume(true);
    dP   = dM/Grid.dx/dt/dt;
    dm   = dM/Grid.dx/Grid.dx/dt;
    df   = dM/Grid.dx/Grid.dx/dt/dt;
    dnu  = Grid.dx*Grid.dx/dt;
    dv   = Grid.dx/dt;

    if (Ncomp > 0)
    {
        dDensity_dc.resize(Ncomp);
        for (size_t n = 0; n < Ncomp; ++n)
        {
            std::stringstream dDensity_dc_str;
            dDensity_dc_str << "DDENSITY_DC_" << ElementNames[n];
            dDensity_dc[n] = FileInterface::ReadParameterD(inp, moduleLocation, dDensity_dc_str.str(), Do_Gravity,0.0);
        }
    }

    // GA Gravitational acceleration
    for (size_t i = 0; i < 3; i++)
    {
        std::stringstream GAstr;
        GAstr << "G0_" << i;
        GA[i] = FileInterface::ReadParameterD(inp, moduleLocation, GAstr.str(), Do_Gravity, 0.0);
    }

    // Read Benzi interaction parameter
    Gb.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        Gb[n].resize(N_Fluid_Comp);
        for (size_t m = 0; m < N_Fluid_Comp; ++m)
        {
            std::stringstream Gbstr;
            Gbstr << "GB_" << n << "_" << m;
            Gb[n][m] = FileInterface::ReadParameterD(inp, moduleLocation, Gbstr.str(), Do_Benzi, (n==m) ? 1 : 0);

        }
    }

    // Read Kupershtokh interaction parameter
    lbGK.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbGK[n].resize(N_Fluid_Comp);
        for (size_t m = 0; m < N_Fluid_Comp; ++m)
        {
            std::stringstream GKstr;
            GKstr << "lbGK_" << n << "_" << m;
            lbGK[n][m] = FileInterface::ReadParameterD(inp, moduleLocation, GKstr.str(), Do_Kupershtokh, (n==m) ? 1 : 0);
        }
    }

    rho_0.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream RHO0str;
        RHO0str << "RHO0_" << n;
        rho_0[n] = FileInterface::ReadParameterD(inp, moduleLocation, RHO0str.str(), false, dRho);
    }

    CriticalDensity.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream RHOCstr;
        RHOCstr << "RHOC_" << n;
        CriticalDensity[n] = FileInterface::ReadParameterD(inp, moduleLocation, RHOCstr.str(), Do_Kupershtokh, dRho);
    }

    lbTemperature.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream TEMPstr;
        TEMPstr << "TEMP_" << n;
        lbTemperature[n] = FileInterface::ReadParameterD(inp, moduleLocation, TEMPstr.str(), Do_Kupershtokh, 0.0);
    }

    lbCriticalTemperature.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
       std::stringstream TEMPCstr;
       TEMPCstr << "TEMPC_" << n;
       lbCriticalTemperature[n] = FileInterface::ReadParameterD(inp, moduleLocation, TEMPCstr.str(), Do_Kupershtokh, 0.0);
    }

    GasParameter.resize(N_Fluid_Comp);
    lbCriticalPressure.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream TEMPCstr;
        TEMPCstr << "CRITP_" << n;
        lbCriticalPressure[n] = FileInterface::ReadParameterD(inp, moduleLocation, TEMPCstr.str(), Do_Kupershtokh, 0.0);
        GasParameter[n]       = lbCriticalPressure[n]/CriticalDensity[n]*dRho*(dt*dt)/(Grid.dx*Grid.dx);
    }

    LiquidDensity .resize(N_Fluid_Comp);
    VaporDensity  .resize(N_Fluid_Comp);
    SurfaceTension.resize(N_Fluid_Comp);
    InterfaceWidth.resize(N_Fluid_Comp);
    if (Do_Benzi)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {

        const BenziGas::EquilibriumValues_t Equilibrium =
            BenziGas::EquilibriumValues(Gb[n][n]/dRho);

        LiquidDensity  [n] = Equilibrium.lbLiquidDensity*dRho;
        VaporDensity   [n] = Equilibrium.lbVaporDensity*dRho;
        SurfaceTension [n] = Equilibrium.lbSurfaceTension*dM/dt/dt;
        InterfaceWidth [n] = 2.0; //TODO calculate

        ConsoleOutput::Write("Vapor density   ["+std::to_string(n)+"]", VaporDensity   [n]);
        ConsoleOutput::Write("Liquid density  ["+std::to_string(n)+"]", LiquidDensity  [n]);
        ConsoleOutput::Write("Surface tension ["+std::to_string(n)+"]", SurfaceTension [n]);
        ConsoleOutput::Write("Interface width ["+std::to_string(n)+"]", InterfaceWidth [n]);
    }
    else if (Do_Kupershtokh)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        const VanDerWaalsGas::EquilibriumValues_t Equilibrium =
            VanDerWaalsGas::EquilibriumValues(lbTemperature[n], lbCriticalTemperature[n]);

        LiquidDensity  [n] = Equilibrium.LiquidDensity*dRho;
        VaporDensity   [n] = Equilibrium.VaporDensity*dRho;
        SurfaceTension [n] = 0.0*dM/dt/dt; //TODO
        InterfaceWidth [n] = 0.402819*std::pow(GasParameter[n],-0.470509);

        ConsoleOutput::Write("Vapor density   ["+std::to_string(n)+"]", VaporDensity   [n]);
        ConsoleOutput::Write("Liquid density  ["+std::to_string(n)+"]", LiquidDensity  [n]);
        ConsoleOutput::Write("Surface tension ["+std::to_string(n)+"]", SurfaceTension [n]);
        ConsoleOutput::Write("Interface width ["+std::to_string(n)+"]", InterfaceWidth [n]);
    }
    else
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        LiquidDensity  [n] = dRho;
        VaporDensity   [n] = dRho;
        InterfaceWidth [n] = 1.0;
    }

    U0X  = FileInterface::ReadParameterD(inp, moduleLocation, std::string("U0X"), false, 0.0);
    U0Y  = FileInterface::ReadParameterD(inp, moduleLocation, std::string("U0Y"), false, 0.0);
    U0Z  = FileInterface::ReadParameterD(inp, moduleLocation, std::string("U0Z"), false, 0.0);
    Pth0 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Pth0"), false, 1.0e5);
    Poutlet = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Poutlet"), false, 1.0e5);
    Pth    = Pth0;
    PthOld = Pth0;
    Poutlet = Poutlet - 1.0e5 + 1.0e-4;

    Wetting.resize(N_Fluid_Comp);
    if (Do_TwoPhase)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        Wetting[n].resize(Nphases);
        for (size_t m = 0; m < Nphases; ++m)
        {
            std::stringstream TEMPCstr;
            TEMPCstr << "Wetting_" << n << "_" << m;
            Wetting[n][m] = FileInterface::ReadParameterD(inp, moduleLocation, TEMPCstr.str());
        }
    }

    FluidMass.resize(N_Fluid_Comp, 0.0);

    // NOTE: the minimum and maximum kinematic viscosity values are based
    // on experience and may be adjusted.
    const double nu_min = 0.05*Grid.dx*Grid.dx/dt;
    const double nu_max = 1.00*Grid.dx*Grid.dx/dt;
    const double nu_opt = 0.5/3.0*Grid.dx*Grid.dx/dt; // optimal value !

    // Set relaxation parameter tau
    lbtau.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        if ((nu[n] < nu_min) or (nu[n] > nu_max))
        {
            std::stringstream message;
            message << "Bad kinematic viscosity"
                    <<            "nu_lb = "<< std::scientific << nu[n]
                    << "\n\t  (Min nu_lb = "<< std::scientific << nu_min
                    <<      "; Opt nu_lb = "<< std::scientific << nu_opt
                    <<      "; Max nu_lb = "<< std::scientific << nu_max
                    <<")!\n\t  Simulation may be unstable! Adjust dx or dt.";
            ConsoleOutput::WriteWarning(message.str(), thisclassname, "Initialize");
        }

        // Set relaxation parameter tau
        lbtau [n] = 3*nu[n]/Grid.dx/Grid.dx*dt + 0.5;
        ConsoleOutput::WriteStandard("lbtau["+std::to_string(n)+"]", lbtau[n]);
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void FlowSolverLBM::Update_dt(double _dt)
{
    dt  = _dt;
    dM  = dRho*Grid.CellVolume(true);
    dP  = dM/Grid.dx/dt/dt;
    dm  = dM/Grid.dx/Grid.dx/dt;
    df  = dM/Grid.dx/Grid.dx/dt/dt;
    dnu = Grid.dx*Grid.dx/dt;
    // NOTE: the minimum and maximum kinematic viscosity values are based
    // on experience and may be adjusted.
    //const double nu_min = 0.05*dx*dx/dt;
    //const double nu_max = 1.00*dx*dx/dt;
    //const double nu_opt = 0.5/3.0*dx*dx/dt; // optimal value !

    // Set relaxation parameter tau
    lbtau.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbtau[n] = 3*nu[n]/Grid.dx/Grid.dx*dt + 0.5;
    }
}

void FlowSolverLBM::Remesh(int newNx, int newNy, int newNz,
                           const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);

    DensityWetting  .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    Obstacle        .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    MomentumDensity .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    ForceDensity    .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    lbPopulations   .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    lbPopulationsTMP.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
}

bool FlowSolverLBM::Read(const Settings& locSettings, const BoundaryConditions& BC, const int tStep)
{
#ifdef MPI_PARALLEL
    std::string FileName = FileInterface::MakeFileName(locSettings.InputRawDataDir, thisclassname+'_'+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    std::string FileName = FileInterface::MakeFileName(locSettings.InputRawDataDir, thisclassname+'_', tStep, ".dat");
#endif

    std::ifstream inp(FileName.c_str(), std::ios::in | std::ios::binary);

    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be opened";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    };

    int locNx = Grid.Nx;
    int locNy = Grid.Ny;
    int locNz = Grid.Nz;
    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));

    if(locNx != Grid.Nx or locNy != Grid.Ny or locNz != Grid.Nz)
    {
        std::stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: ("
                << locNx << ", "
                << locNy << ", "
                << locNz << ") grid points.\n"
                << "Required data dimensions: ("
                << Grid.Nx << ", "
                << Grid.Ny << ", "
                << Grid.Nz << ") grid points.\n";
#ifdef DEBUG
        throw std::runtime_error(message.str());
#else
        return false;
#endif
    }

    size_t locN_Fluid_Comp = N_Fluid_Comp;
    inp.read(reinterpret_cast<char*>(&N_Fluid_Comp), sizeof(size_t));

    if(locN_Fluid_Comp != N_Fluid_Comp)
    {
        std::stringstream message;
        message << "Inconsistent number of fluid components\n"
                << "Input number of fluid components: "
                << locN_Fluid_Comp << "\n"
                << "Required number of fluid components: "
                << N_Fluid_Comp << "\n";
#ifdef DEBUG
        throw std::runtime_error(message.str());
#else
        return false;
#endif
    }

    STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
        {
            inp.read(reinterpret_cast<char*>(&lbPopulations(i,j,k,{n})(ii,jj,kk)), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int m = 0; m < 3; ++m)
        {
            inp.read(reinterpret_cast<char*>(&ForceDensity(i,j,k,{n})[m]), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Obstacle,0)
    {
        inp.read(reinterpret_cast<char*>(&Obstacle(i,j,k)), sizeof(bool));
    }
    STORAGE_LOOP_END

    CalculateDensityAndMomentum();
    SetBoundaryConditions(BC);

    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

void FlowSolverLBM::SetUniformVelocity(const BoundaryConditions& BC,
        const dVector3 U0)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        DensityWetting  (i,j,k,{n}) = dRho;
        MomentumDensity (i,j,k,{n}) = U0*DensityWetting(i,j,k,{n});
        lbPopulations   (i,j,k,{n}) = EquilibriumDistribution(1.0, lbWeights, MomentumDensity(i,j,k,{n})/dm);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

bool FlowSolverLBM::Write(const Settings& locSettings, const int tStep) const
{
#ifdef MPI_PARALLEL
    std::string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname+'_'+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    std::string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname+'_', tStep, ".dat");
#endif

    std::ofstream out(FileName.c_str(), std::ios::out | std::ios::binary);

    if (!out)
    {
        std::stringstream message;
        message << "File \"" << FileName
            << "\" could not be created! Terminating!!!" << std::endl;
        throw std::runtime_error(message.str());
    };

    out.write(reinterpret_cast<const char*>(&Grid.Nx),           sizeof(int));
    out.write(reinterpret_cast<const char*>(&Grid.Ny),           sizeof(int));
    out.write(reinterpret_cast<const char*>(&Grid.Nz),           sizeof(int));
    out.write(reinterpret_cast<const char*>(&N_Fluid_Comp), sizeof(size_t));

    STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
        {
            const double value = lbPopulations(i,j,k,{n})(ii,jj,kk);
            out.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int m = 0; m < 3; ++m)
        {
            const double value = ForceDensity(i,j,k,{n})[m];
            out.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Obstacle,0)
    {
        const bool value = Obstacle(i,j,k);
        out.write(reinterpret_cast<const char*>(&value), sizeof(bool));
    }
    STORAGE_LOOP_END

    out.close();
    return true;
}

void FlowSolverLBM::WriteVTK(const Settings& locSettings, const PhaseField& Phase, const int tStep, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ListOfFields.push_back((VTK::Field_t) {"DensityWetting_" +std::to_string(n), [n, this] (int i,int j,int k){return DensityWetting(i,j,k,{n});}});
    }
    ListOfFields.push_back((VTK::Field_t) {"Density",  [this, &Phase](int i,int j,int k){return Density(Phase,i,j,k);}});
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return Pressure(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t) {"Obstacle", [this](int i,int j,int k){return int(Obstacle(i,j,k));}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+'_', tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void FlowSolverLBM::lbWriteVTK(const Settings& locSettings, const PhaseField& Phase, const int tStep, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ListOfFields.push_back((VTK::Field_t) {"DensityWetting_" +std::to_string(n), [n, this] (int i,int j,int k){return DensityWetting(i,j,k,{n})/dRho;}});
    }
    ListOfFields.push_back((VTK::Field_t) {"Density",  [this, &Phase](int i,int j,int k){return Density(Phase,i,j,k)/dRho;}});
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return Pressure(i,j,k)/dP;}});
    ListOfFields.push_back((VTK::Field_t) {"Obstacle", [this](int i,int j,int k){return int(Obstacle(i,j,k));}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_lattice_units_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void FlowSolverLBM::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetXVector(lbPopulations);
    BC.SetYVector(lbPopulations);
    BC.SetZVector(lbPopulations);

    BC.SetX(DensityWetting);
    BC.SetY(DensityWetting);
    BC.SetZ(DensityWetting);

    BC.SetXVector(MomentumDensity);
    BC.SetYVector(MomentumDensity);
    BC.SetZVector(MomentumDensity);
}

void FlowSolverLBM::CalculateForceTwoPhase(const int i, const int j,
        const int k, PhaseField& Phase)
{
    for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
    for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
    for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
    {
        const dVector3 vel {ii*Grid.dx/dt,jj*Grid.dx/dt,kk*Grid.dx/dt};

        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        for (size_t m = 0; m < N_Fluid_Comp; ++m)
        {
            double tmp = 0.0;
            if (Do_Benzi)
            {
                // Benzi et. al. (2006) - Mesoscopic modeling of a two-phase flow
                // in the presence of boundaries
                tmp = - Gb[n][m]/dRho * psi(DensityWetting(i,j,k,{n}), rho_0[n]) * lbWeights[ii+1][jj+1][kk+1] * psi(DensityWetting(i+ii,j+jj,k+kk,{m}), rho_0[m]);
            }
            else if (Do_Kupershtokh)
            {
                // Calculate phase sepration force according to:
                // Kupershtokh, Medvedev and Karpov, 2009, On equations of state in a lattice Boltzmann method
                // Kupershtokh, Medvedev and Gribanov, 2017, Thermal lattice Boltzmann method for multiphase flows

                // Evaluate equation of state (Van der Waals)
                const double p_000 = VanDerWaalsGas::ReducedPressure(DensityWetting(i   ,j   ,k   ,{n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                const double p_ijk = VanDerWaalsGas::ReducedPressure(DensityWetting(i+ii,j+jj,k+kk,{m})/dRho, lbTemperature[m], CriticalDensity[m]/dRho, lbCriticalTemperature[m]);

                // Calculate interaction potential
                const double Phi_000 = Phi(DensityWetting(i,   j,   k   ,{n})/dRho, p_000, GasParameter[n]);
                const double Phi_ijk = Phi(DensityWetting(i+ii,j+jj,k+kk,{m})/dRho, p_ijk, GasParameter[m]);

                // Interpolate between both way
                tmp = 6.0*lbWeights[ii+1][jj+1][kk+1]*lbGK[n][m]*(ParaKuper*(Phi_ijk*Phi_ijk) + (1.0-2.0*ParaKuper)*(Phi_000*Phi_ijk));
            }

            const dVector3 BenziForceDensity = {tmp*ii*df, tmp*jj*df, tmp*kk*df};

            ForceDensity(i,j,k,{n}) += BenziForceDensity;

            if (!Obstacle(i,j,k) && Obstacle(i+ii,j+jj,k+kk))
            for(auto& it : Phase.Fields(i+ii,j+jj,k+kk))
            {
                Grain& grain = Phase.FieldsProperties[it.index];
                if(grain.State == AggregateStates::Solid)
                {
                    const dVector3 pos = {double(i+ii), double(j+jj), double(k+kk)};
                    dVector3 distanceCM;
                    CalculateDistancePeriodic(pos, grain.Rcm, distanceCM, Grid.Nx, Grid.Ny, Grid.Nz);
                    const dVector3 R = distanceCM * Grid.dx;

                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    {
                        grain.Force  -= BenziForceDensity * Grid.CellVolume(true);
                        grain.Torque -= R.cross(BenziForceDensity) * Grid.CellVolume(true);
                    }
                }
            }
        }
    }
}

double FlowSolverLBM::CalculateLocalIncomingSolidDensityChange(PhaseField& Phase, int i, int j, int k)
{
    double locIncomingDensityChange = 0.0;
    NodeA<double> locPhiDot = Phase.Dot(i,j,k,dt);
    for(auto alpha = locPhiDot.cbegin(); alpha < locPhiDot.cend(); ++alpha)
    {
        if (Phase.FieldsProperties[alpha->index].State  == AggregateStates::Solid)
        {
            locIncomingDensityChange += alpha->value*Phase.FieldsProperties[alpha->index].Density*dt;
        }
    }

    return locIncomingDensityChange;
}

void FlowSolverLBM::AcountForAndLimitPhaseTransformation(PhaseField& Phase, int tStep)
{
    // This method calculates the fluid density change according to phase-
    // transformations like solidification. Further it limits the
    // phase-transformation if not enough liquid is present.
    // It only needs to be calculated if the fluid density is assumed to be
    // inhomogeneous.
  
    if (tStep == 0)
    {
        ConsoleOutput::WriteWarning("Simulation of a two phase fluid and Solidification and Melting is still experimental", thisclassname, "AcountForAndLimitPhaseTransformation");
    }

    if (Do_TwoPhase)
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    if(Phase.Fields(i,j,k).wide_interface())
    {
        double locDensity = DensityWetting(i,j,k,{0});
        double locIncomingDensityChange = CalculateLocalIncomingSolidDensityChange(Phase,i,j,k);
        double locDensityChange = locIncomingDensityChange;

        // Limit incoming density change
        const double locDensityMinium = 2*VaporDensity[0];
        if (locIncomingDensityChange > 0.0) // TODO also consider melting
        {
            double Limiter = 1.0;
            if (locDensity < locDensityMinium)
            {
                locDensityChange = 0.0;
                Limiter = 0.0;
            }
            else if (locDensity - locIncomingDensityChange < locDensityMinium)
            {
                locDensityChange = locDensityMinium-locDensity;
                Limiter = locDensityChange/locIncomingDensityChange;
            }

            // Limit phase-transformation to ensure mass conservation
            for(auto& psi : Phase.FieldsDot(i,j,k))
            {
                const Grain GrainA = Phase.FieldsProperties[psi.indexA];
                const Grain GrainB = Phase.FieldsProperties[psi.indexB];
                if (psi.value1 != 0.0)
                if ((GrainA.is_solid() and GrainB.is_fluid()) or
                    (GrainB.is_solid() and GrainA.is_fluid()))
                {
                    psi.value1 *= Limiter;
                }
            }
        }

        // Apply density change
        DensityWetting(i,j,k,{0}) -= locDensityChange;
        lbPopulations (i,j,k,{0}) = EquilibriumDistribution(DensityWetting(i,j,k,{0})/dRho, lbWeights, MomentumDensity(i,j,k,{0})/dm);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CalculateForceGravity(PhaseField& Phase)
{
    // Calculate gravity on solid bodies
    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if (Phase.FieldsProperties[idx].State == AggregateStates::Solid and
        Phase.FieldsProperties[idx].Mobile)
    {
        Phase.FieldsProperties[idx].Acm += GA;
    }

    // Calculate gravity on fluid nodes
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0,)
    if (not Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ForceDensity(i,j,k,{n}) += GA * DensityWetting(i,j,k,{n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CalculateForceGravity(PhaseField& Phase,
        const Composition& Cx)
{

    // Calculate force correction arising form compositional dependent density.
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0,)
    for(auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); it++)
    {
        double locDeltaDensity = 0.0;
        const Grain& grain = Phase.FieldsProperties[it->index];
        if(grain.State == AggregateStates::Liquid)
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        for (size_t Comp = 0; Comp < Ncomp; ++Comp)
        {
            locDeltaDensity += dDensity_dc[Comp] * (Cx.MoleFractions(i,j,k,{grain.Phase, Comp}) - Cx.Initial({grain.Phase, Comp}));
            ForceDensity(i,j,k,{n}) += GA * locDeltaDensity;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CalculateDensityTC(Temperature& Tx, Composition& Cx, const BoundaryConditions& BC)
{
    double Ts=273.0;
    double mus=1.68e-05; //it's for air
    double S=110.5;
    double AW_air = Cx.AtomicWeightMixture;
    double R= PhysicalConstants::R;     ///<  Universal gas constant  j/mol k
    double Rm=R/AW_air;             ///< gas constant for air j/kg k

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        DensityWetting(i,j,k,{n}) = Pth/(Rm * Tx.Tx(i,j,k));
        double mu=mus*(Tx.Tx(i,j,k)/Ts)*sqrt(Tx.Tx(i,j,k)/Ts)*(Ts+S)/(Tx.Tx(i,j,k)+S);
        nut(i,j,k,{n}) = mu/DensityWetting(i,j,k,{n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double FlowSolverLBM::CalculateSutherlandViscosity(double p, double mw, double Temp)
{
    double Ts=273.0;
    double mus=1.68e-05; //it's for air
    double S=110.5;
    double mu=mus*(Temp/Ts)*sqrt(Temp/Ts)*(Ts+S)/(Temp+S);
    return mu;
}
double FlowSolverLBM::CalculateIdealGasDensity(double p, double Rm, double Temp)
{
    double rho= p/(Rm * Temp);
    return rho;
}

void FlowSolverLBM::SetInitialPopulationsTC(BoundaryConditions& BC, Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,HydroPressure,HydroPressure.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        HydroPressure (i,j,k,{n}) = Poutlet;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        dVector3 lbvel      = Vel.Phase(i,j,k,{n})/dv;
    	dVector3 lbFb       = ForceDensity(i,j,k,{n})/df;
    	dVector3 lbGradRho  = GradRho(i,j,k,{n})*Grid.dx;
        double lbnu         = nut(i,j,k,{n})/dnu;
        double lbrho        = DensityWetting(i,j,k,{n})/dRho;
        double lbph         = HydroPressure(i,j,k,{n})/dP;
        double lbDivVel     = DivVel(i,j,k,{n})*dt;
        lbPopulations(i,j,k,{n}) = EquilibriumDistributionTC( lbvel, lbph, lbnu,  lbDivVel,  lbFb, lbWeights, lbrho, lbGradRho, n);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    BC.SetXVector(lbPopulations);
    BC.SetYVector(lbPopulations);
    BC.SetZVector(lbPopulations);
}
void FlowSolverLBM::CalculateDensityGradient(PhaseField& Phase, Velocities& Vel)
{
    double eps=1e-10;
    std::vector<int> dir(3);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    {
        if(!Obstacle(i,j,k))
        for (int direction=0; direction<3; ++direction)
        {
            dir[0] = (direction == 0)?1*Grid.dNx:0;
            dir[1] = (direction == 1)?1*Grid.dNy:0;
            dir[2] = (direction == 2)?1*Grid.dNz:0;
            int ip=i+dir[0];
            int jp=j+dir[1];
            int kp=k+dir[2];
            int im=i-dir[0];
            int jm=j-dir[1];
            int km=k-dir[2];

            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            {
                double v = 0.0;
                for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); it++)
                if (Phase.FieldsProperties[it->index].State != AggregateStates::Solid)
                {
                    size_t PhaseIdx = Phase.FieldsProperties[it->index].Phase;
                    v = Vel.Phase(i,j,k,{PhaseIdx})[direction] ;
                }
                double Rho   = DensityWetting(i,j,k,{n});
                double Rhop  = DensityWetting(ip,jp,kp,{n});
                double Rhom  = DensityWetting(im,jm,km,{n});

                if(GradRho_Upwind)
                {
                    if(v>=0)
                    {
                        GradRho(i,j,k,{n})[direction] = (Rho-Rhom)/Grid.dx;
                    }
                    else
                    {
                        GradRho(i,j,k,{n})[direction] = (Rhop-Rho)/Grid.dx;
                    }
                }
                else if(GradRho_VanLeer)
                {
                    double Rhohp = (Rho+Rhop)/2.0;
                    double Rhohm = (Rho+Rhom)/2.0;
                    if(v>=0)
                    {
                        int imm=i-2.0*dir[0];
                        int jmm=j-2.0*dir[1];
                        int kmm=k-2.0*dir[2];

                        double Rhomm  = DensityWetting(imm,jmm,kmm,{n});
                        double r  = (Rho-Rhom)/(Rhop-Rho+eps);
                        double phi=(r+fabs(r))/(1.0+fabs(r));
                        Rhohp=Rho+0.5*phi*(Rhop-Rho);

                        double rm  = (Rhom-Rhomm)/(Rho-Rhom+eps);
                        double phim=(rm+fabs(rm))/(1.0+fabs(rm));
                        Rhohm=Rhom+0.5*phim*(Rho-Rhom);
                        GradRho(i,j,k,{n})[direction] =  (Rhohp-Rhohm)/Grid.dx;
                    }
                    else
                    {
                        int ipp=i+2.0*dir[0];
                        int jpp=j+2.0*dir[1];
                        int kpp=k+2.0*dir[2];
                        double Rhopp  = DensityWetting(ipp,jpp,kpp,{n});

                        double r  = (Rhopp-Rhop)/(Rhop-Rho+eps);
                        double phi=(r+fabs(r))/(1.0+fabs(r));
                        Rhohp=Rhop-0.5*phi*(Rhop-Rho);

                        double rm  = (Rhop-Rho)/(Rho-Rhom+eps);
                        double phim=(rm+fabs(rm))/(1.0+fabs(rm));
                        Rhohm=Rho-0.5*phim*(Rho-Rhom);
                        GradRho(i,j,k,{n})[direction] =  (Rhohp-Rhohm)/Grid.dx;
                    }
                }
                else if(GradRho_Central)
                {
                    GradRho(i,j,k,{n})[direction] =  (Rhop-Rhom)/Grid.dx/2.0;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void FlowSolverLBM::CalculateHydrodynamicPressureAndMomentum(Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if(!Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        HydroPressure(i,j,k,{n}) = 0.0;
        MomentumDensity(i,j,k,{n}).set_to_zero();

        for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
        for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
        for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
        {
            double rr = lbPopulations(i,j,k,{n})(ii,jj,kk);
            HydroPressure  (i,j,k,{n})    += rr*dP;
            MomentumDensity(i,j,k,{n})[0] += rr*ii;
            MomentumDensity(i,j,k,{n})[1] += rr*jj;
            MomentumDensity(i,j,k,{n})[2] += rr*kk;
        }
        
        MomentumDensity(i,j,k,{n}) = MomentumDensity(i,j,k,{n})/cs2*pow(Grid.dx,3)/pow(dt,3) + ForceDensity(i,j,k,{n})*dt/2.0;
        dVector3 vel = MomentumDensity(i,j,k,{n}) /DensityWetting(i,j,k,{n});
        HydroPressure (i,j,k,{n}) +=  cs2*dt/2.0 * ( (vel[0]*GradRho(i,j,k,{n})[0] + vel[1]*GradRho(i,j,k,{n})[1] + vel[2]*GradRho(i,j,k,{n})[2]) + DensityWetting(i,j,k,{n}) * DivVel(i,j,k,{n}));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CollisionTC( Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        dVector3 lbvel      = Vel.Phase(i,j,k,{n})/dv;
        dVector3 lbFb       = ForceDensity(i,j,k,{n})/df;
        dVector3 lbGradRho  = GradRho(i,j,k,{n})*Grid.dx;
        double lbu2         = lbvel[0]*lbvel[0]+lbvel[1]*lbvel[1]+lbvel[2]*lbvel[2];
        double lbnu         = nut(i,j,k,{n})/dnu;
        double lbrho        = DensityWetting(i,j,k,{n})/dRho;
        double lbph         = HydroPressure(i,j,k,{n})/dP;
        double lbDivVel     = DivVel(i,j,k,{n})*dt;
        lbtau [n]            = lbnu/lbcs2 + 0.5; 
        double factor       = (1.0-1.0/(2.0*lbtau[n]));

        for(int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
        for(int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
        for(int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
        {
            double CR = (ii-lbvel[0])*lbGradRho[0]+(jj-lbvel[1])*lbGradRho[1]+(kk-lbvel[2])*lbGradRho[2];
            double CF = (ii-lbvel[0])*lbFb[0]+(jj-lbvel[1])*lbFb[1]+(kk-lbvel[2])*lbFb[2];
            double lbw  = lbWeights[ii+1][jj+1][kk+1];
            double lbcu = ii*lbvel[0] + jj*lbvel[1] + kk*lbvel[2];
            double Feq  = lbrho * lbw *(1.0 + lbcu/lbcs2- lbu2/(2.0*lbcs2) + lbcu*lbcu/(2.0*lbcs2*lbcs2) );
            double geq  = Feq*lbcs2 + (lbph - lbrho * lbcs2 ) * lbw; 

            double geqb = geq - 0.5 * (lbrho * lbcs2 * lbw * lbDivVel * factor  + CR*(Feq/lbrho-lbw)*lbcs2*factor + CF*Feq/lbrho);
            double gb  = lbPopulations(i,j,k,{n})(ii,jj,kk);
            lbPopulations(i,j,k,{n})(ii,jj,kk) = gb + 1.0/lbtau[n] * (geqb - gb) + lbcs2 * lbw * lbrho * lbDivVel * factor + CR*(Feq/lbrho-lbw)*lbcs2*factor + CF*Feq/lbrho;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::SetVelocityInletUniform(const BoundaryConditions& BC, const PhaseField& Phase)
{

#ifdef MPI_PARALLEL

    if(MPI_CART_RANK[0]==0)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
        {
            if(i==0)
            {
                if(!Obstacle(i,j,k-1) and !Obstacle(i,j,k+1))
                {
                    for (size_t n = 0; n < N_Fluid_Comp; ++n)
                    {
                        if(Do_ThermalComp)
                        {
                        	for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                            for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                            for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                            {
                                double lbUx=U0X*dt/Grid.dx;
                                if(ii==1)
                                {
                                    double cu = ii*lbUx;
                                    lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                                                                +2.0*lbWeights[ii+1][jj+1][kk+1] * DensityWetting(i-ii,j,k,{n})/dRho * cu;
                                }
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
#endif
}

void FlowSolverLBM::SetCornerCorrection(const BoundaryConditions& BC, const PhaseField& Phase)
{

#ifdef MPI_PARALLEL
    if(MPI_CART_RANK[0]==0)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
        {
            if(i==0)
            {
                if(!Obstacle(i,j,k) and Obstacle(i,j,k-1))   //left corners
                {
                    double lbUx=U0X*dt/Grid.dx;
                    for (size_t n = 0; n < N_Fluid_Comp; ++n)
                    {
                        for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                        {

                            if(ii==1 and kk==0)
                            {
                                double cu = ii*lbUx;
                                lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                                                                    +2.0*lbWeights[ii+1][jj+1][kk+1] * DensityWetting(i-ii,j,k,{n})/dRho * cu;
                            }
                            if(ii==1 and kk==-1)
                            {
                                double cu = ii*lbUx;
                                lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                                        +2.0*lbWeights[ii+1][jj+1][kk+1] * DensityWetting(i-ii,j,k,{n})/dRho * cu;
                            }
                        }
                    }
                }
                if(!Obstacle(i,j,k) and Obstacle(i,j,k+1))   //right corners
                {
                    double lbUx=U0X*dt/Grid.dx;
                    for (size_t n = 0; n < N_Fluid_Comp; ++n)
                    {
                        for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                        {
                            if(ii==1 and kk==0)
                            {
                                double cu = ii*lbUx;
                                lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                                                                    +2.0*lbWeights[ii+1][jj+1][kk+1] * DensityWetting(i-ii,j,k,{n})/dRho * cu;
                            }
                            if(ii==1 and kk==1)
                            {
                                double cu = ii*lbUx;
                                lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                                        +2.0*lbWeights[ii+1][jj+1][kk+1] * DensityWetting(i-ii,j,k,{n})/dRho * cu;
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
#endif
}

void FlowSolverLBM::SetPressureOutlet(Velocities& Vel, const PhaseField& Phase, double KMaxv)
{
    if(Phase.Grid.OffsetX+Grid.Nx==Grid.TotalNx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        {
            if(i==Grid.Nx-1)
            {
                if(!Obstacle(i,j,k))
                for (size_t n = 0; n < N_Fluid_Comp; ++n)
                {
                    //===========================
                    double Pbold   = HydroPressure(i,j,k,{n});
                    double Pbmold  = HydroPressure(i-1,j,k,{n});
                    double Pbmmold = HydroPressure(i-2,j,k,{n});

                    double uxbold  = Vel.Phase(i,j,k,{n})[0];
                    double uxbmold = Vel.Phase(i-1,j,k,{n})[0];
                    double uxbmmold= Vel.Phase(i-2,j,k,{n})[0];

                    double uybold  = Vel.Phase(i,j,k,{n})[1];
                    double uybmold = Vel.Phase(i-1,j,k,{n})[1];
                    double uybmmold= Vel.Phase(i-2,j,k,{n})[1];

                    double uzbold  = Vel.Phase(i,j,k,{n})[2];
                    double uzbmold = Vel.Phase(i-1,j,k,{n})[2];
                    double uzbmmold= Vel.Phase(i-2,j,k,{n})[2];

                    double rhob = DensityWetting(i,j,k,{n});
                    double rhobm= DensityWetting(i-1,j,k,{n});

                    double L1,L3,L4, L5; //charachteristics value
                    double dx=Grid.dx;
                    L1= KMaxv  * (Pbold - Poutlet);
                    L3= uxbold * (uzbold-uzbmold)/dx;
                    L4= uxbold * (uybold-uybmold)/dx;
                    L5=(uxbold+sqrt(cs2)) * ((Pbold-Pbmold)/dx + rhob * sqrt(cs2) * (uxbold-uxbmold)/dx);

                    double uxbn = uxbold - (L5-L1)*dt/(2.0* sqrt(cs2) * rhob );
                    double uybn = uybold - L4 * dt;
                    double uzbn = uzbold - L3 * dt;
                    double pbn  = Pbold - (L5+L1)*dt/2.0;

                    if(uxbmold<0)
                    {
                        L3   = uxbmold * (uzbmold-uzbmmold)/dx;
                        L4   = uxbmold * (uybmold-uybmmold)/dx;
                        L5   =(uxbmold+sqrt(cs2)) * ((Pbmold-Pbmmold)/dx + rhobm * sqrt(cs2) * (uxbmold-uxbmmold)/dx);
                        pbn  = Pbmold -  (L5+L1)*dt/2.0;
                        uxbn = uxbmold - (L5-L1)*dt/(2.0* sqrt(cs2) * rhobm );
                        uzbn = uzbmold - L3 * dt;
                        uybn = uybmold - L4 * dt;
                    }
                    dVector3 lbub  = {uxbn*dt/dx, uybn*dt/dx, uzbn*dt/dx};

                    int im = i-1;
                    double pm = 0.0;
                    dVector3 Mvm = {0.0, 0.0, 0.0};
                    for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                    for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                    for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                    {
                        double rrm = lbPopulations(im,j,k,{n})(ii,jj,kk);
                        pm  += rrm*dP;
                        Mvm[0] += rrm*ii*dm;
                        Mvm[1] += rrm*jj*dm;
                        Mvm[2] += rrm*kk*dm;
                    }
                    dVector3 um = (Mvm * 3.0 + ForceDensity(im,j,k,{n})*dt/2.0)/rhobm;
                    pm +=  cs2/2.0 *dt* (  (um[0]*  GradRho(im,j,k,{n})[0] )
                                         +(um[1]*  GradRho(im,j,k,{n})[1] )
                                         +(um[2]*  GradRho(im,j,k,{n})[2] ) + rhobm * DivVel(i-1,j,k,{n}) );

                    dVector3 lbum = um *dt/Grid.dx;

                    double lbu2b  = lbub[0]*lbub[0]+lbub[1]*lbub[1]+lbub[2]*lbub[2];
                    double lbu2m  = lbum[0]*lbum[0]+lbum[1]*lbum[1]+lbum[2]*lbum[2];

                    //double lbtaub  = nut(i ,j,k,{n})/dnu/lbcs2 + 0.5;
                    //double lbtaum  = nut(im,j,k,{n})/dnu/lbcs2 + 0.5;
                    int ii=-1;
                    //for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                    {
                    for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                    {
                        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                        {
                            //if(!Obstacle(i-ii,j-jj,k-kk))
                            {
                            //dVector3 lbFbb       = ForceDensity(i,j,k,{n})/df;
                            //dVector3 lbFbm       = ForceDensity(im,j,k,{n})/df;
                            //dVector3 lbGradRhob  = GradRho(i,j,k,{n})*Grid.dx;
                            //dVector3 lbGradRhom  = GradRho(im,j,k,{n})*Grid.dx;
                            double lbrhob        = DensityWetting(i,j,k,{n})/dRho;
                            double lbrhom        = DensityWetting(im,j,k,{n})/dRho;
                            double lbphb         = pbn/dP;
                            double lbphm         = pm /dP;
                            //double lbDivVelb     = DivVel(i,j,k,{n})*dt;
                            //double lbDivVelm     = DivVel(im,j,k,{n})*dt;
                            //double factorb       = (1.0+1.0/(2.0*lbtaub));
                            //double factorm       = (1.0+1.0/(2.0*lbtaum));
                            //double CRb           = (ii-lbub[0])*lbGradRhob[0]+(jj-lbub[1])*lbGradRhob[1]+(kk-lbub[2])*lbGradRhob[2];
                            //double CRm           = (ii-lbum[0])*lbGradRhom[0]+(jj-lbum[1])*lbGradRhom[1]+(kk-lbum[2])*lbGradRhom[2];
                            //double CFb           = (ii-lbub[0])*lbFbb[0]+(jj-lbub[1])*lbFbb[1]+(kk-lbub[2])*lbFbb[2];
                            //double CFm           = (ii-lbum[0])*lbFbm[0]+(jj-lbum[1])*lbFbm[1]+(kk-lbum[2])*lbFbm[2];
                            double lbw           = lbWeights[ii+1][jj+1][kk+1];
                            double lbcub         = ii*lbub[0] + jj*lbub[1] + kk*lbub[2];
                            double lbcum         = ii*lbum[0] + jj*lbum[1] + kk*lbum[2];
                            double Feqb          = lbrhob * lbw *(1.0 + lbcub/lbcs2- lbu2b/(2.0*lbcs2) + lbcub*lbcub/(2.0*lbcs2*lbcs2) );
                            double Feqm          = lbrhom * lbw *(1.0 + lbcum/lbcs2- lbu2m/(2.0*lbcs2) + lbcum*lbcum/(2.0*lbcs2*lbcs2) );
                            double geqb          = Feqb*lbcs2 + (lbphb - lbrhob * lbcs2 ) * lbw; 
                            double geqm          = Feqm*lbcs2 + (lbphm - lbrhom * lbcs2 ) * lbw; 

                            //double Ab            = 0.5 * (lbw*lbcs2*lbrhob*lbDivVelb+CRb*lbcs2*(Feqb/lbrhob-lbw)+ CFb*Feqb/lbrhob);
                            //double Am            = 0.5 * (lbw*lbcs2*lbrhom*lbDivVelm+CRm*lbcs2*(Feqm/lbrhom-lbw)+ CFm*Feqm/lbrhom);

                            double gbarm         = lbPopulations(im,j,k,{n})(ii,jj,kk);
                            //double gm            = (gbarm+1/lbtaum/2.0*geqm+Am)/factorm;
                            //double gb            = geqb + gm - geqm;
                            //lbPopulations(i,j,k,{n})(ii,jj,kk) = gb - 1/2.0/lbtaub * (geqb-gb) - Ab;
                            lbPopulations(i,j,k,{n})(ii,jj,kk) = geqb + gbarm - geqm;
                            }
                        }
                    }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void FlowSolverLBM::Set1DNonReflectingPressureOutletatBCXN(Velocities& Vel, const PhaseField& Phase)
{
#ifdef MPI_PARALLEL
    if(MPI_CART_RANK[0] == MPI_CART_SIZE[0]-1)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        {
            if(i==Grid.Nx-1)
            {
                for (size_t n = 0; n < N_Fluid_Comp; ++n)
                {
                    //===========================
                    double Pbold   = HydroPressure(i,j,k,{n});
                    double Pbmold  = HydroPressure(i-1,j,k,{n});
                    double Pbmmold = HydroPressure(i-2,j,k,{n});

                    double uxbold  = Vel.Phase(i,j,k,{n})[0];
                    double uxbmold = Vel.Phase(i-1,j,k,{n})[0];
                    double uxbmmold= Vel.Phase(i-2,j,k,{n})[0];

                    double uzbold  = Vel.Phase(i,j,k,{n})[2];
                    double uzbmold = Vel.Phase(i-1,j,k,{n})[2];
                    double uzbmmold= Vel.Phase(i-2,j,k,{n})[2];

                    double rhob = DensityWetting(i,j,k,{n});
                    double rhobm= DensityWetting(i-1,j,k,{n});

                    double L1, L3, L5; //charachteristics value
                    double KMaxv = 0.2;
                    L1=KMaxv * (Pbold - Poutlet);
                    L3= uxbold * (uzbold-uzbmold)/Grid.dx;
                    L5=(uxbold+sqrt(cs2)) * ((Pbold-Pbmold)/Grid.dx + rhob * sqrt(cs2) * (uxbold-uxbmold)/Grid.dx);

                    double uxbn = uxbold - (L5-L1)*dt/(2.0* sqrt(cs2) * rhob );
                    double uzbn = uzbold - L3 * dt;
                    double pb   = Pbold - (L5+L1)*dt/2.0;

                    if(uxbmold<0)
                    {
                        L3= uxbmold * (uzbmold-uzbmmold)/Grid.dx ;
                        L5=(uxbmold+sqrt(cs2)) * ((Pbmold-Pbmmold)/Grid.dx + rhobm * sqrt(cs2) * (uxbmold-uxbmmold)/Grid.dx);
                        pb = Pbmold - (L5+L1)*dt/2.0;
                        uxbn = uxbmold - (L5-L1)*dt/(2.0* sqrt(cs2) * rhobm );
                        uzbn = uzbmold - L3 * dt;
                    }
                    dVector3 lbub  = {uxbn*dt/Grid.dx, 0.0, uzbn*dt/Grid.dx};

                    int im = i-1;
                    double pm = 0.0;
                    dVector3 Mvm = {0.0, 0.0, 0.0};
                    for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                    for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                    for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                    {
                        double rrm = lbPopulations(im,j,k,{n})(ii,jj,kk);
                        pm  += rrm*dP;
                        Mvm[0] += rrm*ii*dm;
                        Mvm[1] += rrm*jj*dm;
                        Mvm[2] += rrm*kk*dm;
                    }
                    dVector3 um = (Mvm * 3.0 + ForceDensity(im,j,k,{n})*dt/2.0)/rhobm;
                    pm =  cs2/2.0 *dt* (  (um[0]*  GradRho(im,j,k,{n})[0] )
                                         +(um[1]*  GradRho(im,j,k,{n})[1] )
                                         +(um[2]*  GradRho(im,j,k,{n})[2] ) + rhobm * DivVel(i-1,j,k,{n}) ) + pm;

                    dVector3 lbum = um *dt/Grid.dx;

                    double lbu2b  = lbub[0]*lbub[0]+lbub[1]*lbub[1]+lbub[2]*lbub[2];
                    double lbu2m  = lbum[0]*lbum[0]+lbum[1]*lbum[1]+lbum[2]*lbum[2];

                    double lbtaub  = nut(i,j,k,{n})/dnu/lbcs2 + 0.5;
                    double lbtaum  = nut(im,j,k,{n})/dnu/lbcs2 + 0.5;
                    int ii=-1;
                    for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                    {
                        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                        {
                            if(Obstacle(i,j,k-1) and kk==1)
                            {
                            }
                            else if(Obstacle(i,j,k+1) and kk==-1)
                            {

                            }
                            else
                            {
                                dVector3 lbFbb       = ForceDensity(i,j,k,{n})/df;
                                dVector3 lbFbm       = ForceDensity(im,j,k,{n})/df;
                                dVector3 lbGradRhob  = GradRho(i,j,k,{n})*Grid.dx;
                                dVector3 lbGradRhom  = GradRho(im,j,k,{n})*Grid.dx;
                                double lbrhob        = DensityWetting(i,j,k,{n})/dRho;
                                double lbrhom        = DensityWetting(im,j,k,{n})/dRho;
                                double lbphb         = pb/dP;
                                double lbphm         = pm/dP;
                                double lbDivVelb     = DivVel(i,j,k,{n})*dt;
                                double lbDivVelm     = DivVel(im,j,k,{n})*dt;
                                //double factorb       = (1.0+1.0/(2.0*lbtaub));
                                double factorm       = (1.0+1.0/(2.0*lbtaum));
                                double CRb           = (ii-lbub[0])*lbGradRhob[0]+(jj-lbub[1])*lbGradRhob[1]+(kk-lbub[2])*lbGradRhob[2];
                                double CRm           = (ii-lbum[0])*lbGradRhom[0]+(jj-lbum[1])*lbGradRhom[1]+(kk-lbum[2])*lbGradRhom[2];
                                double CFb           = (ii-lbub[0])*lbFbb[0]+(jj-lbub[1])*lbFbb[1]+(kk-lbub[2])*lbFbb[2];
                                double CFm           = (ii-lbum[0])*lbFbm[0]+(jj-lbum[1])*lbFbm[1]+(kk-lbum[2])*lbFbm[2];
                                double lbw           = lbWeights[ii+1][jj+1][kk+1];
                                double lbcub         = ii*lbub[0] + jj*lbub[1] + kk*lbub[2];
                                double lbcum         = ii*lbum[0] + jj*lbum[1] + kk*lbum[2];
                                double Feqb          = lbrhob * lbw *(1.0 + lbcub/lbcs2- lbu2b/(2.0*lbcs2) + lbcub*lbcub/(2.0*lbcs2*lbcs2) );
                                double Feqm          = lbrhom * lbw *(1.0 + lbcum/lbcs2- lbu2m/(2.0*lbcs2) + lbcum*lbcum/(2.0*lbcs2*lbcs2) );
                                double geqb          = Feqb*lbcs2 + (lbphb - lbrhob * lbcs2 ) * lbw; 
                                double geqm          = Feqm*lbcs2 + (lbphm - lbrhom * lbcs2 ) * lbw; 
                                double Ab            = 0.5 * (lbw*lbcs2*lbrhob*lbDivVelb+CRb*lbcs2*(Feqb/lbrhob-lbw)+ CFb*Feqb/lbrhob);
                                double Am            = 0.5 * (lbw*lbcs2*lbrhom*lbDivVelm+CRm*lbcs2*(Feqm/lbrhom-lbw)+ CFm*Feqm/lbrhom);
                                double gbm           = lbPopulations(im,j,k,{n})(ii,jj,kk);
                                double gm            = (gbm+1/lbtaum/2.0*geqm+Am)/factorm;
                                double gb            = geqb + gm - geqm;
                                lbPopulations(i,j,k,{n})(ii,jj,kk) = gb - 1/2.0/lbtaub * (geqb-gm) - Ab;
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
#endif
}

void FlowSolverLBM::SetOutFlow(int order)
{
    int i1 = Grid.Nx-1;
    int i2 = Grid.Nx + Grid.Bcells; 
    int j1 = -Grid.Bcells * Grid.dNy;
    int j2 = Grid.Ny + Grid.Bcells * Grid.dNy; 
    int k1 = -Grid.Bcells * Grid.dNz;
    int k2 = Grid.Nz + Grid.Bcells * Grid.dNz; 

    if(Grid.OffsetX+Grid.Nx==Grid.TotalNx)
    for (int i = i1; i < i2; ++i)
    {
        for (int j = j1; j< j2 ; ++j)
        {
            for (int k = k1; k< k2 ; ++k)
            {
                int ii=-1;
                for (size_t n = 0; n < N_Fluid_Comp; ++n)
                {
                    //for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                    {
                        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                        {
                            for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                            {
                                if(order == 0)
                                {
                                    lbPopulations(i,j,k,{n})(ii,jj,kk) = lbPopulations(i-1,j,k,{n})(ii,jj,kk);
                                }
                                else if (order == 1)
                                {
                                    lbPopulations(i,j,k,{n})(ii,jj,kk) = 2.0*lbPopulations(i-1,j,k,{n})(ii,jj,kk) 
                                                            - lbPopulations(i-2,j,k,{n})(ii,jj,kk) ;
                                }
                            } 
                        }
                    }
                }
            }
        }
    }
}

void FlowSolverLBM::CalculateForceDrag(const int i,const int j,const int k,
        PhaseField& Phase, const Velocities& Vel)
{
    double dx3 = Grid.CellVolume(true);

    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    if(DensityWetting(i,j,k,{n}) != 0.0)
    {
        double locSolidFraction = SolidFraction(i,j,k,Phase);
        const double mu = nu[n]*DensityWetting(i,j,k,{n});  //Dynamic viscosity [kg/m/s]
        for(auto& it : Phase.Fields(i,j,k))
        {
            Grain& grain = Phase.FieldsProperties[it.index];
            if(grain.State == AggregateStates::Solid)
            {
                const double factor =
                    (it.value * it.value * (1.0 - locSolidFraction))/
                    (Phase.Grid.iWidth * Phase.Grid.iWidth);
                const dVector3 locSolidVelocity = Vel.Phase(i,j,k,{grain.Phase});
                const dVector3 DragForceDensity = (FluidVelocity(i,j,k,n) - locSolidVelocity) * mu * factor * h_star;
                const dVector3 pos = {double(i), double(j), double(k)};
                dVector3 distanceCM;
                CalculateDistancePeriodic(pos, grain.Rcm, distanceCM, Grid.Nx, Grid.Ny, Grid.Nz);
                const dVector3 locR = distanceCM * Grid.dx;
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    grain.Force  += DragForceDensity * dx3;
                    grain.Torque += locR.cross(DragForceDensity) * dx3;
                }
            }
        }
    }
}

double FlowSolverLBM::BounceBack(const int i, const int j, const int k,
        const int ii, const int jj, const int kk, const size_t n,
        PhaseField& Phase, const BoundaryConditions& BC, double& lbDensityChange)
{
    double dx3 = Grid.CellVolume(true);

    double NewPopulation = lbPopulations(i,j,k,{n})(-ii,-jj,-kk);
    if (Do_BounceBack)
    {
        for(auto it : Phase.Fields(i-ii,j-jj,k-kk))
        {
            Grain& grain = Phase.FieldsProperties[it.index];
            if(grain.State == AggregateStates::Solid)
            {
                dMatrix3x3 W;
                W(0,0) = 0.0;
                W(1,1) = 0.0;
                W(2,2) = 0.0;
                W(0,1) = -grain.aVel[2];
                W(0,2) =  grain.aVel[1];
                W(1,2) = -grain.aVel[0];
                W(1,0) =  grain.aVel[2];
                W(2,0) = -grain.aVel[1];
                W(2,1) =  grain.aVel[0];

                const dVector3 pos = {double(i-0.5*ii), double(j-0.5*jj), double(k-0.5*kk)};
                dVector3 distanceCM;
                CalculateDistancePeriodic(pos, grain.Rcm, distanceCM, Grid.Nx, Grid.Ny, Grid.Nz);
                const dVector3 locR = distanceCM * Grid.dx;
                const dVector3 Vel = grain.Vcm + W * locR;

                const double tmp = lbWeights[ii+1][jj+1][kk+1] *
                    DensityWetting(i,j,k,{n})/dRho * (Vel[0]*ii + Vel[1]*jj + Vel[2]*kk) * dt/Grid.dx;

                const double lbBBDensity =
                    2.0 * (lbPopulations(i,j,k,{n})(-ii,-jj,-kk) + 3.0 * tmp);

                NewPopulation += it.value * 6.0 * tmp;

                // NOTE: Bounce Back Momentum Density is now in physical units
                dVector3 lbBounceBackForceDensity;
                lbBounceBackForceDensity[0] = - ii * lbBBDensity * it.value;
                lbBounceBackForceDensity[1] = - jj * lbBBDensity * it.value;
                lbBounceBackForceDensity[2] = - kk * lbBBDensity * it.value;

                const dVector3 BounceBackForceDensity = lbBounceBackForceDensity * df;

                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    grain.Force  += BounceBackForceDensity * dx3;
                    grain.Torque += locR.cross(BounceBackForceDensity) * dx3;
                }
            }
        }
    }
    else if (Do_BounceBackElastic)
    {
        // Author: raphael.schiedung@rub.de
        const double rho1    = lbPopulations(i,j,k,{n})(-ii,-jj,-kk);
        const double ci      = std::sqrt(ii*ii+jj*jj+kk*kk);
        const double ici     = 1.0/ci;
        const double p1      = rho1*ci;
        const dVector3 nn    = {-ii*ici, -jj*ici, -kk*ici};
        const double Elastic = 1.00;  // 0.0 models fully inelastic impact

        NewPopulation = 0.0;
        for(auto it : Phase.Fields(i-ii,j-jj,k-kk))
        {
            Grain& grain = Phase.FieldsProperties[it.index];
            if(grain.State == AggregateStates::Solid)
            {
                dMatrix3x3 W;
                W(0,0) = 0.0;
                W(1,1) = 0.0;
                W(2,2) = 0.0;
                W(0,1) = -grain.aVel[2];
                W(0,2) =  grain.aVel[1];
                W(1,2) = -grain.aVel[0];
                W(1,0) =  grain.aVel[2];
                W(2,0) = -grain.aVel[1];
                W(2,1) =  grain.aVel[0];

                const dVector3 pos = {double(i-0.5*ii), double(j-0.5*jj), double(k-0.5*kk)};
                dVector3 distanceCM;
                CalculateDistancePeriodic(pos, grain.Rcm,distanceCM, Grid.Nx, Grid.Ny, Grid.Nz);
                const dVector3 locR   = distanceCM * Grid.dx;
                const dVector3 Vel = grain.Vcm + W * locR;

                const double v2   = Vel*nn*dt/Grid.dx;
                const double rho2 = grain.Density/dRho;
                const double p2   = rho2*v2;
                const double p1p  = (p1*(rho1-rho2*Elastic) + p2*(rho1+rho1*Elastic))/(rho1+rho2);
                const double p2p  = (p1*(rho2+rho2*Elastic) + p2*(rho2-rho1*Elastic))/(rho1+rho2);

                if (p1p < 0)
                {
                    NewPopulation -= it.value*p1p*ici;
                }
                else
                {
                    // The solid bounce back from the fluid (locally).
                    // Check grain Densities !
                    std::cerr << "Warning fluid momentum is too high!"
                              << "Check parameter (eg. Densities).\n";
                }

                dVector3 BounceBackForceDensity = {0,0,0};
                BounceBackForceDensity += nn*it.value*(p2p-p2)*df;
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    grain.Force  += BounceBackForceDensity * dx3;
                    grain.Torque += locR.cross(BounceBackForceDensity) * dx3;
                }
            }
        }
    }
    lbDensityChange += NewPopulation - lbPopulations(i,j,k,{n})(-ii,-jj,-kk);
    return NewPopulation;
}

void FlowSolverLBM::Propagation(PhaseField& Phase, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulationsTMP,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbPopulationsTMP(i,j,k,{n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        double lbDensityChange = 0.0;
        for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
        {
            if (Obstacle(i, j, k))
            {
                lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) =
                    lbPopulations(i,j,k,{n})(ii,jj,kk);
            }
            else if (Obstacle(i-ii, j-jj, k-kk))
            {
                lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) =
                    BounceBack(i, j, k, ii, jj, kk, n, Phase, BC, lbDensityChange);
            }
            else
            {
                lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) =
                    lbPopulations(i-ii, j-jj, k-kk,{n})(ii,jj,kk);
            }
        }
        lbPopulationsTMP(i,j,k,{n})(0,0,0) -= lbDensityChange; //NoSlip
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbPopulations(i,j,k,{n}) = lbPopulationsTMP(i,j,k,{n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END

   // if(dNx) BC.SetXVector(lbPopulations);
   // if(dNy) BC.SetYVector(lbPopulations);
   // if(dNz) BC.SetZVector(lbPopulations);
}

double FlowSolverLBM::SecondOrderBounceBack(const int i, const int j, const int k,
        const int ii, const int jj, const int kk, const size_t n,
        PhaseField& Phase, const BoundaryConditions& BC, double& lbDensityChange)
{
    double NewPopulation = 0.0;
    size_t Gasidx = 0;
    //size_t Solididx = 0;
    //double ep=1e-9;
    for(auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
    if(Phase.FieldsProperties[it->index].State == AggregateStates::Gas)
    {
        Gasidx = it->index;
    }
    dVector3 Norm{0.0,0.0,0.0};
    for (auto it = Phase.Fields(i-ii,j-jj,k-kk).cbegin(); it != Phase.Fields(i-ii,j-jj,k-kk).cend(); ++it)
    {
    	if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
    	{
            //Solididx=it->index;
            NodeAB<dVector3,dVector3> Normals = Phase.Normals(i, j, k);
    		//dVector3 Norm = Normals.get_asym1(Solididx, Gasidx);
    	}
    }

    double Phif   = Phase.Fractions(i,j,k,{Gasidx});
    double DX     = fabs(Phase.Grid.iWidth/Pi * asin(1.0-2.0*Phif));
    //double DX   = fabs(Phase.iWidth/Pi * acos(2.0*Phif-1));
    dVector3 X    = Norm*DX;
    dVector3 Xf   = {double(i), double(j), double(k)};
    dVector3 Xs   = {double(i-ii), double(j-jj), double(k-kk)};
    dVector3 Xsf;
    CalculateDistancePeriodic(Xf, Xs, Xsf, Grid.Nx, Grid.Ny, Grid.Nz);
    double l      = Xsf.length();
    dVector3 esf  = Xsf.normalize();
    double q      = X.abs() / (Norm * esf);
    q=q/l;

    if(q>1.0)
    {
        q=1.0;
    }

    if(q<0.5)
    {
        NewPopulation = 2.0*q*lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                    +(1.0-2.0*q)*lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk);
    }
    else if(q>=0.5)
    {
        NewPopulation = 1.0/2.0/q*lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                       +(2.0*q-1.0)/2.0/q*lbPopulations(i,j,k,{n})(ii,jj,kk);
    }
    
    /*
    if(q<0.5)
    {
        NewPopulation = q*(2.0*q+1.0)*lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                        +(1.0-4.0*q*q)*lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                        -q*(1.0-2.0*q)*lbPopulations(i+2.0*ii,j+2.0*jj,k+2.0*kk,{n})(-ii,-jj,-kk);
    }
    else if(q>=0.5)
    {
        NewPopulation = 1.0/(q*(2.0*q+1))*lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                        + (2.0*q-1.0)/q*lbPopulations(i,j,k,{n})(ii,jj,kk)
                        + (1.0-2.0*q)/(1.0+2.0*q)*lbPopulations(i+ii,j+jj,k+kk,{n})(ii,jj,kk);
    }
    */

    return NewPopulation;
}

void FlowSolverLBM::PropagationSecondOrderBB(PhaseField& Phase, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulationsTMP,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbPopulationsTMP(i,j,k,{n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        double lbDensityChange = 0.0;
        for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
        {
            if (Obstacle(i, j, k))
            {
                lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) =
                    lbPopulations(i,j,k,{n})(ii,jj,kk);
            }
            else if (Obstacle(i-ii, j-jj, k-kk))
            {
                lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) =
                    SecondOrderBounceBack(i, j, k, ii, jj, kk, n, Phase, BC, lbDensityChange);
            }
            else
            {
                lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) =
                    lbPopulations(i-ii, j-jj, k-kk,{n})(ii,jj,kk);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbPopulations(i,j,k,{n}) = lbPopulationsTMP(i,j,k,{n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END

   // if(dNx) BC.SetXVector(lbPopulations);
   // if(dNy) BC.SetYVector(lbPopulations);
   // if(dNz) BC.SetZVector(lbPopulations);
}

void FlowSolverLBM::ApplyForces(PhaseField& Phase, const Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ForceDensity(i,j,k,{n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    {
        if (Do_TwoPhase) CalculateForceTwoPhase(i,j,k, Phase);
        if (Do_Drag and Phase.Fields(i,j,k).interface())
        {
            CalculateForceDrag(i,j,k, Phase, Vel);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (Do_Gravity) CalculateForceGravity(Phase);
}

void FlowSolverLBM::ApplyForces(PhaseField& Phase, const Velocities& Vel,
        const Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ForceDensity(i,j,k,{n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    {
        if (Do_TwoPhase) CalculateForceTwoPhase(i,j,k, Phase);
        if (Do_Drag and Phase.Fields(i,j,k).interface())
        {
            CalculateForceDrag(i,j,k, Phase, Vel);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (Do_Gravity)
    {
        CalculateForceGravity(Phase);
        CalculateForceGravity(Phase, Cx);
    }
}

void FlowSolverLBM::Collision()
{
    // Apply body force
    if (Do_GuoForcing) // Guo, Zheng, and, Shi, Phy. Rev. E (2002)
    {
        // Apply body force with force projection method
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (!Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const dVector3 EqMomentum = MomentumDensity(i,j,k,{n}) + ForceDensity(i,j,k,{n})*dt/2;
            lbPopulations(i,j,k,{n}) += (EquilibriumDistribution(DensityWetting(i,j,k,{n})/dRho, lbWeights, EqMomentum/dm) - lbPopulations(i,j,k,{n}))/lbtau[n];
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // Apply body force with force projection method
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (!Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            assert(DensityWetting(i,j,k,{n})/dRho > DBL_EPSILON);
            assert(DensityWetting(i,j,k,{n}) < DBL_MAX);

            dVector3 lbVel;
            if (Do_ThermalComp)
            {
                lbVel = (MomentumDensity(i,j,k,{n}) * 3.0 + ForceDensity(i,j,k,{n})*dt/2)/DensityWetting(i,j,k,{n})*dt/Grid.dx;
            }
            else
            {
                lbVel = (MomentumDensity(i,j,k,{n}) + ForceDensity(i,j,k,{n})*dt/2)/DensityWetting(i,j,k,{n})*dt/Grid.dx;
            }

            D3Q27 lblocForcing;
            lblocForcing.set_to_zero();
            for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
            for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
            for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
            {
                lblocForcing(ii,jj,kk) =
                    (1.0-0.5/lbtau[n]) * lbWeights[ii+1][jj+1][kk+1]*
                    (3.0*((ii-lbVel[0])*ForceDensity(i,j,k,{n})[0]/df +
                          (jj-lbVel[1])*ForceDensity(i,j,k,{n})[1]/df +
                          (kk-lbVel[2])*ForceDensity(i,j,k,{n})[2]/df) +
                     9.0*(lbVel[0]*ii + lbVel[1]*jj + lbVel[2]*kk)*
                     (ii*ForceDensity(i,j,k,{n})[0]/df +
                      jj*ForceDensity(i,j,k,{n})[1]/df +
                      kk*ForceDensity(i,j,k,{n})[2]/df));
            }

            lbPopulations(i,j,k,{n}) += lblocForcing;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else if (Do_EDForcing)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (!Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const double q = 1.0 - 1.0/lbtau[n];
            lbPopulations(i,j,k,{n}) = (lbPopulations(i,j,k,{n}) -
                                        EquilibriumDistribution(DensityWetting(i,j,k,{n})/dRho, lbWeights,  MomentumDensity(i,j,k,{n})/dm))*q +
                                        EquilibriumDistribution(DensityWetting(i,j,k,{n})/dRho, lbWeights, (MomentumDensity(i,j,k,{n})/dm + ForceDensity(i,j,k,{n})/df*dt));
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

size_t FlowSolverLBM::CountObstacleNodes(void) const
{
    size_t locObstacleNodes = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0, reduction(+:locObstacleNodes))
    {
        if (Obstacle(i,j,k))
        {
            locObstacleNodes++;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locObstacleNodes, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    return locObstacleNodes;
}

size_t FlowSolverLBM::CountFluidNodes(void) const
{
    size_t locFluidNodes = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0, reduction(+:locFluidNodes))
    {
        if (not Obstacle(i,j,k)) locFluidNodes++;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locFluidNodes, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    return locFluidNodes;
}

void FlowSolverLBM::EnforceMassConservation(void)
{
    // Calculate Fluid mass
    std::vector<double> DeltaDensity = CalculateFluidMass();
    size_t FluidNodes = CountFluidNodes();

    if (FluidNodes > 0)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    if (FluidMass[n] != 0)
    {
        DeltaDensity[n] -= FluidMass[n];
        DeltaDensity[n] /= FluidNodes;

        // Fix density
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        {
            if (not Obstacle(i,j,k))
            {
                lbPopulations(i,j,k,{n}) = lbPopulations(i,j,k,{n}) -
                    EquilibriumDistribution(DensityWetting(i,j,k,{n})/dRho + DeltaDensity[n], lbWeights) +
                    EquilibriumDistribution(DensityWetting(i,j,k,{n})/dRho, lbWeights);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    FluidMass = CalculateFluidMass();
}

void FlowSolverLBM::EnforceTotalMassConservation(PhaseField& Phase)
{
    if (N_Fluid_Comp > 1)
    {
        ConsoleOutput::WriteExit("Not yet implemented for multi fluid", thisclassname, "EnforceTotalMassConservation");
        OP_Exit(EXIT_FAILURE);
    }

    double FluidMass = CalculateFluidMass()[0];
    double SolidMass = Phase.CalculateSolidMass();
    double TotalMass = FluidMass + SolidMass;

    if (InitialTotalMass == -1.0) InitialTotalMass = TotalMass;
    double DeltaMass = TotalMass - InitialTotalMass;
    size_t FluidNodes = CountFluidNodes();
    double DensityDelta = DeltaMass/N_Fluid_Comp/FluidNodes/Grid.CellVolume(true);

    if (DeltaMass > std::numeric_limits<double>::epsilon())
    {
        if (FluidNodes > 0)
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        {
            if (not Obstacle(i,j,k))
            {
                lbPopulations(i,j,k,{0})(0,0,0) -= DensityDelta/dRho;
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void FlowSolverLBM::EnforceSolidMomentum(PhaseField& Phase, const dVector3 value)
{
    // Fix Solid momentum
    const double dx = Phase.Grid.dx;
    const double dV = dx*dx*dx;

    dVector3 SolidMomentum = {0.0,0.0,0.0};
    double SolidMass = 0.0;
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if (Phase.FieldsProperties[idx].State == AggregateStates::Solid and
        Phase.FieldsProperties[idx].Mobile)
    {
        const double GrainMass = Phase.FieldsProperties[idx].Volume*dV*
                                 Phase.FieldsProperties[idx].Density;
        SolidMass     += GrainMass;
        SolidMomentum += Phase.FieldsProperties[idx].Vcm * GrainMass;
    }

    const dVector3 VcmFix = (SolidMomentum - value)/SolidMass;

    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if (Phase.FieldsProperties[idx].State == AggregateStates::Solid and
        Phase.FieldsProperties[idx].Mobile)
    {
        Phase.FieldsProperties[idx].Vcm -= VcmFix;
    }
}

bool FlowSolverLBM::SingleSolid(const int i, const int j, const int k, const PhaseField& Phase) const
{
    for(auto it = Phase.Fields(i,j,k).cbegin();
            it != Phase.Fields(i,j,k).cend(); ++it)
    if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
    {
        if(it->value >= 0.95) return true;
    }
    return false;
}

double FlowSolverLBM::SolidFraction(const int i, const int j, const int k, const PhaseField& Phase) const
{
    double SolidFraction = 0.0;
    for(auto it = Phase.Fields(i,j,k).cbegin();
            it != Phase.Fields(i,j,k).cend(); ++it)
    if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
    {
        SolidFraction += it->value;
    }
    return SolidFraction;
}

bool FlowSolverLBM::LocalObstacle(const int i, const int j, const int k, const PhaseField& Phase) const
{
    //return SingleSolid(i,j,k,Phase);
    return (SolidFraction(i,j,k,Phase) >= 0.95);
}

void FlowSolverLBM::DetectObstacles(const PhaseField& Phase)
{
    bool locObstaclesChanged = false;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(), reduction(or:locObstaclesChanged))
    {
        if (LocalObstacle(i,j,k,Phase))
        {
            if (not Obstacle(i,j,k)) locObstaclesChanged = true;
            Obstacle(i,j,k) = true;
        }
        else
        {
            if (Obstacle(i,j,k)) locObstaclesChanged = true;
            Obstacle(i,j,k) = false;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locObstaclesChanged, 1, OP_MPI_CXX_BOOL, OP_MPI_LOR, OP_MPI_COMM_WORLD);
#endif
    ObstaclesChanged = locObstaclesChanged;
}

void FlowSolverLBM::DetectObstaclesAdvection(const PhaseField& Phase,
        const Velocities &Vel, const BoundaryConditions& BC)
{
#ifdef DEBUG
    dVector3 FluidMomentumOld = CalculateFluidMomentum(Phase)[0];
    double FluidMassOld = CalculateFluidMass()[0];
#endif

    size_t locObstaclesAppeared = 0;
    size_t locObstaclesVanished = 0;

    bool locObstaclesChanged = false;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(), reduction (+:locObstaclesAppeared,locObstaclesVanished) reduction(or:locObstaclesChanged))
    {
        ObstacleVanished(i,j,k) = false;
        ObstacleAppeared(i,j,k) = false;
        if (LocalObstacle(i,j,k,Phase))
        {
            if(not Obstacle(i,j,k))
            {
                ObstacleAppeared(i,j,k) = true;
                locObstaclesChanged = true;
                locObstaclesAppeared++;
            }
            Obstacle(i,j,k) = true;
        }
        else
        {
            if(Obstacle(i,j,k))
            {
                ObstacleVanished(i,j,k) = true;
                locObstaclesChanged = true;
                locObstaclesVanished++;
            }
            Obstacle(i,j,k) = false;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locObstaclesAppeared, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locObstaclesVanished, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locObstaclesChanged, 1, OP_MPI_CXX_BOOL, OP_MPI_LOR, OP_MPI_COMM_WORLD);
#endif

    ObstaclesChanged = locObstaclesChanged;

    if (locObstaclesAppeared + locObstaclesVanished > 0)
    {
        SetFluidNodesNearObstacle();
        SetObstacleNodes(Phase,Vel);
        CalculateDensityAndMomentum();
        SetBoundaryConditions(BC);

#ifdef DEBUG
        dVector3 FluidMomentumNew = CalculateFluidMomentum(Phase)[0];
        dVector3 FluidMomentumDelta = FluidMomentumNew-FluidMomentumOld;
        double FluidMomentumDeltaMag = FluidMomentumDelta.abs();
        double FluidMassNew = CalculateFluidMass()[0];
        double FluidMassDelta = FluidMassNew - FluidMassOld;
        if (std::abs(FluidMassDelta)/dRho > 1.0e-10)//DBL_EPSILON)
        {
            ConsoleOutput::WriteWarning("Mass conservation violated! Method needs to be debugged.", thisclassname, "DetectObstaclesAdvection");
            ConsoleOutput::Write("Obstacles Appeared", locObstaclesAppeared);
            ConsoleOutput::Write("Obstacles Vanished", locObstaclesVanished);
            ConsoleOutput::Write("Flud Mass Change [Kg]", FluidMassDelta);
            ConsoleOutput::Write("Flud Mass Change  [1]", FluidMassDelta/dRho);
        }
        if (FluidMomentumDeltaMag/dm/Grid.CellVolume(true) > 1.0e-10)//DBL_EPSILON)
        {
            ConsoleOutput::WriteWarning("Momentum conservation violated! Method may needs to be debugged.", thisclassname, "DetectObstaclesAdvection");
            ConsoleOutput::Write("Obstacles Appeared",locObstaclesAppeared);
            ConsoleOutput::Write("Obstacles Vanished",locObstaclesVanished);
            ConsoleOutput::Write("Fluid Momentum Change [Kg m/s]",FluidMomentumDeltaMag);
            ConsoleOutput::Write("Fluid Momentum Change      [1]",FluidMomentumDeltaMag/dm/Grid.CellVolume(true));
        }
#endif
    }
}

std::vector<dVector3> FlowSolverLBM::CalculateFluidMomentum(const PhaseField& Phase) const
{
    std::vector<dVector3> Momentum(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::vector<double> tmpMomentumX;
        std::vector<double> tmpMomentumY;
        std::vector<double> tmpMomentumZ;

        // Calculate fluid momentum
        #pragma omp parallel
        {
            std::vector<double> tmplocMomentumX;
            std::vector<double> tmplocMomentumY;
            std::vector<double> tmplocMomentumZ;
            #pragma omp for collapse(OMP_COLLAPSE_LOOPS)\
            schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
            for (int i = 0; i < Grid.Nx; i++)
            for (int j = 0; j < Grid.Ny; j++)
            for (int k = 0; k < Grid.Nz; k++)
            if (not Obstacle(i,j,k))
            for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
            for (int kk = -Grid.dNz; kk <= Grid.dNy; kk++)
            {
                double rr = lbPopulations(i,j,k,{n})(ii,jj,kk);
                tmplocMomentumX.push_back(rr*ii);
                tmplocMomentumY.push_back(rr*jj);
                tmplocMomentumZ.push_back(rr*kk);
            }

            std::sort(tmplocMomentumX.begin(),tmplocMomentumX.end(), [] (double a, double b) { return a < b;});
            std::sort(tmplocMomentumY.begin(),tmplocMomentumY.end(), [] (double a, double b) { return a < b;});
            std::sort(tmplocMomentumZ.begin(),tmplocMomentumZ.end(), [] (double a, double b) { return a < b;});

            double locMomentumX = 0.0;
            double locMomentumY = 0.0;
            double locMomentumZ = 0.0;

            for (auto value : tmplocMomentumX) locMomentumX += value;
            for (auto value : tmplocMomentumY) locMomentumY += value;
            for (auto value : tmplocMomentumZ) locMomentumZ += value;

            #pragma omp critical
            {
                tmpMomentumX.push_back(locMomentumX);
                tmpMomentumY.push_back(locMomentumY);
                tmpMomentumZ.push_back(locMomentumZ);
            }
        }

        std::sort(tmpMomentumX.begin(),tmpMomentumX.end(), [] (double a, double b) { return a < b;});
        std::sort(tmpMomentumY.begin(),tmpMomentumY.end(), [] (double a, double b) { return a < b;});
        std::sort(tmpMomentumZ.begin(),tmpMomentumZ.end(), [] (double a, double b) { return a < b;});

        for (auto value : tmpMomentumX) Momentum[n][0] += value;
        for (auto value : tmpMomentumY) Momentum[n][1] += value;
        for (auto value : tmpMomentumZ) Momentum[n][2] += value;

        #ifdef MPI_PARALLEL
        double tmpX = Momentum[n][0];
        double tmpY = Momentum[n][1];
        double tmpZ = Momentum[n][2];
        OP_MPI_Allreduce(&tmpX, &(Momentum[n][0]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&tmpY, &(Momentum[n][1]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&tmpZ, &(Momentum[n][2]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        Momentum[n] *= Grid.dx /dt *dRho*Grid.CellVolume(true);
    }
    return Momentum;
}

void FlowSolverLBM::SetFluidNodesNearObstacle()
{
    const int rrrr = std::max(3,FluidRedistributionRange*FluidRedistributionRange);
    auto InsideRedistributionRange = [rrrr](int ii, int jj, int kk) {return ii*ii + jj*jj + kk*kk <= rrrr;};
    auto redistributable = [this, InsideRedistributionRange](int i, int j, int k, int ii, int jj, int kk)
    {
        if (not Obstacle(i+ii,j+jj,k+kk))
        if (not ObstacleVanished(i+ii,j+jj,k+kk))
        if (InsideRedistributionRange(ii,jj,kk))
        if (ii + jj + kk != 0)
        {
            return true;
        }

        return false;
    };

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-FluidRedistributionRange,)
    {
        ObstacleChangedDensity(i,j,k) = false;
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            lbPopulationsTMP(i,j,k,{n}).set_to_zero();
            if (ObstacleVanished(i,j,k)) lbPopulations(i,j,k,{n}).set_to_zero();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-FluidRedistributionRange,)
    if (ObstacleVanished(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        // If an obstacle vanishes the present cell has to be filled with the
        // surrounding fluid from the neighbour cells

        // Calculate local average population
        double   SumWeights = 0;
        double   avDensity  = 0;
        dVector3 avMomentum = {0.0,0.0,0.0};
        int FluidNeighbours = 0;
        for (int ii = -FluidRedistributionRange*Grid.dNx; ii <= FluidRedistributionRange*Grid.dNx; ++ii)
        for (int jj = -FluidRedistributionRange*Grid.dNy; jj <= FluidRedistributionRange*Grid.dNy; ++jj)
        for (int kk = -FluidRedistributionRange*Grid.dNz; kk <= FluidRedistributionRange*Grid.dNz; ++kk)
        if (redistributable(i,j,k,ii,jj,kk))
        {
            const double weight = DensityWetting(i+ii,j+jj,k+kk,{n});
            SumWeights += weight;
            avDensity  += DensityWetting (i+ii,j+jj,k+kk,{n})*weight;
            avMomentum += MomentumDensity(i+ii,j+jj,k+kk,{n})*weight;
            FluidNeighbours++;
        }

        if (FluidNeighbours)
        {
            avDensity  /= SumWeights;
            avMomentum /= SumWeights;

            DensityWetting  (i,j,k,{n}) = avDensity;
            MomentumDensity (i,j,k,{n}) = avMomentum;
            const D3Q27 DeltaPop = EquilibriumDistribution(avDensity/dRho, lbWeights, avMomentum/dm);
            lbPopulationsTMP(i,j,k,{n}) += DeltaPop;
            ObstacleChangedDensity(i,j,k) = true;

            // Subtract added fluid from neighbouring cells
            for (int ii = -FluidRedistributionRange*Grid.dNx; ii <= FluidRedistributionRange*Grid.dNx; ++ii)
            for (int jj = -FluidRedistributionRange*Grid.dNy; jj <= FluidRedistributionRange*Grid.dNy; ++jj)
            for (int kk = -FluidRedistributionRange*Grid.dNz; kk <= FluidRedistributionRange*Grid.dNz; ++kk)
            if (redistributable(i,j,k,ii,jj,kk))
            {
                #pragma omp critical
                {
                    const double weight = DensityWetting(i+ii,j+jj,k+kk,{n});
                    lbPopulationsTMP(i+ii,j+jj,k+kk,{n}) -= DeltaPop*weight/SumWeights;
                    ObstacleChangedDensity(i+ii,j+jj,k+kk) = true;
                }
            }
        }
        else std::cerr << "Warning no fluid Neighbours\n";
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-FluidRedistributionRange,)
    if (ObstacleAppeared(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        // If an obstacle appears the present fluid has to be moved to the
        // neighbouring cell

        // Count number of fluid neighbours
        double   SumWeights = 0;
        int FluidNeighbours = 0;
        for (int ii = -FluidRedistributionRange*Grid.dNx; ii <= FluidRedistributionRange*Grid.dNx; ++ii)
        for (int jj = -FluidRedistributionRange*Grid.dNy; jj <= FluidRedistributionRange*Grid.dNy; ++jj)
        for (int kk = -FluidRedistributionRange*Grid.dNz; kk <= FluidRedistributionRange*Grid.dNz; ++kk)
        if (redistributable(i,j,k,ii,jj,kk))
        {
            const double weight = DensityWetting(i+ii,j+jj,k+kk,{n});
            SumWeights += weight;
            FluidNeighbours++;
        }

        if (FluidNeighbours)
        {
            // Add fluid from this cell to the neighbouring cells
            const D3Q27 DeltaPop = lbPopulations(i,j,k,{n});
            for (int ii = -FluidRedistributionRange*Grid.dNx; ii <= FluidRedistributionRange*Grid.dNx; ++ii)
            for (int jj = -FluidRedistributionRange*Grid.dNy; jj <= FluidRedistributionRange*Grid.dNy; ++jj)
            for (int kk = -FluidRedistributionRange*Grid.dNz; kk <= FluidRedistributionRange*Grid.dNz; ++kk)
            if (redistributable(i,j,k,ii,jj,kk))
            {
                #pragma omp critical
                {
                    const double weight = DensityWetting(i+ii,j+jj,k+kk,{n});
                    lbPopulationsTMP(i+ii,j+jj,k+kk,{n}) += DeltaPop*weight/SumWeights;
                    ObstacleChangedDensity(i+ii,j+jj,k+kk) = true;
                }
            }
       }
       else std::cerr << "Warning no fluid Neighbours\n";
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    // Add change in fluid nodes
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    {
        if (ObstacleChangedDensity(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            lbPopulations(i,j,k,{n}) += lbPopulationsTMP(i,j,k,{n});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::FixPopulations(void)
{
    // Add change in fluid nodes
    if (Do_FixPopulations)
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells(),)
    {
        if (not Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            bool FixPopulation = false;
            // Check for negative populations and fix them
            for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
            if (lbPopulations(i,j,k,{n})(ii,jj,kk) <= 0.0)
            {
                FixPopulation = true;
            }

            if (DensityWetting(i,j,k,{n}) <= 0.0)
            {
                FixPopulation = true;
                DensityWetting (i,j,k,{n}) = 0.0;
                MomentumDensity(i,j,k,{n}) = {0,0,0};
            }

            if (FixPopulation)
            {
                lbPopulations(i,j,k,{n}) =
                    EquilibriumDistribution(DensityWetting(i,j,k,{n})/dRho, lbWeights, MomentumDensity(i,j,k,{n})/dm);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::SetObstacleNodes(const PhaseField& Phase, const Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    if (Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        MomentumDensity(i,j,k,{n}).set_to_zero();
        if (Do_TwoPhase)
        {
            DensityWetting (i,j,k,{n}) = 0.0;
            for(auto it = Phase.Fields(i,j,k).cbegin();
                    it != Phase.Fields(i,j,k).cend(); it++)
            if (Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
            {
                size_t pIndex = Phase.FieldsProperties[it->index].Phase;

                DensityWetting (i,j,k,{n}) += it->value * Wetting[n][pIndex];

                if (Phase.FieldsProperties[it->index].Mobile)
                {
                    MomentumDensity(i,j,k,{n}) += Vel.Phase(i,j,k,{pIndex}) *
                        DensityWetting(i,j,k,{n}) * it->value;
                }
            }
        }
        else
        {
            DensityWetting (i,j,k,{n}) = dRho;
        }
        lbPopulations(i,j,k,{n}) =
            EquilibriumDistribution(DensityWetting(i,j,k,{n})/dRho, lbWeights, MomentumDensity(i,j,k,{n})/dm);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
D3Q27 FlowSolverLBM::EquilibriumDistributionTC(dVector3 lbvel, double lbph, double lbnut,  double lbDivVel, 
                                               dVector3 lbFb, const double weights[3][3][3], double lbrho,
                                               dVector3 lbGradRho, size_t n)
{
    D3Q27 locPopulations;
    double lbu2         = lbvel[0]*lbvel[0]+lbvel[1]*lbvel[1]+lbvel[2]*lbvel[2];
    lbtau [n]           = lbnut/lbcs2 + 0.5; 
    double factor       = (1.0-1.0/(2.0*lbtau[n]));

    for(int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
    for(int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
    for(int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
    {
        double CR = (ii-lbvel[0])*lbGradRho[0]+(jj-lbvel[1])*lbGradRho[1]+(kk-lbvel[2])*lbGradRho[2];
        double CF = (ii-lbvel[0])*lbFb[0]+(jj-lbvel[1])*lbFb[1]+(kk-lbvel[2])*lbFb[2];
        double lbw  = lbWeights[ii+1][jj+1][kk+1];
        double lbcu = ii*lbvel[0] + jj*lbvel[1] + kk*lbvel[2];
        double Feq  = lbrho * lbw *(1.0 + lbcu/lbcs2- lbu2/(2.0*lbcs2) + lbcu*lbcu/(2.0*lbcs2*lbcs2) );
        double geq  = Feq*lbcs2 + (lbph - lbrho * lbcs2 ) * lbw; 

	    locPopulations(ii,jj,kk) = geq - 0.5 * (lbrho * lbcs2 * lbw * lbDivVel * factor  + CR*(Feq/lbrho-lbw)*lbcs2*factor + CF*Feq/lbrho);
    }
   return locPopulations;
}

void FlowSolverLBM::SetDivUnearObstZero()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (Obstacle(i,j,k))
    {
        for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
        for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
        for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
        {
            DivVel(i+ii,j+jj,k+kk,{0}) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::SetDivVelBCNXZero()
{

#ifdef MPI_PARALLEL
{
    if(MPI_CART_RANK[0]==MPI_CART_SIZE[0]-1)
    {
        for (int j=-DensityWetting.BcellsY(); j<Grid.Ny+DensityWetting.BcellsY(); j++)
        {
            for (int k=-DensityWetting.BcellsZ(); k<Grid.Nz+DensityWetting.BcellsZ(); k++)
            {
                DivVel(Grid.Nx-1,j,k,{0})=0.0;
            }
        }
    }
}
#endif

}

void FlowSolverLBM::SetDivUnearBCsZeroXZPlane()
{

#ifdef MPI_PARALLEL
{
    if(MPI_CART_RANK[0]==0)
    {
        for (int j=-DensityWetting.BcellsY(); j<Grid.Ny+DensityWetting.BcellsY(); j++)
        {
            for (int k=-DensityWetting.BcellsZ(); k<Grid.Nz+DensityWetting.BcellsZ(); k++)
            {
                DivVel(0,j,k,{0})=0.0;
            }
        }
    }

    if(MPI_CART_RANK[0]==MPI_CART_SIZE[0]-1)
    {
        for (int j=-DensityWetting.BcellsY(); j<Grid.Ny+DensityWetting.BcellsY(); j++)
        {
            for (int k=-DensityWetting.BcellsZ(); k<Grid.Nz+DensityWetting.BcellsZ(); k++)
            {
                DivVel(Grid.Nx-1,j,k,{0})=0.0;
            }
        }
    }

    if(Grid.dNy!=0)
    {
        if(MPI_CART_RANK[1]==0)
        {
            for (int k=-DensityWetting.BcellsZ(); k<Grid.Nz+DensityWetting.BcellsZ(); k++)
            {
                for (int i=-DensityWetting.BcellsX(); i<Grid.Nx+DensityWetting.BcellsX(); i++)
                {
                    DivVel(i,0,k,{0})=0.0;
                }
            }
        }
        if(MPI_CART_RANK[1]==MPI_CART_SIZE[1]-1)
        {
            for (int k=-DensityWetting.BcellsZ(); k<Grid.Nz+DensityWetting.BcellsZ(); k++)
            {
                for (int i=-DensityWetting.BcellsX(); i<Grid.Nx+DensityWetting.BcellsX(); i++)
                {
                    DivVel(i,Grid.Ny-1,k,{0})=0.0;
                }
            }
        }
    }

    if(MPI_CART_RANK[2]==0)
    {
        for (int j=-DensityWetting.BcellsY(); j<Grid.Ny+DensityWetting.BcellsY(); j++)
        {
            for (int i=-DensityWetting.BcellsX(); i<Grid.Nx+DensityWetting.BcellsX(); i++)
            {
                DivVel(i,j,0,{0})=0.0;
            }
        }
    }

    if(MPI_CART_RANK[2]==MPI_CART_SIZE[2]-1)
    {
        for (int j=-DensityWetting.BcellsY(); j<Grid.Ny+DensityWetting.BcellsY(); j++)
        {
            for (int i=-DensityWetting.BcellsX(); i<Grid.Nx+DensityWetting.BcellsX(); i++)
            {
                DivVel(i,j,Grid.Nz-1,{0})=0.0;
            }
        }
    }
}
#endif

}

void FlowSolverLBM::Solve(PhaseField& Phase, Velocities& Vel,
        const BoundaryConditions& BC)
{
    DetectObstacles(Phase);
    SetObstacleNodes(Phase, Vel);
    if(ObstaclesChanged) CalculateDensityAndMomentum();

    Propagation(Phase,BC);
    CalculateDensityAndMomentum();

    BC.SetX(DensityWetting);
    BC.SetY(DensityWetting);
    BC.SetZ(DensityWetting);

    ApplyForces(Phase, Vel);
    Collision();
    CalculateDensityAndMomentum(); //Commented by Dmitry // Uncommented by Raphael
    SetBoundaryConditions(BC);

    CalculateFluidVelocities(Vel, Phase, BC);

    #ifdef MPI_PARALLEL
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].Exist and Phase.FieldsProperties[idx].State == AggregateStates::Solid)
    {
        auto locForce = Phase.FieldsProperties[idx].Force;
        OP_MPI_Allreduce(&locForce[0], &(Phase.FieldsProperties[idx].Force[0]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locForce[1], &(Phase.FieldsProperties[idx].Force[1]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locForce[2], &(Phase.FieldsProperties[idx].Force[2]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);

        auto locTorque = Phase.FieldsProperties[idx].Torque;
        OP_MPI_Allreduce(&locTorque[0], &(Phase.FieldsProperties[idx].Torque[0]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locTorque[1], &(Phase.FieldsProperties[idx].Torque[1]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locTorque[2], &(Phase.FieldsProperties[idx].Torque[2]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    }
    #endif
}

void FlowSolverLBM::Solve(PhaseField& Phase, const Composition& Cx,
        Velocities& Vel, const BoundaryConditions& BC)
{
    DetectObstacles(Phase);
    SetObstacleNodes(Phase, Vel);
    if(ObstaclesChanged) CalculateDensityAndMomentum();

    Propagation(Phase, BC);
    CalculateDensityAndMomentum();

    BC.SetX(DensityWetting);
    BC.SetY(DensityWetting);
    BC.SetZ(DensityWetting);

    ApplyForces(Phase, Vel, Cx);
    Collision();
    CalculateDensityAndMomentum(); //Commented by Dmitry // Uncommented by Raphael
    SetBoundaryConditions(BC);

    CalculateFluidVelocities(Vel, Phase, BC);

    #ifdef MPI_PARALLEL
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].Exist and Phase.FieldsProperties[idx].State == AggregateStates::Solid)
    {
        auto locForce = Phase.FieldsProperties[idx].Force;
        OP_MPI_Allreduce(&locForce[0], &(Phase.FieldsProperties[idx].Force[0]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locForce[1], &(Phase.FieldsProperties[idx].Force[1]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locForce[2], &(Phase.FieldsProperties[idx].Force[2]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);

        auto locTorque = Phase.FieldsProperties[idx].Torque;
        OP_MPI_Allreduce(&locTorque[0], &(Phase.FieldsProperties[idx].Torque[0]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locTorque[1], &(Phase.FieldsProperties[idx].Torque[1]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locTorque[2], &(Phase.FieldsProperties[idx].Torque[2]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    }
    #endif
}

double FlowSolverLBM::CalculateLiquidVolume(size_t Fluid_Comp)
{
    double VolumeLiquid = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DensityWetting, 0, reduction(+:VolumeLiquid))
    if (not Obstacle(i,j,k))
    {
        const double DistLiquid = abs(DensityWetting(i,j,k,{Fluid_Comp}) - LiquidDensity[Fluid_Comp]);
        const double DistVapor  = abs(DensityWetting(i,j,k,{Fluid_Comp}) - VaporDensity [Fluid_Comp]);

        if (DistLiquid < DistVapor) VolumeLiquid += 1.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    auto tmp = VolumeLiquid;
    OP_MPI_Allreduce(&tmp, &(VolumeLiquid), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    VolumeLiquid *= Grid.CellVolume(true);

    return VolumeLiquid;
}

double FlowSolverLBM::CalculateLiquidMass(size_t Fluid_Comp)
{
    double MassLiquid = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DensityWetting, 0, reduction(+:MassLiquid))
    if (not Obstacle(i,j,k))
    {
        const double DistLiquid = abs(DensityWetting(i,j,k,{Fluid_Comp}) - LiquidDensity[Fluid_Comp]);
        const double DistVapor  = abs(DensityWetting(i,j,k,{Fluid_Comp}) - VaporDensity [Fluid_Comp]);

        if (DistLiquid < DistVapor)
        {
            for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
            {
                MassLiquid += lbPopulations(i,j,k,{Fluid_Comp})(ii,jj,kk);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    auto tmp = MassLiquid;
    OP_MPI_Allreduce(&tmp, &(MassLiquid), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    MassLiquid *= dRho*Grid.CellVolume(true);
    return MassLiquid;
}

double FlowSolverLBM::CalculateVaporVolume(size_t Fluid_Comp)
{
    double VolumeVapor = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DensityWetting, 0, reduction(+:VolumeVapor))
    if (not Obstacle(i,j,k))
    {
        const double DistLiquid = abs(DensityWetting(i,j,k,{Fluid_Comp}) - LiquidDensity[Fluid_Comp]);
        const double DistVapor  = abs(DensityWetting(i,j,k,{Fluid_Comp}) - VaporDensity [Fluid_Comp]);

        if (DistLiquid > DistVapor) VolumeVapor += 1.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    auto tmp = VolumeVapor;
    OP_MPI_Allreduce(&tmp, &(VolumeVapor), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    VolumeVapor *= Grid.CellVolume(true);
    return VolumeVapor;
}

double FlowSolverLBM::CalculateVapordMass(size_t Fluid_Comp)
{
    double MassVapor = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DensityWetting, 0, reduction(+:MassVapor))
    if (not Obstacle(i,j,k))
    {
        const double DistLiquid = abs(DensityWetting(i,j,k,{Fluid_Comp}) - LiquidDensity[Fluid_Comp]);
        const double DistVapor  = abs(DensityWetting(i,j,k,{Fluid_Comp}) - VaporDensity [Fluid_Comp]);

        if (DistLiquid > DistVapor)
        {
            for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
            {
                MassVapor += lbPopulations(i,j,k,{Fluid_Comp})(ii,jj,kk);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    auto tmp = MassVapor;
    OP_MPI_Allreduce(&tmp, &(MassVapor), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    MassVapor *= dRho*Grid.CellVolume(true);
    return MassVapor;
}

std::vector<double> FlowSolverLBM::CalculateFluidMass(void) const
{
    std::vector<double> Mass(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::vector<double> tmpMass;
        // Calculate fluid momentum
        #pragma omp parallel
        {
            std::vector<double> tmplocMass;
            #pragma omp for collapse(OMP_COLLAPSE_LOOPS)\
            schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
            for (int i = 0; i < Grid.Nx; i++)
            for (int j = 0; j < Grid.Ny; j++)
            for (int k = 0; k < Grid.Nz; k++)
            if (not Obstacle(i,j,k))
            for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
            {
                double rr = lbPopulations(i,j,k,{n})(ii,jj,kk);
                tmplocMass.push_back(rr);
            }

            std::sort(tmplocMass.begin(),tmplocMass.end(), [] (double a, double b) { return a < b;});

            double locMass = 0.0;
            for (auto value : tmplocMass) locMass += value;

            #pragma omp critical
            {
                tmpMass.push_back(locMass);
            }
        }

        std::sort(tmpMass.begin(),tmpMass.end(), [] (double a, double b) { return a < b;});

        for (auto value : tmpMass) Mass[n] += value;

        #ifdef MPI_PARALLEL
        auto tmp = Mass[n];
        OP_MPI_Allreduce(&tmp, &(Mass[n]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        Mass[n] *= dRho*Grid.CellVolume(true);
    }
    return Mass;
}

void FlowSolverLBM::CalculateDensityAndMomentum(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        DensityWetting (i,j,k,{n}) = 0.0;
        MomentumDensity(i,j,k,{n}).set_to_zero();
        for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
        for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
        for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
        {
            const dVector3 vel {ii*Grid.dx/dt,jj*Grid.dx/dt,kk*Grid.dx/dt};
            DensityWetting (i,j,k,{n}) +=     lbPopulations(i,j,k,{n})(ii,jj,kk)*dRho;
            MomentumDensity(i,j,k,{n}) += vel*lbPopulations(i,j,k,{n})(ii,jj,kk)*dRho;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    FixPopulations();
}

void FlowSolverLBM::CalculateFluidVelocities(Velocities& Vel,
        const PhaseField& Phase, const BoundaryConditions& BC) const
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); it++)
    if (Phase.FieldsProperties[it->index].State != AggregateStates::Solid)
    {
        size_t PhaseIdx = Phase.FieldsProperties[it->index].Phase;
        Vel.Phase(i,j,k,{PhaseIdx}) = FluidVelocity(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Vel.CalculateAverage(Phase);
    Vel.SetBoundaryConditions(BC);
}

double FlowSolverLBM::Pressure(const int i, const int j, const int k) const
{
    double Pressure = 0.0;
    if(!Obstacle(i,j,k))
    {
        if (Do_Kupershtokh)
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k,{n}) > DBL_EPSILON)
            {
                Pressure += lbCriticalPressure[n]*dP*VanDerWaalsGas::ReducedPressure(
                        DensityWetting(i,j,k,{n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho,lbCriticalTemperature[n]);
            }
        }
        else if (Do_Benzi)
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k,{n}) > DBL_EPSILON)
            {
                const double& rho   = DensityWetting(i,j,k,{n});
                const double locPsi = psi(rho,rho_0[n]);

                Pressure += cs2*rho + 0.5*dt*cs2*Gb[n][n]*locPsi*locPsi;
            }
        }
        else if (Do_ThermalComp)
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k,{n}) > DBL_EPSILON)
            {
                Pressure += HydroPressure(i,j,k,{n});
            }
        }
        else
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k,{n}) > DBL_EPSILON)
            {
                const double& rho = DensityWetting(i,j,k,{n});

                Pressure += cs2*rho;
            }
        }
    }
    return Pressure;
}

dMatrix3x3 FlowSolverLBM::PressureTensor(const int i, const int j, const int k) const
{
    dMatrix3x3 Pressure({0,0,0,0,0,0,0,0,0});
    if (Do_Kupershtokh)
    {
        dMatrix3x3 lbPressure({0,0,0,0,0,0,0,0,0});
        //Info::WriteWarning("Pressure tensor calculation is inaccurate", thisclassname, "PressureTensor");
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            double LaplacianPhi    = 0.0;
            double GradPhi[3]      = {0,0,0};
            const double weight[3] = {-0.5,0.0,0.5};
            for (int ii = -1; ii <= 1; ii++)
            {
                const double pi = VanDerWaalsGas::ReducedPressure(DensityWetting(i+ii,j,k,{n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                const double pj = VanDerWaalsGas::ReducedPressure(DensityWetting(i,j+ii,k,{n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                const double pk = VanDerWaalsGas::ReducedPressure(DensityWetting(i,j,k+ii,{n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);

                GradPhi[0] += weight[ii+1]*Phi(DensityWetting(i+ii,j,k,{n})/dRho, pi, GasParameter[n]);
                GradPhi[1] += weight[ii+1]*Phi(DensityWetting(i,j+ii,k,{n})/dRho, pj, GasParameter[n]);
                GradPhi[2] += weight[ii+1]*Phi(DensityWetting(i,j,k+ii,{n})/dRho, pk, GasParameter[n]);
            }

            for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
            {
                const double pijk = VanDerWaalsGas::ReducedPressure(DensityWetting(i+ii,j+jj,k+kk,{n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                LaplacianPhi += LaplacianStencil3D_27a[ii+1][jj+1][kk+1]*
                    Phi(DensityWetting(i+ii,j+jj,k+kk,{n})/dRho, pijk , GasParameter[n]);
            }

            const double locPhi  = Phi(DensityWetting(i,j,k,{n})/dRho, VanDerWaalsGas::ReducedPressure(DensityWetting(i,j,k,{n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]), GasParameter[n]);

            const double GradPhi2 =
                GradPhi[0]*GradPhi[0] +
                GradPhi[1]*GradPhi[1] +
                GradPhi[2]*GradPhi[2];

            const double lbp0 = lbcs2*DensityWetting(i,j,k,{n})/dRho
                - 6.0*Gb[n][n]/dRho*0.50*lbcs2*locPhi*locPhi
                - 6.0*Gb[n][n]/dRho*0.50*lbcs2*lbcs2*locPhi*LaplacianPhi
                - 6.0*Gb[n][n]/dRho*0.25*lbcs2*lbcs2*GradPhi2;

            lbPressure(0,0) += lbp0;
            lbPressure(1,1) += lbp0;
            lbPressure(2,2) += lbp0;

            for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
            {
                lbPressure(ii,jj) += 6.0*Gb[n][n]/dRho*0.5*lbcs2*lbcs2*GradPhi[ii]*GradPhi[jj];
            }
        }
        Pressure = lbPressure*dP;
    }
    else if (Do_Benzi)
    {
        //READ: Shan - 2008 - Pressure tensor calculation in a class of nonideal...
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const double& rho = DensityWetting(i,j,k,{n});
            const double  p0  = cs2*rho;
            Pressure(0,0) += p0;
            Pressure(1,1) += p0;
            Pressure(2,2) += p0;

            for (size_t m = 0; m < N_Fluid_Comp; ++m)
            if (std::abs(Gb[n][m]) > DBL_EPSILON)
            for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
            {
                const dVector3 vel {ii*Grid.dx/dt,jj*Grid.dx/dt,kk*Grid.dx/dt};

                Pressure += vel.dyadic(vel) * 0.5 * Gb[n][m] * psi(rho, rho_0[n]) * lbWeights[ii+1][jj+1][kk+1] * psi(DensityWetting(i+ii,j+jj,k,{m}), rho_0[m]);
            }
        }
    }
    else
    {
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const double p0 = cs2*DensityWetting(i,j,k,{n});
            Pressure(0,0) += p0;
            Pressure(1,1) += p0;
            Pressure(2,2) += p0;
        }
    }
    return Pressure;
}

double  FlowSolverLBM::OptimalParaKuper(const double ReducedTemperature,
        const double GasPrameter)
{
    const std::vector<std::vector<double>> OptimalValues =
        {{2.0000000000000018e-01,-1.4092641421153448e-01},
         {3.0000000000000016e-01,-1.4776128416548723e-01},
         {4.0000000000000013e-01,-1.4973625604844859e-01},
         {5.0000000000000011e-01,-1.5026704234216137e-01},
         {6.0000000000000009e-01,-1.5023912552780752e-01},
         {7.0000000000000007e-01,-1.4788430943930206e-01},
         {8.0000000000000004e-01,-1.4180410847718788e-01},
         {9.0000000000000002e-01,-7.8369694891999728e-02}};

    // Return smallest temperature if temperature is too small
    if (OptimalValues.front()[0] > ReducedTemperature) return OptimalValues.front()[1];

    // Search for optimal value
    std::vector<double> PreviousValue = OptimalValues[0];
    std::vector<double> Interpolated  = {0.0,0.0};
    for (auto& Value: OptimalValues)
    {
        if (Value[0] == ReducedTemperature)
        {
            return Value[1];
        }
        else if (ReducedTemperature < Value[0] and
                 ReducedTemperature > PreviousValue[0])
        {
            Interpolated = PreviousValue;

            const double dT = (Value[0]-PreviousValue[0]);
            const double DT = (ReducedTemperature-PreviousValue[0]);
            const double A  = DT/dT;

            Interpolated[0] = ReducedTemperature;
            Interpolated[1] += (Value[1] - PreviousValue[1])*A;

            return Interpolated[1];
        }
        PreviousValue = Value;
    }
    return PreviousValue[1];
}

D3Q27 FlowSolverLBM::EquilibriumDistribution(double lbDensity,
        const double weights[3][3][3], dVector3 lbMomentum)
{
   dVector3 macroVel;
   if(lbDensity != 0)
   {
       macroVel = lbMomentum/lbDensity;
   }
   else
   {
       macroVel.set_to_zero();
   }
   double u2 = (macroVel*macroVel);
   D3Q27 locPopulations;
   for(int x = -1; x <= 1; x++)
   for(int y = -1; y <= 1; y++)
   for(int z = -1; z <= 1; z++)
   {
       double cu = x*macroVel[0] + y*macroVel[1] + z*macroVel[2];
       locPopulations(x,y,z) = lbDensity*weights[x+1][y+1][z+1]*(1.0 - 1.5*u2 + cu*(3.0 + 4.5*cu));
   }
   return locPopulations;
};

double FlowSolverLBM::DensityProfile(const double x, const size_t n) const
{
    if (Do_TwoPhase)
    {
        return 0.5 * (LiquidDensity[n]+VaporDensity[n]) -
               0.5 * (LiquidDensity[n]-VaporDensity[n]) *
               std::tanh(x/InterfaceWidth[n]);
    }
    else
    {
        return 1.0;
    }
}

double FlowSolverLBM::FluidDensity(const int i, const int j, const int k) const
{
    double locFluidDensity = 0.0;
    for (size_t comp = 0; comp < N_Fluid_Comp; ++comp)
    {
        locFluidDensity += FluidDensity(i,j,k,comp);
    }
    return locFluidDensity;
}

double FlowSolverLBM::Density(const PhaseField& Phase,
        const int i, const int j, const int k) const
{

    double SolidDensity = 0.0;
    for (auto it = Phase.Fields(i,j,k).cbegin();
             it != Phase.Fields(i,j,k).cend(); it++)
    if (Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
    {
        SolidDensity += it->value*Phase.FieldsProperties[it->index].Density;
    }
    double locSolidFraction = SolidFraction(i,j,k,Phase);
    //double Density = FluidDensity(i,j,k)*(1.0-locSolidFraction) + SolidDensity*locSolidFraction; //TODO this should be uses
    double Density = DensityWetting(i,j,k,{0})*(1.0-locSolidFraction) + SolidDensity*locSolidFraction;

    return Density;
}

dVector3 FlowSolverLBM::FluidVelocity(const int i, const int j, const int k,
        const size_t comp) const
{
    if ((not Obstacle(i,j,k)) and (DensityWetting(i,j,k,{comp}) > DBL_EPSILON))
    {
        return (MomentumDensity(i,j,k,{comp}))/DensityWetting(i,j,k,{comp});
    }
    else return dVector3{0.0,0.0,0.0};
}

dVector3 FlowSolverLBM::FluidVelocity(const int i, const int j, const int k) const
{
    double locFluidDensity = FluidDensity(i,j,k);
    if ((not Obstacle(i,j,k)) and locFluidDensity > DBL_EPSILON)
    {
        dVector3 locMomentumDensity{0.0,0.0,0.0};
        for (size_t comp = 0; comp < N_Fluid_Comp; ++comp)
            locMomentumDensity += MomentumDensity(i,j,k,{comp});

        return locMomentumDensity/locFluidDensity;
    }
    else return dVector3{0.0,0.0,0.0};
}

void FlowSolverLBM::FinalizeInitialiation(const PhaseField& Phase,
        const Velocities& Vel, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t comp = 0; comp < N_Fluid_Comp; ++comp)
    {
        double   lbDensity  = DensityWetting (i,j,k,{comp})/dRho;
        dVector3 lbMomentum = MomentumDensity(i,j,k,{comp})/dm;
        lbPopulations(i,j,k,{comp}) = FlowSolverLBM::EquilibriumDistribution(lbDensity, lbWeights, lbMomentum);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    DetectObstacles(Phase);
    SetObstacleNodes(Phase,Vel);
    FixPopulations();
    EnforceMassConservation();
    SetBoundaryConditions(BC);
}
}
