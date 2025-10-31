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
#include "BoundaryConditions.h"
#include "GrandPotential/GrandPotentialDensity.h"
#include "GrandPotential/GrandPotentialSolver.h"
#include "Includes.h"
#include "InterfaceProperties.h"
#include "NumericalMethods/RootFindingAlgorithms.h"
#include "NumericalMethods/SystemOfLinearEquationsSolvers.h"
#include "PhaseField.h"
#include "PhysicalConstants.h"
#include "Settings.h"
#include "Temperature.h"
#include "Velocities.h"
#include "VTK.h"

namespace openphase
{
GrandPotentialSolver::GrandPotentialSolver(Settings& locSettings)
{
    Initialize(locSettings);
    ReadInput(DefaultInputFileName);
}
GrandPotentialSolver::GrandPotentialSolver(Settings& locSettings, std::string filename)
{
    Initialize(locSettings);
    ReadInput(filename);
}
void GrandPotentialSolver::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname  = "GrandPotentialSolver";
    thisobjectname = thisclassname + ObjectNameSuffix;

    advectable = true;
    readable   = true;

    Grid = locSettings.Grid;

    ElementNames = locSettings.ElementNames;
    Ncomp        = locSettings.Ncomp;
    Nphases      = locSettings.Nphases;
    PhaseNames   = locSettings.PhaseNames;
    RawDataDir   = locSettings.RawDataDir;

    CInitial               .Allocate({Nphases, Ncomp});
    PhaseMobilities        .Allocate({Nphases, Ncomp});
    PhaseActivationEnergies.Allocate({Nphases, Ncomp});
    InterfaceMobilities    .Allocate({Nphases, Nphases, Ncomp});

    dPhaseMobilities_dConcentration.Allocate({Nphases, Ncomp, Ncomp});

    InitialChemicalPotential.assign(Ncomp, 0.0);
    MolarMasses             .assign(Ncomp, 0.0);
    TOC0                    .assign(Ncomp, 0.0);

    Bcells = Grid.Bcells;
    if (Bcells < 2)
    {
        ConsoleOutput::WriteExit("Too few boundary cell! Two boundary cells are required!", thisclassname, "ReadInput");
        std::exit(EXIT_FAILURE);
    };

    Concentrations      .Allocate(Grid, {Ncomp}, Bcells);
    ConcentrationsDot   .Allocate(Grid, {Ncomp}, Bcells);
    ChemicalPotential   .Allocate(Grid, {Ncomp}, Bcells);
    ChemicalPotentialDot.Allocate(Grid, {Ncomp}, Bcells);
    DiffusionFlux       .Allocate(Grid, {Ncomp}, Bcells);

    initialized = false;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}
void GrandPotentialSolver::ReadInput(const std::string filename)
{
    std::fstream inp(filename.c_str(), std::ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + filename + "\" could not be opened", thisclassname, "ReadInput");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert(thisclassname);
    ConsoleOutput::WriteStandard("Source", filename);

    ReadInput(inp_data);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}
void GrandPotentialSolver::ReadInput(std::stringstream& inp_data)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp_data, thisclassname);

    // Read component masses
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+ElementNames[comp];
        MolarMasses              [comp] = FileInterface::ReadParameterD(inp_data, moduleLocation, "MASS_"+converter, false, 0.0);
        InitialChemicalPotential [comp] = FileInterface::ReadParameterD(inp_data, moduleLocation, "ICP_"+converter, false, 0.0);
    }

    // Read initial concentrations
    UseInitialPressure = FileInterface::ReadParameterB(inp_data, moduleLocation, "UIP",false,false);
    if (UseInitialPressure)
    {
        InitialPressure = FileInterface::ReadParameterD(inp_data, moduleLocation, "IP");
    }
    else
    {
        for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            std::string converter = ""+std::to_string(PhaseIdx)+"_"+ElementNames[comp];
            CInitial ({PhaseIdx, comp}) = FileInterface::ReadParameterD(inp_data, moduleLocation, "CI_"+converter);
        }
    }

    // Read bulk mobilities
    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t comp     = 0; comp     < Ncomp;       comp++)
    {
        std::string converter = ""+std::to_string(PhaseIdx)+"_"+ElementNames[comp];
        PhaseMobilities({PhaseIdx, comp}) = FileInterface::ReadParameterD(inp_data, moduleLocation, "M0_"+converter);
    }

    // Read bulk activation energies
    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t comp     = 0; comp     < Ncomp;       comp++)
    {
        std::string converter = ""+std::to_string(PhaseIdx)+"_"+ElementNames[comp];
        PhaseActivationEnergies({PhaseIdx, comp}) = FileInterface::ReadParameterD(inp_data, moduleLocation, "EA_"+converter, false, 0.0);
    }

    // Read bulk mobilities concentration dependency
    bool PhaseMobilityConcentrationCoupling = FileInterface::ReadParameterB(inp_data, moduleLocation, "dMdc");
    if (PhaseMobilityConcentrationCoupling)
    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t compA    = 0; compA     < Ncomp;     compA++)
    for(size_t compB    = 0; compB     < Ncomp;     compB++)
    {
        std::string converter = ""+std::to_string(PhaseIdx)+"_"+ElementNames[compA]+"_"+ElementNames[compB];
        double valueAB = FileInterface::ReadParameterD(inp_data, moduleLocation, "dMdc_"+converter);
        dPhaseMobilities_dConcentration({PhaseIdx, compA, compB}) = valueAB;
        dPhaseMobilities_dConcentration({PhaseIdx, compB, compA}) = valueAB;
    }

    // Read interface mobilities
    bool InterfaceMobilityConstants = FileInterface::ReadParameterB(inp_data, moduleLocation, "IM");
    if (InterfaceMobilityConstants)
    for(size_t alpha = 0;     alpha < Nphases; alpha++)
    for(size_t beta  = alpha; beta  < Nphases;  beta++)
    for(size_t comp  = 0;     comp  < Ncomp;    comp++)
    {
        std::string converterAB = std::to_string(alpha) + "_" + std::to_string(beta) + "_" + ElementNames[comp];
        //std::string converterBA = std::to_string(beta ) + "_" + std::to_strinf(alpha) + "_" + ElementNames[comp];
        double valueAB = FileInterface::ReadParameterD(inp_data, moduleLocation, "IM_" + converterAB);
        InterfaceMobilities({alpha, beta, comp}) = valueAB;
        InterfaceMobilities({beta, alpha, comp}) = valueAB;
    }

    UseImplicitSolver = FileInterface::ReadParameterB(inp_data, moduleLocation, "Implicit");
    if (UseImplicitSolver)
    {
        ChemicalPotentialAccuracy = FileInterface::ReadParameterD(inp_data, moduleLocation, "ACC");
        MaxIterations             = FileInterface::ReadParameterI(inp_data, moduleLocation, "MAXI");

        ChemicalPotentialOld .Allocate(Grid, {Ncomp}, Bcells);
        ChemicalPotentialDot2.Allocate(Grid, {Ncomp}, Bcells);
    }
    else
    {
        ChemicalPotentialAccuracy = 0.0;
        MaxIterations             = 0;
    }

    ConserveTOC = FileInterface::ReadParameterB(inp_data, moduleLocation, "TOC");
    if (ConserveTOC)
    {
        TOCAccuracy      = FileInterface::ReadParameterD(inp_data, moduleLocation, "TOCACC");
        TOCMaxIterations = FileInterface::ReadParameterI(inp_data, moduleLocation, "TOCMAXI");
    }
    else
    {
        TOCAccuracy      = 0.0;
        TOCMaxIterations = 0;
    }
}
void GrandPotentialSolver::Advect(AdvectionHR& Adv, const Velocities& Vel,
        PhaseField& Phi, const BoundaryConditions& BC, const double dt,
        const double tStep)
{
    //NOTE Advection chemical potential does not work probably
    //if(ChemicalPotentialDot.IsNotAllocated()) ChemicalPotentialDot.Allocate(ChemicalPotential);
    //if(ChemicalPotentialOld.IsNotAllocated()) ChemicalPotentialOld.Allocate(ChemicalPotential);
    //Adv.AdvectField(ChemicalPotential,ChemicalPotentialDot,ChemicalPotentialOld,Vel,BC,dx,dt,tStep);

    Adv.CalculateAdvectionIncrements(Concentrations,ConcentrationsDot,Vel,BC,dt);
}
void GrandPotentialSolver::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(ChemicalPotential);
    BC.SetY(ChemicalPotential);
    BC.SetZ(ChemicalPotential);
    BC.SetX(Concentrations);
    BC.SetY(Concentrations);
    BC.SetZ(Concentrations);
}
void GrandPotentialSolver::Write(const char* filename) const
{
    std::ofstream out;
#ifdef MPI_PARALLEL
    for (int proc = 0; proc < MPI_SIZE; proc++)
    {
        if (MPI_RANK == proc)
        {
            if (MPI_RANK == 0)
            {
                out.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
            }
            else
            {
                out.open(filename, std::ios::out | std::ios::binary | std::ios::app);
            }
            if (MPI_RANK == 0)
            {
                out.write(reinterpret_cast<const char*>(&Grid.TotalNx), sizeof(long int));
                out.write(reinterpret_cast<const char*>(&Grid.Ny     ), sizeof(long int));
                out.write(reinterpret_cast<const char*>(&Grid.Nz     ), sizeof(long int));
                out.write(reinterpret_cast<const char*>(&Ncomp  ), sizeof(size_t));
                for(size_t comp = 0; comp < Ncomp; comp++)
                {
                    out.write(reinterpret_cast<const char*>(&TOC0[comp]), sizeof(size_t));
                }
            }
            for (long int i = 0; i < Grid.Nx;    i++)
            for (long int j = 0; j < Grid.Ny;    j++)
            for (long int k = 0; k < Grid.Nz;    k++)
            for (size_t   n = 0; n < Ncomp; n++)
            {
               out.write(reinterpret_cast<const char*>(&ChemicalPotential(i,j,k)({n})), sizeof(double));
            }
            out.close();
        }
        OP_MPI_Barrier(OP_MPI_COMM_WORLD);
    }
#else
    out.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    out.write(reinterpret_cast<const char*>(&Grid.TotalNx), sizeof(long int));
    out.write(reinterpret_cast<const char*>(&Grid.Ny     ), sizeof(long int));
    out.write(reinterpret_cast<const char*>(&Grid.Nz     ), sizeof(long int));
    out.write(reinterpret_cast<const char*>(&Ncomp  ), sizeof(size_t));
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        out.write(reinterpret_cast<const char*>(&TOC0[comp]), sizeof(size_t));
    }
    for (long int i = 0; i < Grid.Nx;    i++)
    for (long int j = 0; j < Grid.Ny;    j++)
    for (long int k = 0; k < Grid.Nz;    k++)
    for (size_t   n = 0; n < Ncomp; n++)
    {
       out.write(reinterpret_cast<const char*>(&ChemicalPotential(i,j,k,{n})), sizeof(double));
    }
    out.close();
#endif
}
void GrandPotentialSolver::Write(const std::string& filename) const
{
    Write(filename.c_str());
}
bool GrandPotentialSolver::Write(const Settings& locSettings, const int tStep) const
{
    std::string filename = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname+'_', tStep, ".dat");
    Write(filename);
    return true;
}
void GrandPotentialSolver::Read(const BoundaryConditions& BC, const PhaseField& Phase, const GrandPotentialDensity& omega, const char* filename)
// TODO write more to file for restart void GrandPotentialSolver::Read(const BoundaryConditions& BC, const char* filename)
{
    std::ifstream inp;
#ifdef MPI_PARALLEL
    long pos = 0;
    for (int proc = 0; proc < MPI_SIZE; proc++)
    {
        if (MPI_RANK == proc)
        {
            inp.open(filename, std::ios::in | std::ios::binary);
            if (MPI_RANK == 0)
            {
                long int locTotalNx = 0;
                long int locNy      = 0;
                long int locNz      = 0;
                size_t   locNcomp   = 0;
                inp.read(reinterpret_cast<char*>(&locTotalNx), sizeof(long int));
                inp.read(reinterpret_cast<char*>(&locNy     ), sizeof(long int));
                inp.read(reinterpret_cast<char*>(&locNz     ), sizeof(long int));
                inp.read(reinterpret_cast<char*>(&locNcomp  ), sizeof(size_t));
                if (locTotalNx != Grid.TotalNx or locNy != Grid.Ny or locNz != Grid.Nz or locNcomp != Ncomp)
                {
                    std::stringstream message;
                    message << thisclassname << " Inconsistent system dimensions!\n"
                            << "Filename: " << filename << "\n"
                            << "Input data dimensions:    (" << locTotalNx << "," << locNy << "," << locNz << ","  << locNcomp << ").\n"
                            << "Required data dimensions: (" <<    Grid.TotalNx << "," <<    Grid.Ny << "," <<    Grid.Nz << ","  <<    Ncomp << ").\n";
                   throw std::runtime_error(message.str());
                }
                for(size_t comp = 0; comp < Ncomp; comp++)
                {
                    inp.read(reinterpret_cast<char*>(&TOC0[comp]), sizeof(size_t));
                    assert(TOC0[comp] != 0);
                }
            }
            else inp.seekg(pos);

            for (long int i = 0; i < Grid.Nx;    i++)
            for (long int j = 0; j < Grid.Ny;    j++)
            for (long int k = 0; k < Grid.Nz;    k++)
            for (size_t   n = 0; n < Ncomp; n++)
            {
                inp.read(reinterpret_cast<char*>(&ChemicalPotential(i,j,k)({n})), sizeof(double));
            }
            pos = inp.tellg();
            inp.close();
        }
        OP_MPI_Barrier(OP_MPI_COMM_WORLD);
        OP_MPI_Bcast(&pos,1,OP_MPI_LONG,proc,OP_MPI_COMM_WORLD);
        if (proc == 0)
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                OP_MPI_Bcast(&(TOC0[comp]),1,OP_MPI_DOUBLE,0,OP_MPI_COMM_WORLD);
                assert(TOC0[comp] != 0);
            }
        }
    }
#else
    inp.open(filename, std::ios::in | std::ios::binary);
    long int locTotalNx = 0;
    long int locNy      = 0;
    long int locNz      = 0;
    size_t   locNcomp   = 0;
    inp.read(reinterpret_cast<char*>(&locTotalNx), sizeof(long int));
    inp.read(reinterpret_cast<char*>(&locNy     ), sizeof(long int));
    inp.read(reinterpret_cast<char*>(&locNz     ), sizeof(long int));
    inp.read(reinterpret_cast<char*>(&locNcomp  ), sizeof(size_t));
    if (locTotalNx != Grid.TotalNx or locNy != Grid.Ny or locNz != Grid.Nz or locNcomp != Ncomp)
    {
        std::stringstream message;
        message << thisclassname << " Inconsistent system dimensions!\n"
                << "Filename: " << filename << "\n"
                << "Input data dimensions:    (" << locTotalNx << "," << locNy << "," << locNz << ","  << locNcomp << ").\n"
                << "Required data dimensions: (" <<    Grid.TotalNx << "," <<    Grid.Ny << "," <<    Grid.Nz << ","  <<    Ncomp << ").\n";
       throw std::runtime_error(message.str());
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        inp.read(reinterpret_cast<char*>(&TOC0[comp]), sizeof(size_t));
    }
    for (long int i = 0; i < Grid.Nx;    i++)
    for (long int j = 0; j < Grid.Ny;    j++)
    for (long int k = 0; k < Grid.Nz;    k++)
    for (size_t   n = 0; n < Ncomp; n++)
    {
        inp.read(reinterpret_cast<char*>(&ChemicalPotential(i,j,k,{n})), sizeof(double));
    }
    inp.close();
#endif

    SetBoundaryConditions(BC);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-1,)
    {
        CalculateLocalConcentrations (i,j,k,Phase,omega);
        //CaclulateLocalDiffusionFlux  (i,j,k,{0.0,0.0,0.0},Phase,Temp);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
}
void GrandPotentialSolver::Read(const BoundaryConditions& BC, const PhaseField& Phase, const GrandPotentialDensity& omega, const std::string& filename)
{
    Read(BC,Phase,omega,filename.c_str());
}
void GrandPotentialSolver::Read(const BoundaryConditions& BC, const PhaseField& Phase, const GrandPotentialDensity& omega, const  long tStep)
{
    std::string filename = FileInterface::MakeFileName(RawDataDir, thisclassname+'_', tStep, ".dat");
    Read(BC,Phase,omega,filename);
}
double GrandPotentialSolver::GrandPotential(const PhaseField& Phase, const GrandPotentialDensity& omega) const
{
    return CalculateVolumeIntegral([&Phase, &omega](long i, long j, long k){return omega(i,j,k,Phase);},Grid.dx);
}
void GrandPotentialSolver::SetInitialConcentration(const PhaseField& Phase, const GrandPotentialDensity& omega)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            // Calculate initial concentration
            double FinalConcentration = 0;
            for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
                FinalConcentration += alpha->value*CInitial({PhaseIdx,comp});
            }

            ChemicalPotential(i,j,k,{comp}) = 0.0;

            // Because an initial chemical potential is needed as input
            // a chemical potential which maps into the desired initial
            // concentration is calculated using newtons method
            double precision = DBL_EPSILON;
            double f,df, dmu = 0;
            int n = 0;
            do
            {
                CalculateLocalConcentrations(i,j,k,Phase,omega);
                f  = Concentrations(i,j,k,{comp}) - FinalConcentration;
                df = Susceptibility(i,j,k,comp,Phase,omega);
                if (df == 0.0) break;
                dmu = -f/df/2;
                ChemicalPotential(i,j,k,{comp})+=dmu;
                n++;
            }
            while (std::abs(dmu) != 0.0 and std::abs(dmu) > std::abs(ChemicalPotential(i,j,k,{comp}))*precision and n < 1000);
        }
        CalculateLocalConcentrations (i,j,k,Phase,omega);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void GrandPotentialSolver::SetInitialPressure(const PhaseField& Phase, const GrandPotentialDensity& omega)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,)
    {
        // Calculate initial concentration
        ChemicalPotential(i,j,k,{0}) = InitialChemicalPotential[0];

        // Because an initial chemical potential is needed as input
        // a chemical potential which maps into the desired initial
        // concentration is calculated using newtons method
        double precision = 0.1;
        auto PressureDelta = [this, &Phase, &omega, i, j, k](double mu)
        {
            //TODO mu should not be unused!!!!
            return (-omega(i,j,k,Phase)) - InitialPressure;
        };
        auto dPressureDelta_dChemicalPotential = [ &Phase, &omega, i, j, k]([[maybe_unused]] double mu)
        {
            double tmp = 0;
            for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
                tmp -= alpha->value*omega(PhaseIdx).dChemicalPotential(i,j,k,0);
            }
            return tmp;
        };

        try
        {
            RootFindingAlgorithms::Newton(PressureDelta, dPressureDelta_dChemicalPotential, ChemicalPotential(i,j,k,{0}), precision, 1000);
        }
        catch (std::runtime_error& ecep)
        {
            ConsoleOutput::WriteWarning(ecep.what(),thisclassname,"SetInitialPressure");
        }
        if (std::abs(-omega(i,j,k,Phase)-InitialPressure) > 1.0 )
        {
            std::cerr << i << " " << j << " " << k << " " << (-omega(i,j,k,Phase)) <<  "\n";
            OP_Exit(EXIT_FAILURE);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    ConsoleOutput::Write("Set Pressure sucessufully!");
}
void GrandPotentialSolver::CalculateCenterOfMassVelocitiesOfGrains(const PhaseField& Phase)
{
    size_t Nthreads = 1;
    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif
    CenterOfMassVelocityGrain.assign(Phase.FieldsProperties.size(),{0.,0.,0.});
    std::vector<std::vector<dVector3>> CenterOfMassVelocityGrainThread(Nthreads,std::vector<dVector3>(Phase.FieldsProperties.size(),{0.,0.,0.}));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,)
    {
        dVector3 locDensityFlux = {0.,0.,0.};
        double locDensity = 0.;
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            locDensityFlux += DiffusionFlux(i,j,k,{comp});
            locDensity += Concentrations(i,j,k,{comp});
        }

        size_t idx = omp_get_thread_num();
        for (auto alpha = Phase.Fields(i,j,k).cbegin();
                  alpha != Phase.Fields(i,j,k).cend(); alpha++)
        if(locDensity > 0.)
        if(Phase.FieldsProperties[alpha->index].is_solid())
        {
            CenterOfMassVelocityGrainThread[idx][alpha->index] += locDensityFlux/locDensity*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    // OpenMP reduction
    for(size_t idx = 0; idx < Nthreads; idx++)
    for(size_t alpha = 0; alpha < CenterOfMassVelocityGrain.size(); alpha++)
    {
        CenterOfMassVelocityGrain[alpha] += CenterOfMassVelocityGrainThread[idx][alpha];
    }

    for(size_t alpha = 0; alpha < CenterOfMassVelocityGrain.size(); alpha++)
    if(Phase.FieldsProperties[alpha].Volume > 0.)
    {
        CenterOfMassVelocityGrain[alpha] /= Phase.FieldsProperties[alpha].Volume;
    }

    // TODO MPI reduction
}
double GrandPotentialSolver::TotalAmountOfComponent(size_t comp) const
{
    return CalculateVolumeIntegral([this,comp](long i, long j, long k){return Concentrations(i,j,k,{comp});},Grid.dx);
}
double GrandPotentialSolver::TotalAmountOfComponentChange(size_t comp) const
{
    return (TotalAmountOfComponent(comp) - TOC0[comp])/TOC0[comp];
}
void GrandPotentialSolver::SetInitial(const PhaseField& Phase, const GrandPotentialDensity& omega, const BoundaryConditions& BC, double Temp)
{
    //TODO use either (mu0_0, mu0_1,..) or (rho0, c_0, c_1, ...) to determine initial state
    if   (UseInitialPressure) SetInitialPressure      (Phase,omega);
    else                      SetInitialConcentration (Phase,omega);

    SetBoundaryConditions(BC);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-1,)
    {
        CalculateLocalConcentrations(i,j,k,Phase,omega);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        TOC0[comp] = TotalAmountOfComponent(comp);
        if (std::abs(TOC0[comp]) < DBL_EPSILON)
        {
            std::stringstream message;
            message << "There is no  amount of " << ElementNames[comp] << " present!\n"
                    << "This may be unintentional and hence an error!";
            ConsoleOutput::WriteWarning(message.str(), thisclassname, "SetInitial");
        }
    }
}
void GrandPotentialSolver::EnforceConservationOfTOC(const PhaseField& Phase, const GrandPotentialDensity& omega)
{
    for (std::size_t comp = 0; comp < Ncomp; comp++)
    {
        auto residual = [&] (double& delta)
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
            {
                ChemicalPotential(i,j,k,{comp}) += delta;
                CalculateLocalConcentrations(i,j,k,comp,Phase,omega);
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            double r = (TotalAmountOfComponent(comp) - TOC0[comp])/TOC0[comp];

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,comp)
            {
                ChemicalPotential(i,j,k,{comp}) -= delta;
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            return r;
        };

        double delta0 = 0;
        double delta1 = 0.1;
        try
        {
            RootFindingAlgorithms::Secant(residual, delta0, delta1, TOCAccuracy, TOCMaxIterations);
        }
        catch (std::runtime_error& ecep)
        {
            ConsoleOutput::WriteWarning(ecep.what(),thisclassname,"EnforceConservation");
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
        {
            ChemicalPotential(i,j,k,{comp}) += delta0;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}
void GrandPotentialSolver::CaclulateLocalDiffusionFlux(long i, long j, long k, dVector3 locVelocity, const PhaseField& Phase, const Temperature& Temp)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        // Calculate chemical potential gradient
        DiffusionFlux(i,j,k,{comp}) = {0.0,0.0,0.0};
        dVector3 locForce = {0.0,0.0,0.0};
        if (Grid.dNx) locForce[0] -= (ChemicalPotential(i+1,j  ,k   ,{comp}) - ChemicalPotential(i-1,j  ,k  ,{comp}))/Grid.dx/2.0;
        if (Grid.dNy) locForce[1] -= (ChemicalPotential(i  ,j+1,k   ,{comp}) - ChemicalPotential(i  ,j-1,k  ,{comp}))/Grid.dx/2.0;
        if (Grid.dNz) locForce[2] -= (ChemicalPotential(i  ,j  ,k +1,{comp}) - ChemicalPotential(i  ,j  ,k-1,{comp}))/Grid.dx/2.0;

        DiffusionFlux(i,j,k,{comp}) += locForce*Mobility(i,j,k,comp,Phase,Temp);
        DiffusionFlux(i,j,k,{comp}) -= locVelocity*Concentrations(i,j,k,{comp});
    }
}
dVector3 GrandPotentialSolver::CenterOfMassVelocity(long i, long j, long k, const PhaseField& Phase) const
{
    dVector3 locCenterOfMassVelocity = {0.,0.,0.};
    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        locCenterOfMassVelocity = (CenterOfMassVelocityGrain[alpha->index] + Phase.FieldsProperties[alpha->index].Vcm)*alpha->value;
    }
    return locCenterOfMassVelocity;
}
double GrandPotentialSolver::MassDensity(long i, long j, long k, const PhaseField& Phase, const GrandPotentialDensity& omega) const
{
    double locMassGrandPotentialDensity = 0;
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        locMassGrandPotentialDensity += Concentration(i,j,k,comp,Phase,omega)*MolarMasses[comp];
    }
    return locMassGrandPotentialDensity;
}
double GrandPotentialSolver::Susceptibility(long i, long j, long k, size_t comp, const PhaseField& Phase, const GrandPotentialDensity& omega) const
{
    // Calculate calculate local bulk susceptibility coefficient
    double locSusceptibility = 0.0;
    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
        double locPhaseSusceptibility = -omega(PhaseIdx).dChemicalPotential2(i,j,k,comp);
        locSusceptibility += alpha->value*locPhaseSusceptibility;
    }
    return locSusceptibility;
}
double GrandPotentialSolver::Concentration(long i, long j, long k, size_t comp, const PhaseField& Phase, const GrandPotentialDensity& omega) const
{
    double locConcentration = 0;
    for (auto alpha = Phase.Fields(i,j,k).cbegin();
              alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
        locConcentration -= alpha->value*omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
    }
    return locConcentration;
}
void GrandPotentialSolver::CalculateLocalConcentrations(long i, long j, long k, size_t comp, const PhaseField& Phase, const GrandPotentialDensity& omega)
{
    Concentrations(i,j,k,{comp}) = 0.0;
    for (auto alpha = Phase.Fields(i,j,k).cbegin();
              alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
        Concentrations(i,j,k,{comp}) -= alpha->value*omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
    }
}
void GrandPotentialSolver::CalculateLocalConcentrations(long i, long j, long k, const PhaseField& Phase, const GrandPotentialDensity& omega)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        CalculateLocalConcentrations(i,j,k,comp,Phase,omega);
    }
}
void GrandPotentialSolver::CalculateLocalIncrements(long i, long j, long k, const PhaseField& Phase, const GrandPotentialDensity& omega, double dt)
{
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        // Calculate change of concentration due to diffusion
        double ConcentrationDot = ConcentrationsDot(i,j,k,{comp});
        if (Grid.dNx) ConcentrationDot -= (DiffusionFlux(i+1,j  ,k   ,{comp})[0]-DiffusionFlux(i-1,j  ,k   ,{comp})[0])/Grid.dx/2.0;
        if (Grid.dNy) ConcentrationDot -= (DiffusionFlux(i  ,j+1,k   ,{comp})[1]-DiffusionFlux(i  ,j-1,k   ,{comp})[1])/Grid.dx/2.0;
        if (Grid.dNz) ConcentrationDot -= (DiffusionFlux(i  ,j  ,k +1,{comp})[2]-DiffusionFlux(i  ,j  ,k -1,{comp})[2])/Grid.dx/2.0;

        double locInverseSusceptibility = 1/Susceptibility(i,j,k,comp,Phase,omega);
        ChemicalPotentialDot(i,j,k,{comp}) += ConcentrationDot*locInverseSusceptibility;

        // Calculate change of concentration due to phase-transformation
        NodeA locPhaseDot = Phase.Dot(i,j,k, dt);
        for (auto alpha = Phase.Fields(i,j,k).cbegin();
                  alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
            double locPhaseConcentration = -omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
            ChemicalPotentialDot(i,j,k,{comp}) -= locPhaseConcentration*locPhaseDot.get_value(alpha->index)*locInverseSusceptibility;
        }
        ConcentrationsDot(i,j,k,{comp}) = 0.0;
    }
}
void GrandPotentialSolver::CalculateLocalIncrements1(long i, long j, long k, const PhaseField& Phase, const GrandPotentialDensity& omega, double dt)
{
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        // Calculate change of concentration due to diffusion
        double ConcentrationDot = 0.0;
        if (Grid.dNx) ConcentrationDot -= (DiffusionFlux(i+1,j  ,k   ,{comp})[0]-DiffusionFlux(i-1,j  ,k   ,{comp})[0])/Grid.dx/2.0;
        if (Grid.dNy) ConcentrationDot -= (DiffusionFlux(i  ,j+1,k   ,{comp})[1]-DiffusionFlux(i  ,j-1,k   ,{comp})[1])/Grid.dx/2.0;
        if (Grid.dNz) ConcentrationDot -= (DiffusionFlux(i  ,j  ,k +1,{comp})[2]-DiffusionFlux(i  ,j  ,k -1,{comp})[2])/Grid.dx/2.0;

        double locInverseSusceptibility = 1/Susceptibility(i,j,k,comp,Phase,omega);
        ChemicalPotentialDot(i,j,k,{comp}) += ConcentrationDot*locInverseSusceptibility;

        // Calculate change of concentration due to phase-transformation due to curvature
        for (auto alpha = Phase.Fields(i,j,k).cbegin();
                  alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
            double locPhaseConcentration = -omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
            double locPhaseDot1 = Phase.Dot1(i,j,k, dt).get_value(alpha->index);
            ChemicalPotentialDot(i,j,k,{comp}) -= locPhaseConcentration*locPhaseDot1*locInverseSusceptibility;
        }
    }
}
void GrandPotentialSolver::CalculateLocalIncrements2(long i, long j, long k, const PhaseField& Phase, const GrandPotentialDensity& omega, double dt)
{
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        double locInverseSusceptibility = 1/Susceptibility(i,j,k,comp,Phase,omega);
        // Calculate change of concentration due to phase-transformation due to pressure difference
        for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
            double locPhaseConcentration = -omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
            double locPhaseDot2 = Phase.Dot2(i,j,k,dt).get_value(alpha->index);
            ChemicalPotentialDot2(i,j,k,{comp}) -= locPhaseConcentration*locPhaseDot2*locInverseSusceptibility;
        }
    }
}
void GrandPotentialSolver::CalculateLocalPhaseFieldIncrements (long i, long j, long k, PhaseField& Phase, const GrandPotentialDensity& omega, const InterfaceProperties& IP)
{
    for (auto it = Phase.FieldsDot(i,j,k).begin(); it != Phase.FieldsDot(i,j,k).end(); ++it)
    {
        Phase.FieldsDot(i,j,k).set_sym2(it->indexA,it->indexB,0.0);
    }

    double Prefactor = Pi*Pi/(8.0*Phase.Grid.Eta*Phase.LocalN(Phase.Fields(i, j, k)));

    for(auto alpha = Phase.Fields(i, j, k).cbegin();
             alpha != Phase.Fields(i, j, k).cend(); ++alpha)
    {
        size_t AlphaIdx = Phase.FieldsProperties[alpha->index].Phase;
        double locOmegaAlpha = omega(AlphaIdx)(i,j,k);
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            size_t BetaIdx = Phase.FieldsProperties[beta->index].Phase;
            if (AlphaIdx != BetaIdx)
            {
                double locOmegaBeta  = omega(BetaIdx)(i,j,k);
                double locOmegaDelta = locOmegaBeta - locOmegaAlpha ;
                double loc_dPhi_dt   = locOmegaDelta*IP.Properties(i,j,k).get_mobility(alpha->index, beta->index)*Prefactor;
                Phase.FieldsDot(i,j,k).add_asym2(alpha->index, beta->index, loc_dPhi_dt);
            }
        }
    }
}
void GrandPotentialSolver::MergeLocalIncrements(long i, long j, long k, double dt)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ChemicalPotential   (i,j,k,{comp}) += dt*ChemicalPotentialDot(i,j,k,{comp});
        ChemicalPotentialDot(i,j,k,{comp}) = 0.0;
    }
}
void GrandPotentialSolver::MergeLocalIncrements1(long i, long j, long k, double dt)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ChemicalPotentialOld(i,j,k,{comp})  = ChemicalPotential(i,j,k,{comp});
        ChemicalPotential   (i,j,k,{comp}) += dt*ChemicalPotentialDot(i,j,k,{comp});
        ChemicalPotentialDot(i,j,k,{comp})  = 0.0;
    }
}
void GrandPotentialSolver::MergeLocalIncrements2(long i, long j, long k, double dt)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ChemicalPotentialDot (i,j,k,{comp})  = dt*ChemicalPotentialDot2(i,j,k,{comp});
        ChemicalPotential    (i,j,k,{comp}) += ChemicalPotentialDot(i,j,k,{comp});
        ChemicalPotentialDot2(i,j,k,{comp})  = 0.0;
    }
}
void GrandPotentialSolver::MergeLocalIncrements2Implicit(long i, long j, long k, PhaseField& Phase, const GrandPotentialDensity& omega, const InterfaceProperties& IP, const double dt)
{
    if (Ncomp == 1)
    {
        // This method uses the secant method (Newton's method with numeric calculation of the gradient) to solve the semi implicit time-step
        double ChemicalPotential_Curvature = ChemicalPotential(i,j,k,{0});
        auto residual = [this,i,j,k,&Phase,&omega,&IP,dt,ChemicalPotential_Curvature](double mu)
        {
            double mu0 = ChemicalPotential(i,j,k,{0});
            ChemicalPotential(i,j,k,{0}) = mu;
            ChemicalPotentialDot2(i,j,k,{0}) = 0.0;

            CalculateLocalConcentrations       (i,j,k,Phase,omega);
            CalculateLocalPhaseFieldIncrements (i,j,k,Phase,omega,IP);
            CalculateLocalIncrements2          (i,j,k,Phase,omega,dt);

            double res = ChemicalPotential(i,j,k,{0}) - ChemicalPotential_Curvature - ChemicalPotentialDot2(i,j,k,{0})*dt;
            ChemicalPotential(i,j,k,{0}) = mu0;
            return res;
        };

        // Use chemical potential of the previous time step as starting point for the root fining algorithm
        ChemicalPotential(i,j,k,{0}) = ChemicalPotentialOld(i,j,k,{0});

        try
        {
            double& x0 = ChemicalPotential(i,j,k,{0});
            double  x1 = ChemicalPotential(i,j,k,{0})*1.01;
            RootFindingAlgorithms::Secant(residual, x0, x1, ChemicalPotentialAccuracy, MaxIterations);
            //int error = RootFindingAlgorithms::Broyden(residual, ChemicalPotential(i,j,k)({0}), ChemicalPotentialAccuracy, MaxIterations);
        }
        catch (std::runtime_error& ecep)
        {
            ConsoleOutput::WriteWarning(ecep.what(),thisclassname,"MergeLocalIncrements2Implicit");
        }
    }
    else
    {
        // This method uses the Broyden's method to solve the semi implicit
        // time-step (Newton's method for multiple dimensions where the
        // inverse Jacobian is calculated numerically only for the first iteration).
        typedef std::vector<double> vector;
        vector ChemicalPotential_Curvature(Ncomp);
        for (std::size_t comp = 0; comp < Ncomp; comp++)
        {
            ChemicalPotential_Curvature[comp] = ChemicalPotential(i,j,k,{comp});
        }
        auto resiudal = [this,i,j,k,&Phase,&omega,&IP,dt,&ChemicalPotential_Curvature] (vector& x0)
        {
            std::vector<double> mu (Ncomp);
            std::vector<double> f0 (Ncomp);
            for (std::size_t comp = 0; comp < Ncomp; comp++)
            {
                mu[comp] = ChemicalPotential (i,j,k,{comp});
                ChemicalPotential     (i,j,k,{comp}) = x0[comp];
                ChemicalPotentialDot2 (i,j,k,{comp}) = 0;
            }

            CalculateLocalConcentrations       (i,j,k,Phase,omega);
            CalculateLocalPhaseFieldIncrements (i,j,k,Phase,omega,IP);
            CalculateLocalIncrements2          (i,j,k,Phase,omega,dt);

            for (std::size_t comp = 0; comp < Ncomp; comp++)
            {
                f0[comp] = ChemicalPotential(i,j,k,{comp}) - ChemicalPotential_Curvature[comp] - ChemicalPotentialDot2(i,j,k,{comp})*dt;
                ChemicalPotential (i,j,k,{comp}) = mu[comp];
            }
            return f0;
        };

        // Use chemical potential of the previous time step as starting point for the root fining algorithm
        std::vector<double> mu (Ncomp);
        for (std::size_t comp = 0; comp < Ncomp; comp++)
        {
            mu[comp] = ChemicalPotentialOld(i,j,k,{comp});
        }

        try
        {
            RootFindingAlgorithms::Broyden(resiudal, mu, ChemicalPotentialAccuracy, MaxIterations);
        }
        catch (std::runtime_error& ecep)
        {
            ConsoleOutput::WriteWarning(ecep.what(),thisclassname,"MergeLocalIncrements2Implicit");
        }

        for (std::size_t comp = 0; comp < Ncomp; comp++)
        {
            ChemicalPotential(i,j,k,{comp}) = mu[comp] ;
        }
    }
}
double GrandPotentialSolver::Mobility(long i, long j, long k, size_t comp, const PhaseField& Phase, const Temperature& Temp) const
{
    assert(Temp.Tx(i,j,k) > 0.0 && "Negative Temperature");
    double locMobility = 0;

    for (auto alpha = Phase.Fields(i,j,k).cbegin();
              alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
        locMobility += alpha->value*PhaseMobilities({PhaseIdx,comp})
            *std::exp(-PhaseActivationEnergies({PhaseIdx,comp})/PhysicalConstants::R/Temp.Tx(i,j,k));
    }

    if (Use_InterfaceMobilities)
    for (auto alpha = Phase.Fields(i,j,k).cbegin();
              alpha != Phase.Fields(i,j,k).cend(); alpha++)
    for (auto beta = alpha+1;
              beta  != Phase.Fields(i,j,k).cend();  beta++)
    {
        size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
        {
            size_t BetaIdx = Phase.FieldsProperties[beta->index].Phase;
            locMobility += 4*alpha->value*beta->value*InterfaceMobilities({PhaseIdx,BetaIdx,comp});
        }
    }

    if (Use_dPhaseMobility_dConcentration)
    for(size_t compB = 0; compB < Ncomp; compB++)
    for (auto alpha = Phase.Fields(i,j,k).cbegin();
              alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsProperties[alpha->index].Phase;
        locMobility += alpha->value*dPhaseMobilities_dConcentration({PhaseIdx,comp,compB});
    }

    return locMobility;
}
double GrandPotentialSolver::MaximumTimeStep(const PhaseField& Phase, const Temperature& Temp, const GrandPotentialDensity& omega)
{
    const double factor = (Grid.Active() > 1) ? ((Grid.Active() > 2) ? 6 : 4) : 2; // factor is derived stability criterion (div(grad(mu))
    double dt_max = std::numeric_limits<double>::max();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,reduction(min:dt_max))
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            double M = std::abs(Mobility(i,j,k,comp,Phase,Temp)/Susceptibility(i,j,k,comp,Phase,omega));
            double loc_dt_max = Grid.dx*Grid.dx/factor/M;
            dt_max = std::min(loc_dt_max, dt_max);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return dt_max;
}
void GrandPotentialSolver::SolveExplicit(PhaseField& Phase, const GrandPotentialDensity& omega, const BoundaryConditions&  BC, const InterfaceProperties& IP, const Temperature& Temp, const double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        CalculateLocalIncrements(i,j,k,Phase,omega,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        MergeLocalIncrements(i,j,k,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void GrandPotentialSolver::SolveImplicit(PhaseField& Phase, const GrandPotentialDensity& omega, const BoundaryConditions& BC, const InterfaceProperties& IP, const Temperature& Temp, const double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        CalculateLocalIncrements1(i,j,k,Phase,omega,dt);
        //CalculateLocalIncrements2(i,j,k,Phase,omega,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        MergeLocalIncrements1(i,j,k,dt);
        MergeLocalIncrements2Implicit(i,j,k,Phase,omega,IP,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void GrandPotentialSolver::Solve(PhaseField& Phase, const GrandPotentialDensity& omega, const BoundaryConditions& BC, const InterfaceProperties& IP, const Temperature& Temp, double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-1,)
    {
        CalculateLocalPhaseFieldIncrements(i,j,k,Phase,omega,IP);
        CaclulateLocalDiffusionFlux       (i,j,k,{0.0,0.0,0.0},Phase,Temp);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (UseImplicitSolver) SolveImplicit(Phase, omega, BC, IP, Temp, dt);
    else                   SolveExplicit(Phase, omega, BC, IP, Temp, dt);

    if (ConserveTOC) EnforceConservationOfTOC(Phase,omega);
    else
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,)
        {
            CalculateLocalConcentrations (i,j,k,Phase,omega);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    SetBoundaryConditions(BC);
}
double GrandPotentialSolver::MolarVolume(long i, long j, long k, const PhaseField& Phase, const GrandPotentialDensity& omega) const
{
    double locTotalConcentration = 0.0;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        locTotalConcentration += Concentration(i,j,k,comp,Phase,omega);
    }
    if (locTotalConcentration < 0 ) locTotalConcentration = 0;
    double locMolarVolume = (locTotalConcentration > 0.0) ? 1.0/locTotalConcentration : std::numeric_limits<float>::max();
    return locMolarVolume;
}
double GrandPotentialSolver::MoleFraction(long i, long j, long k, size_t comp, const PhaseField& Phase, const GrandPotentialDensity& omega) const
{
     double value = MolarVolume(i,j,k,Phase,omega)*Concentration(i,j,k,comp,Phase,omega);
     if (value > 1) value = 1;
     if (value < 0) value = 0;
     return value;
}
void GrandPotentialSolver::WriteVTK(const Settings& locSettings, const PhaseField& Phase, const GrandPotentialDensity& omega, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){"Chemical Potential "+ElementNames[comp]+" [J]", [this, comp                ] (long i, long j, long k) {return ChemicalPotential (i,j,k,{comp});}});
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+ " [kg/m^3]"                , [this, comp, &Phase, &omega] (long i, long j, long k) {return MolarMasses[comp]*Concentrations(i,j,k,{comp});}});
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+ " [mol/mol]"               , [this, comp, &Phase, &omega] (long i, long j, long k) {return MoleFraction      (i,j,k,comp,Phase,omega);}});
    }
    ListOfFields.push_back((VTK::Field_t){"Molar Volume [m^3/mol]",    [this, &Phase, &omega] (long i, long j, long k) {return MolarVolume(i,j,k,Phase,omega);}});
    ListOfFields.push_back((VTK::Field_t){"Grand Potential GrandPotentialDensity [J/m^3]", [&Phase, &omega] (long i, long j, long k) {return omega(i,j,k,Phase);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + "_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void GrandPotentialSolver::WriteVTKChemicalPotential(const Settings& locSettings, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){"Chemical Potential "+ElementNames[comp]+" [J]", [this, comp] (long i, long j, long k) {return ChemicalPotential (i,j,k,{comp});}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ChemicalPotential_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void GrandPotentialSolver::WriteVTKConcentration(const Settings& locSettings, const PhaseField& Phase, const GrandPotentialDensity& omega, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+" [1/m^3]" , [this, comp, &Phase, &omega] (long i, long j, long k) {return Concentration(i,j,k,comp,Phase,omega);}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Concentration_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void GrandPotentialSolver::WriteVTKMoleFraction(const Settings& locSettings, const PhaseField& Phase, const GrandPotentialDensity& omega, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+ " [mol/mol]", [this, comp, &Phase, &omega] (long i, long j, long k) {return MoleFraction(i,j,k,comp,Phase,omega);}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "MoleFraction_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void GrandPotentialSolver::WriteVTKMolarVolume(const Settings& locSettings, const PhaseField& Phase, const GrandPotentialDensity& omega, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Molar Volume [mole/m^3]" , [this, &Phase, &omega] (long i, long j, long k) {return MolarVolume(i,j,k,Phase,omega);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "MolarVolume_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
//void GrandPotentialSolver::WriteVTKMobility(const Settings& locSettings, long tStep, long precision) const
//{
//    std::vector<VTK::Field_t> ListOfFields;
//    for(size_t comp = 0; comp < Ncomp; comp++)
//    {
//        ListOfFields.push_back((VTK::Field_t){"Mobility_"+ElementNames[comp]+" [(m^3 x)/kg]", [this,comp] (long i, long j, long k) {return Mobilities(i,j,k)({comp});}});
//    }
//    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Mobility_", tStep, ".vts");
//    VTK::Write(Filename, locSettings, ListOfFields, precision);
//}
void GrandPotentialSolver::WriteVTKDiffusionFlux(const Settings& locSettings, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"DiffusionFlux_"+ElementNames[comp]+" [1/(m^2 s)]", [this,comp] (long i, long j, long k) {return DiffusionFlux(i,j,k,{comp});}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "DiffusionFlux_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void GrandPotentialSolver::WriteVTKVelocity(const Settings& locSettings, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"Velocity_"+ElementNames[comp]+" [1/(m^2 s)]", [this,comp] (long i, long j, long k) {return (Concentrations(i,j,k,{comp}) > 0.0) ? DiffusionFlux(i,j,k,{comp})/Concentrations(i,j,k,{comp}) : dVector3({0.0,0.0,0.0});}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Velocity_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void GrandPotentialSolver::WriteVTKSusceptibility(const Settings& locSettings, const PhaseField& Phase, const GrandPotentialDensity& omega, long tStep, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"Susceptibility_"+ElementNames[comp]+" [(m^3 x)/kg]", [this,comp,&Phase,&omega] (long i, long j, long k) {return Susceptibility(i,j,k,comp,Phase,omega);}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Susceptibility_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
//void GrandPotentialSolver::WriteVTKDiffusivity(const Settings& locSettings, const PhaseField& Phase, const GrandPotentialDensity& omega, long tStep, long precision) const
//{
//    std::vector<VTK::Field_t> ListOfFields;
//    for(size_t comp = 0; comp < Ncomp; comp++)
//    {
//        ListOfFields.push_back((VTK::Field_t){"Diffusivity_"+ElementNames[comp]+" [(m^3 x)/kg]", [this,comp,&Phase,&omega] (long i, long j, long k) {return Mobilities(i,j,k)({comp})/Susceptibility(i,j,k,comp,Phase,omega);}});
//    }
//    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Diffusivity_", tStep, ".vts");
//    VTK::Write(Filename, locSettings, ListOfFields, precision);
//}
void GrandPotentialSolver::WriteVTKCenterOfMassVelocity(const Settings& locSettings, const PhaseField& Phase, long tStep, long precision)
{
    CalculateCenterOfMassVelocitiesOfGrains(Phase);
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Rigid Body Velocity [m/s]", [this, &Phase] (long i, long j, long k) {return CenterOfMassVelocity(i,j,k,Phase);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "CenterOfMassVelocity_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
}// namespace openphase::GrandPotential
