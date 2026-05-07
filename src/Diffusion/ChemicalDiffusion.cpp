/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 *  File created :   2025
 *  Main contributors :   Oleg Shchyglo; Matthias Stratmann; Efim Borukhovich
 *
 */

#include "Diffusion/ChemicalDiffusion.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Composition.h"
#include "Temperature.h"
#include "BoundaryConditions.h"
#include "Diffusion/DiffusionProperties.h"
#include "Thermodynamics/ThermodynamicPropertiesEQP.h"
#include "Tools/TimeInfo.h"

namespace openphase
{

using namespace std;

enum class DiffusionModels
{
    Diagonal,
    DICTRA,
    Pairwise
};

enum class PositiveDiffusivityModes
{
    None,
    Diagonal,
    CrossTerms,
    All
};

inline pair<int, int> three_state_bounds_selector(int selector, int lower_bound, int upper_bound)
{
    pair<int, int> result;
    switch (selector)
    {
        case -1: // lower bound only
        {
            result = make_pair(lower_bound, lower_bound + 1);
            break;
        }
        case 0: // both bounds
        {
            result = make_pair(lower_bound, upper_bound);
            break;
        }
        case 1: // upper bound only
        {
            result = make_pair(upper_bound - 1, upper_bound);
            break;
        }
    }
    return result;
}

inline bool model_component_switch(DiffusionModels model, size_t comp_idx, size_t grad_idx)
{
    bool answer = false;

    switch(model)
    {
        case DiffusionModels::Diagonal:
        {
            answer = (comp_idx == grad_idx);
            break;
        }
        case DiffusionModels::DICTRA:
        {
            answer = true;
            break;
        }
        case DiffusionModels::Pairwise:
        {
            answer = (comp_idx != grad_idx);
            break;
        }
    }
    return answer;
}

inline bool positive_diffusivity_switch(PositiveDiffusivityModes mode, size_t comp_idx, size_t grad_idx)
{
    bool answer = false;

    switch(mode)
    {
        case PositiveDiffusivityModes::None:
        {
            answer = false;
            break;
        }
        case PositiveDiffusivityModes::Diagonal:
        {
            answer = (comp_idx == grad_idx);
            break;
        }
        case PositiveDiffusivityModes::CrossTerms:
        {
            answer = (comp_idx != grad_idx);
            break;
        }
        case PositiveDiffusivityModes::All:
        {
            answer = true;
            break;
        }
    }
    return answer;
}

class OP_EXPORTS ChemicalDiffusion_Impl : public OPObject                       ///<  Chemical diffusion solver
{
 public:
    ChemicalDiffusion_Impl();
    ~ChemicalDiffusion_Impl();

    ChemicalDiffusion_Impl(Settings& locSettings,
                         const std::string InputFileName = DefaultInputFileName)///< Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from the input stream

    void Solve(PhaseField& Phi, Composition& Cx,Temperature& Tx,
               DiffusionProperties& DP,
               BoundaryConditions& BC,double dt);                               ///< Solves chemical diffusion problem
    void SolveSF(PhaseField& Phi, Composition& Cx,Temperature& Tx,
               DiffusionProperties& DP,
               BoundaryConditions& BC,double dt);                               ///< Solves chemical diffusion problem

    void Solve(PhaseField& Phi, Composition& Cx,Temperature& Tx,
               DiffusionProperties& DP, ThermodynamicPropertiesEQP& TP,
               BoundaryConditions& BC,double dt);                               ///< Solves chemical diffusion problem

    void PrintStatistics(const Settings& locSettings, const int tStep);         ///< Prints solver execution statistics

 protected:
 private:

    // Iterative solver
    void SolveIteratively(PhaseField& Phi, Composition& Cx,Temperature& Tx,
                              DiffusionProperties& DP, ThermodynamicPropertiesEQP& TP,
                              BoundaryConditions& BC,double dt);
    void SolveIterativelySF(PhaseField& Phi, Composition& Cx,Temperature& Tx,
                              DiffusionProperties& DP, 
                              BoundaryConditions& BC,double dt);
    void SolveIteratively(PhaseField& Phi, Composition& Cx,Temperature& Tx,
                              DiffusionProperties& DP, 
                              BoundaryConditions& BC,double dt);
    // Calculation of increments
    void ClearIncrements(Composition& Cx);
    void HomogenizeFastDiffusors(PhaseField& Phi, Composition& Cx);
    void CalculateDiffusionIncrements(PhaseField& Phi,
                                      Composition& Cx,
                                      Temperature& Tx,
                                      DiffusionProperties& DP,
                                      double dt);
    void CalculateDiffusionIncrementsSF(PhaseField& Phi,
                                      Composition& Cx,
                                      Temperature& Tx,
                                      DiffusionProperties& DP,
                                      double dt);


    void CalculateChemicalPotentialContribution(PhaseField& Phi,
                                                Composition& Cx,
                                                Temperature& Tx,
                                                DiffusionProperties& DP,
                                                ThermodynamicPropertiesEQP& TP,
                                                double dt);

    void Solve1Dextension(Composition& Cx,
                          Temperature& Tx,
                          DiffusionProperties& DP,
                          Composition1Dextension& CxExt,
                          BoundaryConditionTypes& lextBC,
                          double dt);

    void CalculateAntitrappingIncrements(PhaseField& Phi,
                                         Composition& Cx,
                                         Temperature& Tx,
                                         DiffusionProperties& DP,
                                         double dt);
    // Increments limiting to prevent out of range compositions
    bool CalculateLimits(PhaseField& Phi, Composition& Cx, double dt);

    void LimitDiffusionIncrements(PhaseField& Phi,
                                  Composition& Cx,
                                  Temperature& Tx,
                                  DiffusionProperties& DP,
                                  double dt);

    void LimitChemicalPotentialContribution(PhaseField& Phi,Composition& Cx, Temperature& Tx,
                                 DiffusionProperties& DP,
                                 ThermodynamicPropertiesEQP& TP,
                                 double dt);

    void LimitAntitrappingIncrements(PhaseField& Phi,
                                  Composition& Cx, Temperature& Tx,
                                  DiffusionProperties& DP,
                                  double dt);
    // Updating composition
    void ApplyIncrements(PhaseField& Phi, Composition& Cx, double dt);
    void ApplyIncrementsSF(PhaseField& Phi, Composition& Cx, double dt);

    bool OPDIFF;
    bool ConsiderChemicalPotential;                                             ///< True if chemical potential contribution should be considered
    bool EnableExecutionTiming;                                                 ///< True if diffusion solver execution timing is enabled
    bool EnableAntitrapping;                                                    ///< True if antitrapping current calculation is enabled
    bool EnableLimitsChecking;                                                  ///< True if composition limits checking is enabled

    GridParameters Grid;                                                        ///< Simulation grid parameters

    DiffusionModels Model;                                                      ///< Diffusion model selector
    PositiveDiffusivityModes Mode;                                              ///< Positive diffusion coefficients mode selector
    LaplacianStencils DiffusionStencil;                                         ///< Diffusion stencil selector
    LaplacianStencil DStencil;                                                  ///< Diffusion stencil. Uses Laplacian stencil as the basis

    TimeInfo Timer;
};

ChemicalDiffusion::ChemicalDiffusion() : impl_(new ChemicalDiffusion_Impl) {}
ChemicalDiffusion::~ChemicalDiffusion() = default;

ChemicalDiffusion::ChemicalDiffusion(Settings& locSettings,
                                     const std::string InputFileName)
                                     : impl_(new ChemicalDiffusion_Impl)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void ChemicalDiffusion::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    impl_->Initialize(locSettings, ObjectNameSuffix);
    thisclassname = impl_->thisclassname;
    thisobjectname = impl_->thisobjectname;
    initialized = impl_->initialized;
}

void ChemicalDiffusion::ReadInput(const std::string InputFileName)
{
    impl_->ReadInput(InputFileName);
}

void ChemicalDiffusion::ReadInput(std::stringstream& inp)
{
    impl_->ReadInput(inp);
}

void ChemicalDiffusion::Solve(PhaseField& Phi, Composition& Cx,
            Temperature& Tx,
            DiffusionProperties& DP, 
            BoundaryConditions& BC, double dt)
{
    impl_->Solve(Phi, Cx, Tx, DP, BC, dt);
}

void ChemicalDiffusion::SolveSF(PhaseField& Phi, Composition& Cx,
            Temperature& Tx,
            DiffusionProperties& DP, 
            BoundaryConditions& BC, double dt)
{
    impl_->SolveSF(Phi, Cx, Tx, DP, BC, dt);
}

void ChemicalDiffusion::Solve(PhaseField& Phi, Composition& Cx,
            Temperature& Tx,
            DiffusionProperties& DP, ThermodynamicPropertiesEQP& TP,
            BoundaryConditions& BC, double dt)
{
    impl_->Solve(Phi, Cx, Tx, DP, TP, BC, dt);
}

void ChemicalDiffusion::PrintStatistics(const Settings& locSettings,
                             const int tStep)
{
    impl_->PrintStatistics(locSettings,tStep);
}

ChemicalDiffusion_Impl::ChemicalDiffusion_Impl()
{
    OPDIFF = false;
    EnableAntitrapping = false;
    ConsiderChemicalPotential = false;
    EnableExecutionTiming = false;
    EnableLimitsChecking = false;
    Mode = PositiveDiffusivityModes::None;
    Model = DiffusionModels::DICTRA;
    DiffusionStencil = LaplacianStencils::Isotropic;
};

ChemicalDiffusion_Impl::~ChemicalDiffusion_Impl()
{

};

void ChemicalDiffusion_Impl::PrintStatistics(const Settings& locSettings,
                             const int tStep)
{
    if(EnableExecutionTiming)
    {
        Timer.PrintWallClockSummary();
    }
}

void ChemicalDiffusion_Impl::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ChemicalDiffusion";
    thisobjectname = thisclassname + ObjectNameSuffix;

    EnableAntitrapping = false;
    ConsiderChemicalPotential = false;
    EnableExecutionTiming = false;
    Mode = PositiveDiffusivityModes::None;
    Model = DiffusionModels::DICTRA;
    DiffusionStencil = LaplacianStencils::Isotropic;

    Grid = locSettings.Grid;

    switch(Grid.Active())
    {
        case 1:
        {
            DStencil.SetNoCenter(LaplacianStencil1D_3, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            break;
        }
        case 2:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil2D_5, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil2D_9, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            break;
        }
        case 3:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil3D_7, Grid.dx);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil3D_27a, Grid.dx);
            }
            break;
        }
        default:
        {

        }
    }

    Timer.Initialize(locSettings, thisclassname+"::Solve() execution time statistics");
    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ChemicalDiffusion_Impl::ReadInput(const string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    ReadInput(data);

    ConsoleOutput::WriteLine();
}

void ChemicalDiffusion_Impl::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    string tmpModel = FileInterface::ReadParameterK(inp, moduleLocation, "DiffusionModel", false, "DICTRA");

    if (tmpModel == "DIAGONAL")
    {
        Model = DiffusionModels::Diagonal;
    }
    else if (tmpModel == "DICTRA")
    {
        Model = DiffusionModels::DICTRA;
    }
    else if (tmpModel == "PAIRWISE")
    {
        Model = DiffusionModels::Pairwise;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong diffusion model specified!\nThe default \"DICTRA\" model is used!", thisclassname, "ReadInput()");
    }

    string tmpDiff = FileInterface::ReadParameterK(inp, moduleLocation, "PositiveDiffusivity", false, "NONE");

    if (tmpDiff == "NONE")
    {
        Mode = PositiveDiffusivityModes::None;
    }
    else if (tmpDiff == "DIAGONAL")
    {
        Mode = PositiveDiffusivityModes::Diagonal;
    }
    else if (tmpDiff == "CROSSTERMS")
    {
        Mode = PositiveDiffusivityModes::CrossTerms;
    }
    else if (tmpDiff == "ALL")
    {
        Mode = PositiveDiffusivityModes::All;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong positive diffusion coefficients setting specified!\nThe default \"NONE\" model is used!", thisclassname, "ReadInput()");
    }

    string tmp1 = FileInterface::ReadParameterK(inp, moduleLocation, string("DiffusionStencil"), false, string("ISOTROPIC"));
    if(tmp1 == "SIMPLE")
    {
        DiffusionStencil = LaplacianStencils::Simple;
    }
    else if(tmp1 == "ISOTROPIC")
    {
        DiffusionStencil = LaplacianStencils::Isotropic;
    }
    else if(tmp1 == "LB")
    {
        DiffusionStencil = LaplacianStencils::LB;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong diffusion stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadInput()");
    }

    switch(Grid.Active())
    {
        case 1:
        {
            DStencil.SetNoCenter(LaplacianStencil1D_3, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            break;
        }
        case 2:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil2D_5, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil2D_9, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            break;
        }
        case 3:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil3D_7, Grid.dx);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil3D_27a, Grid.dx);
            }
            break;
        }
        default:
        {

        }
    }

    EnableAntitrapping = FileInterface::ReadParameterB(inp,moduleLocation, string("EnableAntitrapping"),false,false);
    ConsiderChemicalPotential = FileInterface::ReadParameterB(inp, moduleLocation, "ConsiderChemicalPotential", false, false);
    EnableExecutionTiming = FileInterface::ReadParameterB(inp, moduleLocation, "EnableExecutionTiming", false, false);
    EnableLimitsChecking = FileInterface::ReadParameterB(inp, moduleLocation, "EnableLimitsChecking", false, false);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void ChemicalDiffusion_Impl::Solve(PhaseField& Phi,
                                  Composition& Cx,
                                  Temperature& Tx,
                                  DiffusionProperties& DP,
                                  BoundaryConditions& BC,
                                  double dt)
{
    SolveIteratively(Phi,Cx,Tx,DP,BC,dt);
}

void ChemicalDiffusion_Impl::SolveSF(PhaseField& Phi,
                                  Composition& Cx,
                                  Temperature& Tx,
                                  DiffusionProperties& DP,
                                  BoundaryConditions& BC,
                                  double dt)
{

    SolveIterativelySF(Phi,Cx,Tx,DP,BC,dt);
}

void ChemicalDiffusion_Impl::Solve(PhaseField& Phi,
                                  Composition& Cx,
                                  Temperature& Tx,
                                  DiffusionProperties& DP,
                                  ThermodynamicPropertiesEQP& TP,
                                  BoundaryConditions& BC,
                                  double dt)
{
    SolveIteratively(Phi,Cx,Tx,DP,TP,BC,dt);
}

void ChemicalDiffusion_Impl::SolveIteratively(PhaseField& Phi,
                                              Composition& Cx,
                                              Temperature& Tx,
                                              DiffusionProperties& DP,
                                              BoundaryConditions& BC,
                                              double dt)
{
    /* Solves chemical diffusion in all phases iteratively*/

    // Calculate new (safe) time step

    double new_dt = dt;
    int Niterations = 1;

    if(DP.MaxDiffusionCoefficient != 0.0)
    {
        double max_time_step = 0.16*Phi.Grid.dx*Phi.Grid.dx/DP.MaxDiffusionCoefficient;
        int number_substeps  = max(1.0,ceil(dt/max_time_step));

        new_dt  = dt/double(number_substeps);
        Niterations = number_substeps;
    }

    if(Niterations > 10)
    {
        string message = "Too many internal diffusion iterations -> " + std::to_string(Niterations) + ".\n"
                       + "Adjust simulation time step to avoid simulation instability.";
        ConsoleOutput::WriteWarning(message, thisclassname, "Solve()");
    }

    // Iterate

    for(int iteration = 0; iteration < Niterations; iteration++)
    {
        if(EnableExecutionTiming) Timer.SetStart();

        // Calculate Increments

        HomogenizeFastDiffusors(Phi,Cx);
        if(EnableExecutionTiming) Timer.SetTimeStamp("HomogenizeFastDiffusors");

        Cx.SetBoundaryConditions(BC);
        if(EnableExecutionTiming) Timer.SetTimeStamp("SetBoundaryConditions");

        CalculateDiffusionIncrements(Phi,Cx,Tx,DP,new_dt);
        if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateDiffusionIncrements");

        if(EnableAntitrapping)
        {
            CalculateAntitrappingIncrements(Phi,Cx,Tx,DP,new_dt);
        }
        if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateAntitrappingIncrements");

        // Limit Increments

        if(EnableLimitsChecking)
        {
            bool limiting_needed = CalculateLimits(Phi, Cx, new_dt);
            if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateLimits");

            if(limiting_needed)
            {
                Cx.SetLimitsBoundaryConditions(BC);
                if(EnableExecutionTiming) Timer.SetTimeStamp("SetLimitsBoundaryConditions");

                LimitDiffusionIncrements(Phi, Cx, Tx, DP, new_dt);
                if(EnableExecutionTiming) Timer.SetTimeStamp("LimitDiffusionIncrements");


                if(EnableAntitrapping)
                {
                    LimitAntitrappingIncrements(Phi, Cx, Tx, DP, new_dt);
                }
                if(EnableExecutionTiming) Timer.SetTimeStamp("LimitAntitrappingIncrements");
            }
        }
        // Apply Increments

        ApplyIncrements(Phi,Cx,new_dt);
        if(EnableExecutionTiming) Timer.SetTimeStamp("ApplyIncrements");

        ClearIncrements(Cx);
        if(EnableExecutionTiming) Timer.SetTimeStamp("ClearIncrements");

        Cx.SetBoundaryConditions(BC);
        if(EnableExecutionTiming) Timer.SetTimeStamp("SetBoundaryConditions");

        // Solve diffusion in 1D extensions

    #ifdef MPI_PARALLEL
        if(MPI_RANK == 0)
    #endif
        if(Cx.ExtensionX0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionX0, BC.BC0X, new_dt);}
    #ifdef MPI_PARALLEL
        if(MPI_RANK == MPI_SIZE - 1)
    #endif
        if(Cx.ExtensionXN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionXN, BC.BCNX, new_dt);}

        if(Cx.ExtensionY0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionY0, BC.BC0Y, new_dt);}
        if(Cx.ExtensionYN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionYN, BC.BCNY, new_dt);}
        if(Cx.ExtensionZ0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionZ0, BC.BC0Z, new_dt);}
        if(Cx.ExtensionZN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionZN, BC.BCNZ, new_dt);}
        if(EnableExecutionTiming) Timer.SetTimeStamp("Solve1Dextensions");
    }
}

void ChemicalDiffusion_Impl::SolveIterativelySF(PhaseField& Phi,
                                              Composition& Cx,
                                              Temperature& Tx,
                                              DiffusionProperties& DP,
                                              BoundaryConditions& BC,
                                              double dt)
{
    /* Solves chemical diffusion in all phases iteratively*/

    // Calculate new (safe) time step

    double new_dt = dt;
    int Niterations = 1;

    if(DP.MaxDiffusionCoefficient != 0.0)
    {
        double max_time_step = 0.16*Phi.Grid.dx*Phi.Grid.dx/DP.MaxDiffusionCoefficient;
        int number_substeps  = max(1.0,ceil(dt/max_time_step));

        new_dt  = dt/double(number_substeps);
        Niterations = number_substeps;
    }

    if(Niterations > 10)
    {
        string message = "Too many internal diffusion iterations -> " + std::to_string(Niterations) + ".\n"
                       + "Adjust simulation time step to avoid simulation instability.";
        ConsoleOutput::WriteWarning(message, thisclassname, "Solve()");
    }

    // Iterate

    for(int iteration = 0; iteration < Niterations; iteration++)
    {
        if(EnableExecutionTiming) Timer.SetStart();

        // Calculate Increments

        HomogenizeFastDiffusors(Phi,Cx);
        if(EnableExecutionTiming) Timer.SetTimeStamp("HomogenizeFastDiffusors");

        Cx.SetBoundaryConditions(BC);
        if(EnableExecutionTiming) Timer.SetTimeStamp("SetBoundaryConditions");

        CalculateDiffusionIncrementsSF(Phi,Cx,Tx,DP,new_dt);
        if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateDiffusionIncrements");

        if(EnableAntitrapping)
        {
            CalculateAntitrappingIncrements(Phi,Cx,Tx,DP,new_dt);
        }
        if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateAntitrappingIncrements");

        // Limit Increments

        if(EnableLimitsChecking)
        {
            bool limiting_needed = CalculateLimits(Phi, Cx, new_dt);
            if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateLimits");

            if(limiting_needed)
            {
                Cx.SetLimitsBoundaryConditions(BC);
                if(EnableExecutionTiming) Timer.SetTimeStamp("SetLimitsBoundaryConditions");

                LimitDiffusionIncrements(Phi, Cx, Tx, DP, new_dt);
                if(EnableExecutionTiming) Timer.SetTimeStamp("LimitDiffusionIncrements");


                if(EnableAntitrapping)
                {
                    LimitAntitrappingIncrements(Phi, Cx, Tx, DP, new_dt);
                }
                if(EnableExecutionTiming) Timer.SetTimeStamp("LimitAntitrappingIncrements");
            }
        }
        // Apply Increments

        ApplyIncrementsSF(Phi,Cx,new_dt);
        if(EnableExecutionTiming) Timer.SetTimeStamp("ApplyIncrements");

        ClearIncrements(Cx);
        if(EnableExecutionTiming) Timer.SetTimeStamp("ClearIncrements");

        Cx.SetBoundaryConditions(BC);
        if(EnableExecutionTiming) Timer.SetTimeStamp("SetBoundaryConditions");

        // Solve diffusion in 1D extensions

    #ifdef MPI_PARALLEL
        if(MPI_RANK == 0)
    #endif
        if(Cx.ExtensionX0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionX0, BC.BC0X, new_dt);}
    #ifdef MPI_PARALLEL
        if(MPI_RANK == MPI_SIZE - 1)
    #endif
        if(Cx.ExtensionXN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionXN, BC.BCNX, new_dt);}

        if(Cx.ExtensionY0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionY0, BC.BC0Y, new_dt);}
        if(Cx.ExtensionYN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionYN, BC.BCNY, new_dt);}
        if(Cx.ExtensionZ0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionZ0, BC.BC0Z, new_dt);}
        if(Cx.ExtensionZN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionZN, BC.BCNZ, new_dt);}
        if(EnableExecutionTiming) Timer.SetTimeStamp("Solve1Dextensions");
    }
}

void ChemicalDiffusion_Impl::SolveIteratively(PhaseField& Phi,
                                              Composition& Cx,
                                              Temperature& Tx,
                                              DiffusionProperties& DP,
                                              ThermodynamicPropertiesEQP& TP,
                                              BoundaryConditions& BC,
                                              double dt)
{
    /* Solves chemical diffusion in all phases iteratively*/

    // Calculate new (safe) time step

    double new_dt = dt;
    int Niterations = 1;

    if(DP.MaxDiffusionCoefficient != 0.0)
    {
        double max_time_step = 0.16*Phi.Grid.dx*Phi.Grid.dx/DP.MaxDiffusionCoefficient;
        int number_substeps  = max(1.0,ceil(dt/max_time_step));

        new_dt  = dt/double(number_substeps);
        Niterations = number_substeps;
    }

    if(Niterations > 10)
    {
        string message = "Too many internal diffusion iterations -> " + std::to_string(Niterations) + ".\n"
                       + "Adjust simulation time step to avoid simulation instability.";
        ConsoleOutput::WriteWarning(message, thisclassname, "Solve()");
    }

    // Iterate

    for(int iteration = 0; iteration < Niterations; iteration++)
    {
        if(EnableExecutionTiming) Timer.SetStart();

        // Calculate Increments

        HomogenizeFastDiffusors(Phi,Cx);
        if(EnableExecutionTiming) Timer.SetTimeStamp("HomogenizeFastDiffusors");

        Cx.SetBoundaryConditions(BC);
        if(EnableExecutionTiming) Timer.SetTimeStamp("SetBoundaryConditions");

        CalculateDiffusionIncrements(Phi,Cx,Tx,DP,new_dt);
        if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateDiffusionIncrements");

        if(ConsiderChemicalPotential)
        {
            CalculateChemicalPotentialContribution(Phi,Cx,Tx,DP,TP,new_dt);
        }
        if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateChemicalPotentialContribution");

        if(EnableAntitrapping)
        {
            CalculateAntitrappingIncrements(Phi,Cx,Tx,DP,new_dt);
        }
        if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateAntitrappingIncrements");

        // Limit Increments

        if(EnableLimitsChecking)
        {
            bool limiting_needed = CalculateLimits(Phi, Cx, new_dt);
            if(EnableExecutionTiming) Timer.SetTimeStamp("CalculateLimits");

            if(limiting_needed)
            {
                Cx.SetLimitsBoundaryConditions(BC);
                if(EnableExecutionTiming) Timer.SetTimeStamp("SetLimitsBoundaryConditions");

                LimitDiffusionIncrements(Phi, Cx, Tx, DP, new_dt);
                if(EnableExecutionTiming) Timer.SetTimeStamp("LimitDiffusionIncrements");

                if(ConsiderChemicalPotential)
                {
                    LimitChemicalPotentialContribution(Phi, Cx, Tx, DP, TP, new_dt);
                }
                if(EnableExecutionTiming) Timer.SetTimeStamp("LimitChemicalPotentialContribution");

                if(EnableAntitrapping)
                {
                    LimitAntitrappingIncrements(Phi, Cx, Tx, DP, new_dt);
                }
                if(EnableExecutionTiming) Timer.SetTimeStamp("LimitAntitrappingIncrements");
            }
        }
        // Apply Increments

        ApplyIncrements(Phi,Cx,new_dt);
        if(EnableExecutionTiming) Timer.SetTimeStamp("ApplyIncrements");

        ClearIncrements(Cx);
        if(EnableExecutionTiming) Timer.SetTimeStamp("ClearIncrements");

        if(ConsiderChemicalPotential)
        {
            TP.ClearChemicalPotential();
        }
        if(EnableExecutionTiming) Timer.SetTimeStamp("ClearChemicalPotential");

        Cx.SetBoundaryConditions(BC);
        if(EnableExecutionTiming) Timer.SetTimeStamp("SetBoundaryConditions");

        // Solve diffusion in 1D extensions

    #ifdef MPI_PARALLEL
        if(MPI_RANK == 0)
    #endif
        if(Cx.ExtensionX0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionX0, BC.BC0X, new_dt);}
    #ifdef MPI_PARALLEL
        if(MPI_RANK == MPI_SIZE - 1)
    #endif
        if(Cx.ExtensionXN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionXN, BC.BCNX, new_dt);}

        if(Cx.ExtensionY0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionY0, BC.BC0Y, new_dt);}
        if(Cx.ExtensionYN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionYN, BC.BCNY, new_dt);}
        if(Cx.ExtensionZ0.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionZ0, BC.BC0Z, new_dt);}
        if(Cx.ExtensionZN.isActive()) {Solve1Dextension(Cx, Tx, DP, Cx.ExtensionZN, BC.BCNZ, new_dt);}
        if(EnableExecutionTiming) Timer.SetTimeStamp("Solve1Dextensions");
    }
}

void ChemicalDiffusion_Impl::HomogenizeFastDiffusors(PhaseField& Phi, Composition& Cx)
{
    /** This function homogenizes the composition profile for certain elements
    in each phase, but not in each grain separately! The local difference is
    applied to the composition of the reference element. */

    for(size_t n = 0; n < Cx.Nphases; n++)
    {
        bool homogenize = false;

        for(size_t comp = 0; comp < Cx.Ncomp; comp++)
        if(!Cx.Phase[n].Component[comp].isReference and
            Cx.Phase[n].Component[comp].isFastDiffusor)
        {
            homogenize = true;
        }

        if(homogenize)
        {
            size_t ref = Cx.Phase[n].ReferenceElement;                                  // Reference element, against all other elements 'diffuse'
            Tensor<double,1> freec({Cx.Phase[n].Ncomp});                                // Mean value for mole fraction per phase
            double freep = 0.0;                                                         // Volume of phase 'n'

            Tensor<double,1> todistribute({Cx.Phase[n].Ncomp});                         // Composition that couldn't be distributed

            freec.set_to_zero();
            todistribute.set_to_zero();

            /* Calculate the mean composition of element 'comp' in phase 'n'. */

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,reduction(TensorD1Sum:freec) reduction(+:freep))
            if(Phi.Fractions(i,j,k,{n}) > DBL_EPSILON)
            {
                freep += Phi.Fractions(i,j,k,{n});

                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                if(Cx.Phase[n].Component[comp].isFastDiffusor and
                  !Cx.Phase[n].Component[comp].isReference)
                {
                    freec({comp}) += Cx.MoleFractions(i,j,k,{n,comp})*Phi.Fractions(i,j,k,{n});
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            Tensor<double,1> meancomp = freec;                                          // Mean composition to set each local composition value to

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                meancomp({comp}) /= freep;
            }

            /* Distribute the mean composition to every point of the same phase. */

            if(freep > DBL_EPSILON)
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,reduction(TensorD1Sum:todistribute))
            if(Phi.Fractions(i,j,k,{n}) > DBL_EPSILON)
            {
                Tensor<double,1> delta({Cx.Phase[n].Ncomp});
                delta.set_to_zero();

                /* Calculate local reference composition with all local
                composition changes. */

                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                if(Cx.Phase[n].Component[comp].isFastDiffusor and
                  !Cx.Phase[n].Component[comp].isReference)
                {
                    delta({comp}) = meancomp({comp}) - Cx.MoleFractions(i,j,k,{n,comp});
                }

                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                if(comp != ref)
                {
                    delta({ref}) -= delta({comp});
                }

                double refcomp = Cx.MoleFractions(i,j,k,{n,ref}) + delta({ref});

                /* Check if new reference composition is still within limits,
                otherwise change composition deltas. */

                if((refcomp < Cx.Phase[n].Component[ref].Min) and
                   (refcomp > Cx.Phase[n].Component[ref].Max))
                {
                    /* The composition difference can NOT be applied to the
                    reference element. The composition of component 'comp' can
                    therefore not be set to the mean value. */

                    Tensor<double,1> prevDelta = delta;

                    double deltaxref = delta({ref});
                    double xrefmax   = Cx.Phase[n].Component[ref].Max;
                    double xrefmin   = Cx.Phase[n].Component[ref].Min;
                    double xref      = Cx.MoleFractions(i,j,k,{n,ref});

                    if(refcomp < Cx.Phase[n].Component[ref].Min - DBL_EPSILON)
                    {
                        if(deltaxref < DBL_EPSILON and xref > xrefmin)
                        {
                            delta *= (xrefmin-xref)/deltaxref;
                        }
                        else
                        {
                            delta *= 0.0;
                        }
                    }
                    else if(refcomp > Cx.Phase[n].Component[ref].Max + DBL_EPSILON)
                    {
                        if(deltaxref > DBL_EPSILON and xref < xrefmax)
                        {
                            delta *= (xrefmax-xref)/deltaxref;
                        }
                        else
                        {
                            delta *= 0.0;
                        }
                    }

                    todistribute += (prevDelta - delta) * Phi.Fractions(i,j,k,{n});
                }

                /* The composition difference can (now) be applied to the
                reference element. */

                Cx.MoleFractionsTotal(i,j,k) += delta*Phi.Fractions(i,j,k,{n});

                /* Mass conservation error is stored, so it can later be
                applied to all other points again. */

                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                if(fabs(todistribute({comp})) > DBL_EPSILON)
                {
                    cout << "Mass conservation error for phase " << n << " and element "
                         << todistribute({comp}) << endl;
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
        }
    }
}

void ChemicalDiffusion_Impl::ClearIncrements(Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,Cx.MoleFractionsTotal.Bcells(),)
    {
        Cx.MoleFractionsDotIn(i,j,k).set_to_zero();
        Cx.MoleFractionsDotOut(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ChemicalDiffusion_Impl::CalculateDiffusionIncrements(PhaseField& Phi,
                          Composition& Cx, Temperature& Tx,
                          DiffusionProperties& DP, double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractions,0,)
    {
        for(size_t n = 0; n < Cx.Nphases; n++)
        {
            double locPhaseFraction = Phi.Fractions(i,j,k,{n});
            if (locPhaseFraction != 0.0)
            {
                for (size_t comp = 0; comp < Cx.Ncomp; comp++)
                if (!Cx.Phase[n].Component[comp].isReference and
                    !Cx.Phase[n].Component[comp].isFastDiffusor)
                {
                    double localDelta = 0.0;

                    for (size_t grad = 0; grad < Cx.Ncomp; grad++)
                    if(!Cx.Phase[n].Component[grad].isReference and
                       !Cx.Phase[n].Component[grad].isFastDiffusor and
                       model_component_switch(Model,comp,grad))
                    {
                        double locDiffusionCoefficient = DP.DiffusionCoefficients(i,j,k,{n,comp,grad});
                        double locMoleFraction = Cx.MoleFractions(i,j,k,{n,grad});

                        for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                        {
                            int di = ds->di;
                            int dj = ds->dj;
                            int dk = ds->dk;

                            double xyzPhaseFraction = Phi.Fractions(i+di,j+dj,k+dk,{n});
                            if (xyzPhaseFraction != 0.0)
                            {
                                double xyzDiffusionCoefficient = DP.DiffusionCoefficients(i+di,j+dj,k+dk,{n,comp,grad});
                                double xyzMoleFraction = Cx.MoleFractions(i+di,j+dj,k+dk,{n,grad});

                                // geometric mean between adjacent and local phase fraction
                                double meanPhi = 0.5*(locPhaseFraction + xyzPhaseFraction);

                                // arithmetic mean between adjacent and local diffusion coefficient
                                double meanDC = 0.5*(locDiffusionCoefficient + xyzDiffusionCoefficient);

                                // absolute values for diffusion coefficients based on the selected mode
                                if(positive_diffusivity_switch(Mode,comp,grad))
                                {
                                    meanDC = fabs(meanDC);
                                }

                                // difference between local and adjacent composition
                                double diffCx = locMoleFraction - xyzMoleFraction;

                                localDelta += -ds->weight*meanPhi*meanDC*diffCx;
                            }
                        }
                    }
                    if (localDelta < 0.0)
                    {
                        Cx.MoleFractionsDotOut(i, j, k, {n,comp}) += localDelta;
                    }
                    if (localDelta > 0.0)
                    {
                        Cx.MoleFractionsDotIn(i, j, k, {n,comp}) += localDelta;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ChemicalDiffusion_Impl::CalculateDiffusionIncrementsSF(
    PhaseField& Phi,
    Composition& Cx,
    Temperature& Tx,
    DiffusionProperties& DP,
    double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractions,0,)
    {
        // ---- compute phase sum once ----
        double phiSum = 0.0;
        for (size_t n = 0; n < Cx.Nphases; ++n)
            phiSum += Phi.Fractions(i,j,k,{n});

        if (phiSum < 1e-12) continue;

        for(size_t n = 0; n < Cx.Nphases; ++n)
        {
            double locPhaseFraction = Phi.Fractions(i,j,k,{n});
            if (locPhaseFraction == 0.0) continue;

            for (size_t comp = 0; comp < Cx.Ncomp; ++comp)
            if (!Cx.Phase[n].Component[comp].isReference &&
                !Cx.Phase[n].Component[comp].isFastDiffusor)
            {
                double localDelta = 0.0;

                for (size_t grad = 0; grad < Cx.Ncomp; ++grad)
                if (!Cx.Phase[n].Component[grad].isReference &&
                    !Cx.Phase[n].Component[grad].isFastDiffusor &&
                    model_component_switch(Model,comp,grad))
                {
                    double locDC = DP.DiffusionCoefficients(i,j,k,{n,comp,grad});
                    double locCx = Cx.MoleFractions(i,j,k,{n,grad});

                    for (auto ds = DStencil.begin(); ds != DStencil.end(); ++ds)
                    {
                        int di = ds->di, dj = ds->dj, dk = ds->dk;

                        double phi_nb = Phi.Fractions(i+di,j+dj,k+dk,{n});
                        if (phi_nb == 0.0) continue;

                        double nbDC = DP.DiffusionCoefficients(i+di,j+dj,k+dk,{n,comp,grad});
                        double nbCx = Cx.MoleFractions(i+di,j+dj,k+dk,{n,grad});

                        double meanPhi = 0.5*(locPhaseFraction + phi_nb) / phiSum;
                        double meanDC  = 0.5*(locDC + nbDC);

                        if (positive_diffusivity_switch(Mode,comp,grad))
                            meanDC = fabs(meanDC);

                        localDelta += -ds->weight * meanPhi * meanDC * (locCx - nbCx);
                    }
                }

                if (localDelta < 0.0)
                    Cx.MoleFractionsDotOut(i,j,k,{n,comp}) += localDelta;
                else
                    Cx.MoleFractionsDotIn (i,j,k,{n,comp}) += localDelta;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}



void ChemicalDiffusion_Impl::CalculateChemicalPotentialContribution(
                                              PhaseField& Phi,
                                              Composition& Cx,
                                              Temperature& Tx,
                                              DiffusionProperties& DP,
                                              ThermodynamicPropertiesEQP& TP,
                                              double dt)
{
    double R = PhysicalConstants::R;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractions,0,)
    {
        for(size_t n = 0; n < Cx.Nphases; n++)
        {
            double locPhaseFraction = Phi.Fractions(i,j,k,{n});
            if (locPhaseFraction != 0.0)
            {
                for (size_t comp = 0; comp < Cx.Ncomp; comp++)
                if (!Cx.Phase[n].Component[comp].isReference and
                    !Cx.Phase[n].Component[comp].isFastDiffusor)
                {
                    double localDelta = 0.0;

                    double locDiffusionCoefficient = DP.DiffusionCoefficients(i,j,k,{n,comp,comp});

                    for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                    {
                        int di = ds->di;
                        int dj = ds->dj;
                        int dk = ds->dk;

                        double xyzPhaseFraction = Phi.Fractions(i+di,j+dj,k+dk,{n});
                        if (xyzPhaseFraction != 0.0)
                        {
                            double xyzDiffusionCoefficient = DP.DiffusionCoefficients(i+di,j+dj,k+dk,{n,comp,comp});

                            // geometric mean between adjacent and local phase fraction
                            double meanPhi = 0.5*(locPhaseFraction + xyzPhaseFraction);

                            double invRT    = Cx.TotalMolarVolume/(R*Tx(i,j,k));
                            double invRTxyz = Cx.TotalMolarVolume/(R*Tx(i+di,j+dj,k+dk));

                            // arithmetic mean between adjacent and local chemical mobilities
                            double meanDC = 0.5*(invRT*locDiffusionCoefficient
                                            + invRTxyz*xyzDiffusionCoefficient);

                            // absolute values for diffusion coefficients based on the selected mode
                            if(positive_diffusivity_switch(Mode,comp,comp))
                            {
                                meanDC = fabs(meanDC);
                            }
                            // difference between adjacent and local chemical potential
                            double diffMu = TP.dMu(i,j,k,{n,comp})
                                           -TP.dMu(i+di,j+dj,k+dk,{n,comp});

                            localDelta += -ds->weight*meanPhi*meanDC*diffMu;
                        }
                    }
                    if (localDelta < 0.0)
                    {
                        Cx.MoleFractionsDotOut(i, j, k, {n,comp}) += localDelta;
                    }
                    if (localDelta > 0.0)
                    {
                        Cx.MoleFractionsDotIn(i, j, k, {n,comp}) += localDelta;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ChemicalDiffusion_Impl::Solve1Dextension(Composition& Cx,
                                              Temperature& Tx,
                                              DiffusionProperties& DP,
                                              Composition1Dextension& CxExt,
                                              BoundaryConditionTypes& extBC,
                                              double dt)
{
    std::pair<int, int> Xbounds = three_state_bounds_selector(CxExt.Direction[0], 0, Cx.Grid.Nx);
    std::pair<int, int> Ybounds = three_state_bounds_selector(CxExt.Direction[1], 0, Cx.Grid.Ny);
    std::pair<int, int> Zbounds = three_state_bounds_selector(CxExt.Direction[2], 0, Cx.Grid.Nz);

    Tensor<double,2> boundaryDC({Cx.Ncomp,Cx.Ncomp});
    Tensor<double,1> boundaryValue({Cx.Ncomp});
    double area = 0.0;

    // calculate boundary values of composition and diffusion coefficients
    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        for (size_t comp = 0; comp < Cx.Ncomp; comp++)
        {
            boundaryValue({comp}) += Cx.MoleFractions(i,j,k,{CxExt.PhaseIndex,comp});

            for (size_t grad = 0; grad < Cx.Ncomp; grad++)
            {
                boundaryDC({comp,grad}) += fabs(DP.DiffusionCoefficients(i,j,k,{CxExt.PhaseIndex,comp,grad}));
            }
        }
        area++;
    }

#ifdef MPI_PARALLEL
    if(CxExt.Direction[0] == 0)
    {
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, boundaryValue.data(), Cx.Ncomp, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, boundaryDC.data(), Cx.Ncomp*Cx.Ncomp, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &area, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    }
#endif

    // set boundary value of composition and diffusion coefficients in the extension
    CxExt.Data(-1) = boundaryValue/area;
    boundaryDC /= area;
    double dx_2 = 1.0/pow(Cx.Grid.dx,2);

    // calculate diffusion increments
    for(long int x = 0; x < (long int)CxExt.size(); x++)
    {
        CxExt.DataDot(x).set_to_zero();

        for (size_t comp = 0; comp < Cx.Ncomp; comp++)
        if (!Cx.Phase[CxExt.PhaseIndex].Component[comp].isReference)
        {
            double localDelta = 0.0;

            for (size_t grad = 0; grad < Cx.Ncomp; grad++)
            if(!Cx.Phase[CxExt.PhaseIndex].Component[grad].isReference and
               model_component_switch(Model,comp,grad))
            {
                localDelta += dx_2*boundaryDC({comp,grad})*(CxExt.Data(x-1,{grad}) - 2.0*CxExt.Data(x,{grad}) + CxExt.Data(x+1,{grad}));
            }
            CxExt.DataDot(x,{comp}) += localDelta;
        }
    }

    // apply diffusion increments
    for(long int x = 0; x < (long int)CxExt.size(); x++)
    {
        for (size_t comp = 0; comp < Cx.Ncomp; comp++)
        if (!Cx.Phase[CxExt.PhaseIndex].Component[comp].isReference)
        {
            CxExt.Data(x,{comp}) += CxExt.DataDot(x,{comp})*dt;
        }
    }

    // calculate mole fraction of reference element
    for(long int x = 0; x < (long int)CxExt.size(); x++)
    {
        double ref_MF = 1.0;
        for (size_t comp = 0; comp < Cx.Ncomp; comp++)
        if (!Cx.Phase[CxExt.PhaseIndex].Component[comp].isReference)
        {
            ref_MF -= CxExt.Data(x)({comp});
        }

        for (size_t comp = 0; comp < Cx.Ncomp; comp++)
        if (Cx.Phase[CxExt.PhaseIndex].Component[comp].isReference)
        {
            CxExt.Data(x,{comp}) = ref_MF;
        }
    }

    // set extension boundary conditions
    CxExt.setBC(extBC);

    // set boundary values in the simulation domain
    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    for (size_t comp = 0; comp < Cx.Ncomp; comp++)
    {
        Cx.MoleFractions(i+CxExt.Direction[0],j+CxExt.Direction[1],k+CxExt.Direction[2],{CxExt.PhaseIndex,comp}) = CxExt.Data(0)({comp});
    }
}

void ChemicalDiffusion_Impl::CalculateAntitrappingIncrements(PhaseField& Phi,
          Composition& Cx, Temperature& Tx,
          DiffusionProperties& DP, double dt)

{
    /** This function calculates the Anti-Trapping flux as described in
    10.1088/0965-0393/17/7/073001 equation 47 and 10.1016/j.actamat.2013.03.042
    equation 6. */

    double koef = Phi.Grid.Eta/Pi;
    double dx = Phi.Grid.dx;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    if(Phi.Fields(i,j,k).interface())
    {
        for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
        {
            int di = ds->di;
            int dj = ds->dj;
            int dk = ds->dk;

            dVector3 stencil_direction(std::array<int,3>{di,dj,dk});
            stencil_direction.normalize();

            NodePF locPhaseFields = Phi.Fields(i,j,k);
            locPhaseFields.add_existing_all(Phi.Fields(i+di, j+dj, k+dk));
            locPhaseFields *= 0.5;

            if(locPhaseFields.size() > 1)
            for(auto alpha = locPhaseFields.cbegin();
                     alpha != locPhaseFields.cend() - 1; ++alpha)
            {
                size_t gIndexA = alpha->index;
                size_t pIndexA = Phi.FieldsProperties[gIndexA].Phase;

                double scaleA = Phi.FieldsProperties[gIndexA].VolumeRatio;

                for(auto beta = alpha + 1;
                         beta != locPhaseFields.cend(); ++beta)
                {
                    size_t gIndexB = beta->index;
                    size_t pIndexB = Phi.FieldsProperties[gIndexB].Phase;

                    if(pIndexA != pIndexB)
                    {
                        double scaleB = Phi.FieldsProperties[gIndexB].VolumeRatio;

                        const dVector3& Normal = Phi.Normal(alpha,beta);
                        double Projection = (Normal*stencil_direction);

                        double locPsi = 0.5
                        *(Phi.FieldsDot(i   ,j   ,k   ).get_asym1(gIndexA,gIndexB)
                         +Phi.FieldsDot(i+di,j+dj,k+dk).get_asym1(gIndexA,gIndexB));

                        for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                        if (!(Cx.Phase[pIndexA].Component[comp].isReference) and
                            !(Cx.Phase[pIndexB].Component[comp].isReference))
                        {
                            double locDCalpha
                                 = 0.5*(DP.DiffusionCoefficients(i,j,k,{pIndexA,comp,comp})
                                       +DP.DiffusionCoefficients(i+di,j+dj,k+dk,{pIndexA,comp,comp}));

                            double locDCbeta
                                 = 0.5*(DP.DiffusionCoefficients(i,j,k,{pIndexB,comp,comp})
                                       +DP.DiffusionCoefficients(i+di,j+dj,k+dk,{pIndexB,comp,comp}));

                            // absolute values for diffusion coefficients based on the selected mode
                            if(positive_diffusivity_switch(Mode,comp,comp))
                            {
                                locDCalpha = fabs(locDCalpha);
                                locDCbeta  = fabs(locDCbeta);
                            }

                            double localDelta = -dx*ds->weight*
                                        sqrt(alpha->value * beta->value)
                                       *((Cx.MoleFractions(i,j,k,{pIndexA,comp})
                                        + Cx.MoleFractions(i+di,j+dj,k+dk,{pIndexA,comp}))
                                        -(Cx.MoleFractions(i,j,k,{pIndexB,comp})
                                        + Cx.MoleFractions(i+di,j+dj,k+dk,{pIndexB,comp})))
                                         *locPsi*Projection*koef*scaleA*scaleB
                                        *(locDCalpha - locDCbeta)/(locDCalpha + locDCbeta);

                            if (localDelta < -DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotOut(i,j,k,{pIndexA,comp}) += localDelta;
                            }
                            if (localDelta > DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotIn(i,j,k,{pIndexA,comp}) += localDelta;
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

bool ChemicalDiffusion_Impl::CalculateLimits(PhaseField& Phi,
                                    Composition& Cx, double dt)
{
    bool LimitingNeeded = false;

    // Clear limits
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,Cx.MoleFractionsTotal.Bcells(),)
    {
        for(size_t alpha = 0; alpha < Cx.Nphases; ++alpha)
        for(size_t comp = 0; comp < Cx.Ncomp; comp++)
        {
            Cx.NormIn(i,j,k,{alpha,comp})  = 1.0;
            Cx.NormOut(i,j,k,{alpha,comp}) = 1.0;
        }
        Cx.Limiting(i,j,k) = false;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    // Calculate limits
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,reduction(||:LimitingNeeded))
    {
        bool out_of_limits = false;

        for(size_t n = 0; n < Cx.Nphases; n++)
        if (Phi.Fractions(i,j,k,{n}) > 0.0)
        for(size_t comp = 0; comp < Cx.Ncomp; comp++)
        {
            double dPhaseIn  = Cx.MoleFractionsDotIn(i,j,k,{n,comp})*dt
                               /Phi.Fractions(i,j,k,{n});
            double dPhaseOut = Cx.MoleFractionsDotOut(i,j,k,{n,comp})*dt
                               /Phi.Fractions(i,j,k,{n});

            if(dPhaseIn  <=  DBL_EPSILON) dPhaseIn  = 0.0;
            if(dPhaseOut >= -DBL_EPSILON) dPhaseOut = 0.0;

            double PhaseOld = Cx.MoleFractions(i,j,k,{n,comp});

            if((PhaseOld + dPhaseIn > Cx.Phase[n].Component[comp].Max) and
               (dPhaseIn != 0.0))
            {
                Cx.NormIn(i,j,k,{n,comp})
                = (Cx.Phase[n].Component[comp].Max - PhaseOld)/dPhaseIn;

                out_of_limits = true;
            }

            if((PhaseOld + dPhaseOut < Cx.Phase[n].Component[comp].Min) and
               (dPhaseOut != 0.0))
            {
                Cx.NormOut(i,j,k,{n,comp})
                = (Cx.Phase[n].Component[comp].Min - PhaseOld)/dPhaseOut;

                out_of_limits = true;
            }
        }

        if(out_of_limits)
        {
            LimitingNeeded = true;

            for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
            {
                int di = ds->di;
                int dj = ds->dj;
                int dk = ds->dk;

                Cx.Limiting(i+di,j+dj,k+dk) = true;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE,&LimitingNeeded, 1, OP_MPI_CXX_BOOL, OP_MPI_LOR, OP_MPI_COMM_WORLD);
    #endif
    return LimitingNeeded;
}

void ChemicalDiffusion_Impl::LimitDiffusionIncrements(
                   PhaseField& Phi, Composition& Cx, Temperature& Tx,
                   DiffusionProperties& DP, double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractions,0,)
    if(Cx.Limiting(i,j,k))
    {
        for(size_t n = 0; n < Cx.Nphases; n++)
        if (Phi.Fractions(i,j,k,{n}) > 0.0)
        {
            for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
            {
                int di = ds->di;
                int dj = ds->dj;
                int dk = ds->dk;

                if (Phi.Fractions(i+di,j+dj,k+dk,{n}) > 0.0)
                {
                    // arithmetic mean between adjacent and local phase fraction
                    double meanPhi = 0.5*(Phi.Fractions(i,j,k,{n})
                                        + Phi.Fractions(i+di,j+dj,k+dk,{n}));

                    for (size_t comp = 0; comp < Cx.Ncomp; comp++)
                    if (!Cx.Phase[n].Component[comp].isReference and
                        !Cx.Phase[n].Component[comp].isFastDiffusor)
                    {
                        double localDelta = 0.0;
                        for (size_t grad = 0; grad < Cx.Ncomp; grad++)
                        if(!Cx.Phase[n].Component[grad].isReference and
                           !Cx.Phase[n].Component[grad].isFastDiffusor and
                           model_component_switch(Model,comp,grad))
                        {
                            // arithmetic mean between adjacent and local diffusion coefficient
                            double meanDC = 0.5*(DP.DiffusionCoefficients(i,j,k,{n,comp,grad})
                                               + DP.DiffusionCoefficients(i+di,j+dj,k+dk,{n,comp,grad}));

                            // absolute values for diffusion coefficients based on the selected mode
                            if(positive_diffusivity_switch(Mode,comp,grad))
                            {
                                meanDC = fabs(meanDC);
                            }
                            // difference between adjacent and local composition
                            double diffCx = Cx.MoleFractions(i,j,k,{n,grad})
                                          - Cx.MoleFractions(i+di,j+dj,k+dk,{n,grad});

                            localDelta += -ds->weight*meanPhi*meanDC*diffCx;
                        }

                        if (Cx.NormIn(i,j,k,{n,comp}) <= Cx.NormOut(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormIn(i, j, k, {n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement > DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotIn(i,j,k,{n,comp}) -= localIncrement;
                            }
                        }
                        if (Cx.NormIn(i,j,k,{n,comp}) > Cx.NormOut(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormOut(i+di,j+dj,k+dk,{n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement > DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotIn(i,j,k,{n,comp}) -= localIncrement;
                            }
                        }

                        if (Cx.NormOut(i,j,k,{n,comp}) < Cx.NormIn(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormOut(i, j, k, {n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement < -DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotOut(i,j, k, {n,comp}) -= localIncrement;
                            }
                        }
                        if (Cx.NormOut(i,j,k,{n,comp}) >= Cx.NormIn(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormIn(i+di, j+dj, k+dk, {n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement < -DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotOut(i,j,k,{n,comp}) -= localIncrement;
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ChemicalDiffusion_Impl::LimitChemicalPotentialContribution(
                   PhaseField& Phi, Composition& Cx, Temperature& Tx,
                   DiffusionProperties& DP, ThermodynamicPropertiesEQP& TP,
                   double dt)
{
    double R = PhysicalConstants::R;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractions,0,)
    if(Cx.Limiting(i,j,k))
    {
        for(size_t n = 0; n < Cx.Nphases; n++)
        if (Phi.Fractions(i,j,k,{n}) > 0.0)
        {
            double invRT = Cx.TotalMolarVolume/(R*Tx(i,j,k));

            for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
            {
                int di = ds->di;
                int dj = ds->dj;
                int dk = ds->dk;

                if (Phi.Fractions(i+di,j+dj,k+dk,{n}) > 0.0)
                {
                    // arithmetic mean between adjacent and local phase fraction
                    double meanPhi = 0.5*(Phi.Fractions(i,j,k,{n})
                                        + Phi.Fractions(i+di,j+dj,k+dk,{n}));

                    for (size_t comp = 0; comp < Cx.Ncomp; comp++)
                    if (!Cx.Phase[n].Component[comp].isReference and
                        !Cx.Phase[n].Component[comp].isFastDiffusor)
                    {
                        double invRTxyz  = Cx.TotalMolarVolume/(R*Tx(i+di,j+dj,k+dk));

                        // arithmetic mean between adjacent and local diffusion coefficient
                        double meanDC = 0.5*(invRT*DP.DiffusionCoefficients(i,j,k,{n,comp,comp})
                                        + invRTxyz*DP.DiffusionCoefficients(i+di,j+dj,k+dk,{n,comp,comp}));

                        // absolute values for diffusion coefficients based on selected mode
                        if(positive_diffusivity_switch(Mode,comp,comp))
                        {
                            meanDC = fabs(meanDC);
                        }

                        double diffMu = TP.dMu(i,j,k,{n,comp})
                                       -TP.dMu(i+di,j+dj,k+dk,{n,comp});

                        double localDelta = -ds->weight*meanPhi*meanDC*diffMu;


                        if (Cx.NormIn(i,j,k,{n,comp}) <= Cx.NormOut(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormIn(i, j, k, {n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement > DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotIn(i,j,k,{n,comp}) -= localIncrement;
                            }
                        }
                        if (Cx.NormIn(i,j,k,{n,comp}) > Cx.NormOut(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormOut(i+di,j+dj,k+dk,{n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement > DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotIn(i,j,k,{n,comp}) -= localIncrement;
                            }
                        }

                        if (Cx.NormOut(i,j,k,{n,comp}) < Cx.NormIn(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormOut(i, j, k, {n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement < -DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotOut(i,j, k, {n,comp}) -= localIncrement;
                            }
                        }
                        if (Cx.NormOut(i,j,k,{n,comp}) >= Cx.NormIn(i+di,j+dj,k+dk,{n,comp}))
                        {
                            double localNorm = Cx.NormIn(i+di, j+dj, k+dk, {n,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement < -DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotOut(i,j,k,{n,comp}) -= localIncrement;
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ChemicalDiffusion_Impl::LimitAntitrappingIncrements(
                   PhaseField& Phi, Composition& Cx, Temperature& Tx,
                   DiffusionProperties& DP,
                   double dt)
{
    double koef = Phi.Grid.Eta/Pi;
    double dx = Phi.Grid.dx;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    if(Cx.Limiting(i, j, k))
    if(Phi.Fields(i,j,k).interface())
    for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
    {
        int di = ds->di;
        int dj = ds->dj;
        int dk = ds->dk;
        dVector3 stencil_direction(std::array<int,3>{di,dj,dk});
        stencil_direction.normalize();

        NodePF locPhaseFields = Phi.Fields(i,j,k);
        locPhaseFields.add_existing_all(Phi.Fields(i+di, j+dj, k+dk));
        locPhaseFields *= 0.5;

        for(auto alpha = locPhaseFields.cbegin();
                 alpha != locPhaseFields.cend() - 1; ++alpha)
        {
            size_t gIndexA = alpha->index;
            size_t pIndexA = Phi.FieldsProperties[gIndexA].Phase;

            double scaleA = Phi.FieldsProperties[gIndexA].VolumeRatio;

            for(auto beta = alpha + 1;
                     beta != locPhaseFields.cend(); ++beta)
            {
                size_t gIndexB = beta->index;
                size_t pIndexB = Phi.FieldsProperties[gIndexB].Phase;

                if(pIndexA != pIndexB)
                {
                    double scaleB = Phi.FieldsProperties[gIndexB].VolumeRatio;

                    dVector3 Normal = Phi.Normal(alpha, beta);
                    double Projection = (Normal*stencil_direction);

                    double locPsi = 0.5
                    *(Phi.FieldsDot(i   ,j   ,k   ).get_asym1(gIndexA,gIndexB)
                     +Phi.FieldsDot(i+di,j+dj,k+dk).get_asym1(gIndexA,gIndexB));

                    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                    if (!(Cx.Phase[pIndexA].Component[comp].isReference) and
                        !(Cx.Phase[pIndexB].Component[comp].isReference))
                    {
                        double locDCalpha
                             = 0.5*(DP.DiffusionCoefficients(i,j,k,{pIndexA,comp,comp})
                                   +DP.DiffusionCoefficients(i+di,j+dj,k+dk,{pIndexA,comp,comp}));

                        double locDCbeta
                             = 0.5*(DP.DiffusionCoefficients(i,j,k,{pIndexB,comp,comp})
                                   +DP.DiffusionCoefficients(i+di,j+dj,k+dk,{pIndexB,comp,comp}));

                        if(positive_diffusivity_switch(Mode,comp,comp))
                        {
                            locDCalpha = fabs(locDCalpha);
                            locDCbeta  = fabs(locDCbeta);
                        }

                        double localDelta = dx*(-ds->weight)*
                                    sqrt(alpha->value * beta->value)
                                   *((Cx.MoleFractions(i,j,k,{pIndexA,comp})
                                    + Cx.MoleFractions(i+di,j+dj,k+dk,{pIndexA,comp}))
                                    -(Cx.MoleFractions(i,j,k,{pIndexB,comp})
                                    + Cx.MoleFractions(i+di,j+dj,k+dk,{pIndexB,comp})))
                                     *locPsi*Projection*koef*scaleA*scaleB
                                    *(locDCalpha - locDCbeta)/(locDCalpha + locDCbeta);

                        if (Cx.NormIn (i   , j   , k   , {pIndexA,comp}) <=
                            Cx.NormOut(i+di, j+dj, k+dk, {pIndexA,comp}) and
                            localDelta > DBL_EPSILON)
                        {
                            double localNorm = Cx.NormIn(i,j,k,{pIndexA,comp});
                            double localIncrement = localDelta*(1.0 - localNorm);

                            if(localIncrement > DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotIn (i   , j   , k   , {pIndexA,comp})
                                -= localIncrement;

                                Cx.MoleFractionsDotOut(i+di, j+dj, k+dk, {pIndexA,comp})
                                += localIncrement;
                            }
                        }
                        if (Cx.NormOut(i   , j   , k   , {pIndexA,comp}) <
                            Cx.NormIn (i+di, j+dj, k+dk, {pIndexA,comp}) and
                            localDelta < -DBL_EPSILON)
                        {
                            double localNorm = Cx.NormOut(i,j,k,{pIndexA, comp});
                            double localIncrement = localDelta*(1.0 - localNorm);

                            if(localIncrement < -DBL_EPSILON)
                            {
                                Cx.MoleFractionsDotOut(i   , j   , k   , {pIndexA, comp})
                                -= localIncrement;

                                Cx.MoleFractionsDotIn (i+di, j+dj, k+dk, {pIndexA, comp})
                                += localIncrement;
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ChemicalDiffusion_Impl::ApplyIncrements(PhaseField& Phi,
                                    Composition& Cx, double dt)
{
    /** This function applies the previously calculated phase composition fluxes
    to the total composition.*/

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        for(size_t n = 0; n < Cx.Nphases; n++)
        if (Phi.Fractions(i,j,k,{n}) != 0.0)
        {
            double refdot = 0.0;

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            if(!(Cx.Phase[n].Component[comp].isReference))
            {
                double dPhase = (Cx.MoleFractionsDotIn(i,j,k,{n,comp}) +
                                 Cx.MoleFractionsDotOut(i,j,k,{n,comp}))*dt;

                Cx.MoleFractionsTotal(i,j,k, {comp}) += dPhase;
                refdot -= dPhase;
            }

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            if(Cx.Phase[n].Component[comp].isReference)
            {
                Cx.MoleFractionsTotal(i,j,k, {comp}) += refdot;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ChemicalDiffusion_Impl::ApplyIncrementsSF(
    PhaseField& Phi,
    Composition& Cx,
    double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        for (size_t n = 0; n < Cx.Nphases; ++n)
        if (Phi.Fractions(i,j,k,{n}) != 0.0)
        {
            double refDot = 0.0;

            for (size_t comp = 0; comp < Cx.Ncomp; ++comp)
            if (!Cx.Phase[n].Component[comp].isReference)
            {
                double dot =
                    Cx.MoleFractionsDotIn (i,j,k,{n,comp}) +
                    Cx.MoleFractionsDotOut(i,j,k,{n,comp});

                Cx.MoleFractionsTotal   (i,j,k,{comp}) += dot * dt;
                Cx.MoleFractionsTotalDot(i,j,k,{comp}) += dot;

                refDot -= dot;
            }

            for (size_t comp = 0; comp < Cx.Ncomp; ++comp)
            if (Cx.Phase[n].Component[comp].isReference)
            {
                Cx.MoleFractionsTotal   (i,j,k,{comp}) += refDot * dt;
                Cx.MoleFractionsTotalDot(i,j,k,{comp}) += refDot;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


}// namespace openphase
