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

 *   File created :   2015
 *   Main contributors :   Philipp Engels; Marvin Tegeler; Raphael Schiedung
 *
 */

#include "AdvectionHR.h"
#include "FileInterface.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "PhaseField.h"
#include <set>

namespace openphase
{

AdvectionHR::AdvectionHR(Settings& locSettings, std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void AdvectionHR::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "AdvectionHR";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Scheme = AdvectionSchemes::Upwind;

    if(locSettings.Grid.Bcells < 2)
    {
        std::string message  = "AdvectionHD requires at least 2 boundary cells for correct operation!\n";
                    message += "The actually selected number of boundary cells is ";
                    message += std::to_string(locSettings.Grid.Bcells);
                    message += ".\n";
                    message += "Adjust parameter $Bcells in @Settings section of the project input file accordingly!";
        ConsoleOutput::WriteExit(message, thisclassname, "Initialize()");
        OP_Exit(EXIT_FAILURE);
    }
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void AdvectionHR::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::fstream inpF(InputFileName.c_str(), std::ios::in | std::ios_base::binary);
    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    ConsoleOutput::WriteLineInsert("AdvectionHR input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    ReadInput(inp);

    ConsoleOutput::WriteLine();
}

void AdvectionHR::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);
    std::string schemeString = FileInterface::ReadParameterK(inp, moduleLocation, "scheme", false, "UPWIND");

    if (schemeString == "MINMOD")
    {
        Scheme = AdvectionSchemes::Minmod;
    }
    else if (schemeString == "MC" or schemeString == "MONOTONISEDCENTRAL")
    {
        Scheme = AdvectionSchemes::MonotonizedCentral;
    }
    else if (schemeString == "SUPERBEE")
    {
        Scheme = AdvectionSchemes::Superbee;
    }
    /*else if (schemeString == "LAXWENDROFF")
    {
        Scheme = AdvectionSchemes::LaxWendroff;
    }*/
    else if (schemeString == "UPWIND")
    {
        Scheme = AdvectionSchemes::Upwind;
    }
    else
    {
        ConsoleOutput::WriteWarning("Wrong advection scheme specified!\nThe default \"Upwind\" scheme is used!", thisclassname, "ReadInput()");
    }
    switch(Scheme)
    {
        case AdvectionSchemes::Minmod:
        {
            Limiter = &Slope_Limiter_Minmod;
            break;
        }
        case AdvectionSchemes::MonotonizedCentral:
        {
            Limiter = &Slope_Limiter_MC;
            break;
        }
        case AdvectionSchemes::Superbee:
        {
            Limiter = &Slope_Limiter_Superbee;
            break;
        }
        /*case AdvectionSchemes::LaxWendroff:
        {
            Limiter = &Slope_Limiter_LaxWendroff;
            break;
        }*/
        case AdvectionSchemes::Upwind:
        default:
        {
            Limiter = &Slope_Limiter_Upwind;
            break;
        }
    }
}

void AdvectionHR::CalculateLocalAdvectionPhaseField(
        const int i, const int j, const int k,
        const int ii, const int jj, const int kk,
        PhaseField& Phase,
        const BoundaryConditions& BC,
        const Velocities& Vel,
        int& CFL,
        const int direction,
        const double dt,
        double (*Limiter)(const double,const double),
        size_t FluidIdx)
{
    if (Phase.Fields(i,j,k).wide_interface())
    {
        const auto& q   = Phase.Fields(i     ,j     ,k     );
        const auto& qp  = Phase.Fields(i+  ii,j+  jj,k+  kk);
        const auto& qm  = Phase.Fields(i-  ii,j-  jj,k-  kk);
        const auto& qpp = Phase.Fields(i+2*ii,j+2*jj,k+2*kk);
        const auto& qmm = Phase.Fields(i-2*ii,j-2*jj,k-2*kk);

        // Determine indices of close by grains
        std::set<size_t> GrainIdxs;
        for (auto alpha = q  .cbegin(); alpha != q  .cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qp .cbegin(); alpha != qp .cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qm .cbegin(); alpha != qm .cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qpp.cbegin(); alpha != qpp.cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qmm.cbegin(); alpha != qmm.cend(); alpha++) GrainIdxs.insert(alpha->index);

        bool SolidInterface = false;
        for (auto idx = GrainIdxs.cbegin(); idx != GrainIdxs.cend(); idx++)
        {
            const Grain& grain = Phase.FieldsProperties[*idx];
            if (grain.State == AggregateStates::Solid)
            {
                SolidInterface = true;
                const double v  = InteractionSolidFluid::LocalGrainVelocity(i,j,k,grain,BC, Phase)[direction];
                const double dphi_dt = AdvectionKernel<false>(CFL,v,v,v,q.get_value(*idx),qp.get_value(*idx),qm.get_value(*idx),qpp.get_value(*idx),qmm.get_value(*idx),dt,Vel.Grid.dx,Limiter);
                Phase.FieldsAdvectionDot(i,j,k).add_value(    *idx, dphi_dt);
                Phase.FieldsAdvectionDot(i,j,k).add_value(FluidIdx,-dphi_dt);
            }
        }

        if (!SolidInterface)
        {
            const double v  = Vel.Average(i   ,j    ,k  )[direction];
            const double vp = Vel.Average(i+ii,j+jj,k+kk)[direction];
            const double vm = Vel.Average(i-ii,j-jj,k-kk)[direction];

            for (auto idx = GrainIdxs.cbegin(); idx != GrainIdxs.cend(); idx++)
            {
                const double dphi_dt = AdvectionKernel<true>(CFL,v,vp,vm,q.get_value(*idx),qp.get_value(*idx),qm.get_value(*idx),qpp.get_value(*idx),qmm.get_value(*idx),dt,Vel.Grid.dx,Limiter);
                Phase.FieldsAdvectionDot(i,j,k).set_value(*idx, dphi_dt);
            }
        }
    }
}

void AdvectionHR::CalculateLocalAdvectionPhaseFieldALE(
        const int i, const int j, const int k,
        const int ii, const int jj, const int kk,
        PhaseField& Phase,
        const BoundaryConditions& BC,
        const Velocities& Vel,
        int& CFL,
        const int direction,
        const double dt,
        double (*Limiter)(const double,const double))
{
    if (Phase.Fields(i,j,k).wide_interface())
    {
        const auto& q   = Phase.Fields(i     ,j     ,k     );
        const auto& qp  = Phase.Fields(i+  ii,j+  jj,k+  kk);
        const auto& qm  = Phase.Fields(i-  ii,j-  jj,k-  kk);
        const auto& qpp = Phase.Fields(i+2*ii,j+2*jj,k+2*kk);
        const auto& qmm = Phase.Fields(i-2*ii,j-2*jj,k-2*kk);

        // Determine indices of close by grains
        std::set<size_t> GrainIdxs;
        for (auto alpha = q  .cbegin(); alpha != q  .cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qp .cbegin(); alpha != qp .cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qm .cbegin(); alpha != qm .cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qpp.cbegin(); alpha != qpp.cend(); alpha++) GrainIdxs.insert(alpha->index);
        for (auto alpha = qmm.cbegin(); alpha != qmm.cend(); alpha++) GrainIdxs.insert(alpha->index);

        const double v  = Vel.Average(i   ,j    ,k  )[direction];
        const double vp = Vel.Average(i+ii,j+jj,k+kk)[direction];
        const double vm = Vel.Average(i-ii,j-jj,k-kk)[direction];

        for (auto idx = GrainIdxs.cbegin(); idx != GrainIdxs.cend(); idx++)
        {
            const double dphi_dt = AdvectionKernel<false>(CFL,v,vp,vm,q.get_value(*idx),qp.get_value(*idx),qm.get_value(*idx),qpp.get_value(*idx),qmm.get_value(*idx),dt,Vel.Grid.dx,Limiter);
            Phase.FieldsAdvectionDot(i,j,k).set_value(*idx, dphi_dt);
        }
    }
}

void AdvectionHR::CalculateAdvectionPhaseField(
        PhaseField& Phase,
        const BoundaryConditions& BC,
        const Velocities& Vel,
        const int direction,
        double dt,
        double (*Limiter)(const double,const double),
        size_t FluidIdx,
        bool& warning)
{
    Phase.FieldsAdvectionBackup = Phase.Fields;

    const int ii = (direction == 0) ? 1 : 0;
    const int jj = (direction == 1) ? 1 : 0;
    const int kk = (direction == 2) ? 1 : 0;

    unsigned int iterations = 1;
    int CFL;
    do
    {
        CFL = 0;
        for (unsigned int it = 0; it < iterations; ++it)
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:CFL))
            {
                Phase.FieldsAdvectionDot(i,j,k).clear();
                CalculateLocalAdvectionPhaseField(i,j,k,ii,jj,kk,Phase,BC,Vel,CFL,direction,dt,Limiter, FluidIdx);
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:CFL))
            {
                for (auto alpha  = Phase.FieldsAdvectionDot(i,j,k).begin();
                          alpha != Phase.FieldsAdvectionDot(i,j,k).end(); alpha++)
                {
                    Phase.Fields(i,j,k).add_value(alpha->index, alpha->value*dt);
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            BC.Set(Phase.Fields,direction);
        }
        if (CFL != 0)
        {
            Phase.Fields = Phase.FieldsAdvectionBackup;

            dt *= 0.5;
            warning = true;
        }
        iterations *= 2;
    }
    while (CFL > 0);
}

void AdvectionHR::CalculateAdvectionPhaseFieldALE(
        PhaseField& Phase,
        const BoundaryConditions& BC,
        const Velocities& Vel,
        const int direction,
        double dt,
        double (*Limiter)(const double,const double),
        bool& warning)
{
    Phase.FieldsAdvectionBackup = Phase.Fields;

    const int ii = (direction == 0) ? 1 : 0;
    const int jj = (direction == 1) ? 1 : 0;
    const int kk = (direction == 2) ? 1 : 0;

    unsigned int iterations = 1;
    int CFL;
    do
    {
        CFL = 0;
        for (unsigned int it = 0; it < iterations; ++it)
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:CFL))
            {
                Phase.FieldsAdvectionDot(i,j,k).clear();
                CalculateLocalAdvectionPhaseFieldALE(i,j,k,ii,jj,kk,Phase,BC,Vel,CFL,direction,dt,Limiter);
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:CFL))
            {
                for (auto alpha  = Phase.FieldsAdvectionDot(i,j,k).begin();
                          alpha != Phase.FieldsAdvectionDot(i,j,k).end(); alpha++)
                {
                    Phase.Fields(i,j,k).add_value(alpha->index, alpha->value*dt);
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            BC.Set(Phase.Fields,direction);
        }
        if (CFL != 0)
        {
            Phase.Fields = Phase.FieldsAdvectionBackup;

            dt *= 0.5;
            warning = true;
        }
        iterations *= 2;
    }
    while (CFL > 0);
}

size_t AdvectionHR::GetFluidIndex(const PhaseField& Phase, long tStep)
{
    const auto methodname = "GetFluidIndex()";
    size_t num_gas_grain = 0;
    std::set<size_t> gas_phases;
    for (size_t idx = 0; idx != Phase.FieldsProperties.size(); idx++)
    {
        if (Phase.FieldsProperties[idx].State == AggregateStates::Gas)
        {
            gas_phases.insert(Phase.FieldsProperties[idx].Phase);
            num_gas_grain++;
        }
    }

    if (gas_phases.size() > 1)
    {
        ConsoleOutput::WriteExit("More than one gas phase found! Immiscible gas phases can not be considered!", thisclassname, methodname);
        OP_Exit(EXIT_FAILURE);
    }
    else if (num_gas_grain > 1)
    {
        ConsoleOutput::WriteExit("More than one gas grain found! Immiscible gas grain can not be considered!", thisclassname, methodname);
        OP_Exit(EXIT_FAILURE);
    }
    else if (gas_phases.size() == 1 and num_gas_grain == 1)
    {
        for (size_t idx = 0; idx != Phase.FieldsProperties.size(); idx++)
        {
            if (Phase.FieldsProperties[idx].State == AggregateStates::Gas) return idx;
        }
    }

    size_t num_liquid_grain = 0;
    std::set<size_t> liquid_phases;
    for (size_t idx = 0; idx != Phase.FieldsProperties.size(); idx++)
    {
        if (Phase.FieldsProperties[idx].State == AggregateStates::Liquid)
        {
            liquid_phases.insert(Phase.FieldsProperties[idx].Phase);
            num_liquid_grain++;
        }
    }

    if (liquid_phases.size() > 1)
    {
        ConsoleOutput::WriteExit("Multiple liquid phases found. Add vapor phase to calculate advection!", thisclassname, methodname);
        OP_Exit(EXIT_FAILURE);
    }
    else if (num_liquid_grain > 1)
    {
        ConsoleOutput::WriteExit("More than one liquid grain found! Add vapor phase to calculate advection!", thisclassname, methodname);
        OP_Exit(EXIT_FAILURE);
    }
    else if (liquid_phases.size() == 1 and num_liquid_grain == 1)
    {
        if (tStep == 0) ConsoleOutput::WriteWarning("Consider to add a vapor phase to calculate advection! To allow modeling of cavitation with phase-field.", thisclassname, methodname);
        for (size_t idx = 0; idx != Phase.FieldsProperties.size(); idx++)
        {
            if (Phase.FieldsProperties[idx].State == AggregateStates::Liquid) return idx;
        }
    }

    //NOTE that there must be a unique Fluid index for the current algorithm
    // to function. If you want to consider multiple fluid grains consider
    // to implement a nucleation algorithm which may nucleate gas grain in
    // the solid-solid interface due to advection.
    ConsoleOutput::WriteExit("No fluid phase found! Can not calculate phase-field advection.", thisclassname, methodname);
    OP_Exit(EXIT_FAILURE);
    return 0;
}

void AdvectionHR::AdvectPhaseField(
        PhaseField& Phase,
        const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dt,
        const int tStep,
        const bool finalize)
{
    assert(Phase.Fields.Bcells() >= 2 && "Number of Bcells for storage PhaseField::Fields needs to be 2 or higher.");

    size_t FluidIdx = GetFluidIndex(Phase, tStep);

    bool warning = false;
    if (tStep % 2 == 0)
    {
        if (Vel.Grid.dNx) CalculateAdvectionPhaseField(Phase, BC, Vel, 0, dt, Limiter, FluidIdx, warning);
        if (Vel.Grid.dNy) CalculateAdvectionPhaseField(Phase, BC, Vel, 1, dt, Limiter, FluidIdx, warning);
        if (Vel.Grid.dNz) CalculateAdvectionPhaseField(Phase, BC, Vel, 2, dt, Limiter, FluidIdx, warning);
    }
    else
    {
        if (Vel.Grid.dNz) CalculateAdvectionPhaseField(Phase, BC, Vel, 2, dt, Limiter, FluidIdx, warning);
        if (Vel.Grid.dNy) CalculateAdvectionPhaseField(Phase, BC, Vel, 1, dt, Limiter, FluidIdx, warning);
        if (Vel.Grid.dNx) CalculateAdvectionPhaseField(Phase, BC, Vel, 0, dt, Limiter, FluidIdx, warning);
    }

    if (warning)
    {
        ConsoleOutput::WriteWarning("Warning reduce time step, cfl > 0.5!", "AdvectionHR", "CalculateAdvection");
    }

    Phase.Finalize(BC,finalize);
}

void AdvectionHR::AdvectPhaseFieldALE(
        PhaseField& Phase,
        const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dt,
        const int tStep,
        const bool finalize)
{
    assert(Phase.Fields.Bcells() >= 2 && "Number of Bcells for storage PhaseField::Fields needs to be 2 or higher.");

    bool warning = false;
    if (tStep % 2 == 0)
    {
        if (Vel.Grid.dNx) CalculateAdvectionPhaseFieldALE(Phase, BC, Vel, 0, dt, Limiter, warning);
        if (Vel.Grid.dNy) CalculateAdvectionPhaseFieldALE(Phase, BC, Vel, 1, dt, Limiter, warning);
        if (Vel.Grid.dNz) CalculateAdvectionPhaseFieldALE(Phase, BC, Vel, 2, dt, Limiter, warning);
    }
    else
    {
        if (Vel.Grid.dNz) CalculateAdvectionPhaseFieldALE(Phase, BC, Vel, 2, dt, Limiter, warning);
        if (Vel.Grid.dNy) CalculateAdvectionPhaseFieldALE(Phase, BC, Vel, 1, dt, Limiter, warning);
        if (Vel.Grid.dNx) CalculateAdvectionPhaseFieldALE(Phase, BC, Vel, 0, dt, Limiter, warning);
    }

    if (warning)
    {
        ConsoleOutput::WriteWarning("Warning reduce time step, cfl > 0.5!", "AdvectionHR", "CalculateAdvection");
    }

    Phase.Finalize(BC,finalize);
}
}// namespace openphase
