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
 *  File created :   2014
 *  Main contributors :   Oleg Shchyglo; Matthias Stratmann
 *
 */

#include "Diffusion/DiffusionProperties.h"
#include "Settings.h"
#include "Includes.h"
#include "ConsoleOutput.h"
#include "FileInterface.h"
#include "PhysicalConstants.h"
#include "BoundaryConditions.h"
#include "PhaseField.h"
#include "Composition.h"
#include "Temperature.h"
#include "VTK.h"

namespace openphase
{

using namespace std;

void DiffusionProperties::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "DiffusionProperties";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;
    Ncomp = locSettings.Ncomp;
    PhaseNames = locSettings.PhaseNames;
    ElementNames = locSettings.ElementNames;

    MaxDiffusionCoefficient  = 0.0;
    MaxMoleFractionDeviation = 1.0e-3;
    MaxTemperatureDeviation  = 1.0;

    ExtrapolationMode = ExtrapolationModes::None;

    size_t Bcells = Grid.Bcells;

    GlobalExtrapolationData.Allocate(std::array<size_t,1>{Nphases});
    for(size_t m = 0; m < Nphases; m++)
    {
        GlobalExtrapolationData({m}).Allocate(Ncomp);
    }
    LocalExtrapolationData.Allocate(Grid,{Nphases},Bcells);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,LocalExtrapolationData,LocalExtrapolationData.Bcells(),)
    {
        for(size_t m = 0; m < Nphases; m++)
        {
            LocalExtrapolationData(i,j,k,{m}).Allocate(Ncomp);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    DiffusionCoefficients.Allocate(Grid,{Nphases,Ncomp,Ncomp},Bcells);

    ScalingFactor.Allocate({Nphases,Ncomp,Ncomp});
    GrainBoundaryFactor.Allocate({Nphases,Nphases});

    locSettings.AddForRemeshing(*this);
    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void DiffusionProperties::ReadInput(string InputFileName)
{
    ConsoleOutput::WriteLineInsert("DiffusionProperties input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    ReadInput(inp);
    ConsoleOutput::WriteLine();
}

void DiffusionProperties::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);
    
    string extp = FileInterface::ReadParameterK(inp, moduleLocation, "ExtrapolationMode");

    if(extp == "NONE")
    {
        ExtrapolationMode = ExtrapolationModes::None;
    }
    else if (extp == "LOCAL")
    {
        ExtrapolationMode = ExtrapolationModes::Local;
    }
    else if (extp == "GLOBAL")
    {
        ExtrapolationMode = ExtrapolationModes::Global;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong extrapolation mode specified!", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }

    if(ExtrapolationMode != ExtrapolationModes::None)
    {
        MaxMoleFractionDeviation = FileInterface::ReadParameterD(inp, moduleLocation, "MDev", false, MaxMoleFractionDeviation);
        MaxTemperatureDeviation  = FileInterface::ReadParameterD(inp, moduleLocation, "TDev", false, MaxTemperatureDeviation);
    }

    for(size_t phase = 0; phase < Nphases; phase++)
    for(size_t  comp = 0;  comp <   Ncomp;  comp++)
    for(size_t  grad = 0;  grad <   Ncomp;  grad++)
    {
        string key = "ScalingFactor_" + to_string(phase) + "_" + ElementNames[comp] + "_" + ElementNames[grad];

        ScalingFactor({phase,comp,grad}) = FileInterface::ReadParameterD(inp, moduleLocation, key, false, 1.0);
    }

    for(size_t alpha =     0; alpha < Nphases; alpha++)
    for(size_t beta  = alpha;  beta < Nphases;  beta++)
    {
        string key = "GrainBoundaryFactor_" + to_string(alpha) + "_" + to_string(beta);

        GrainBoundaryFactor({alpha,beta}) = FileInterface::ReadParameterD(inp, moduleLocation, key, false, 1.0);
        GrainBoundaryFactor({beta,alpha}) = GrainBoundaryFactor({alpha,beta});
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void DiffusionProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(DiffusionCoefficients);
    BC.SetY(DiffusionCoefficients);
    BC.SetZ(DiffusionCoefficients);
}

void DiffusionProperties::Remesh(const int newNx, const int newNy, const int newNz, const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);
    DiffusionCoefficients.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    SetBoundaryConditions(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void DiffusionProperties::CheckRange(PhaseField& Phi,
                                     Composition& Cx,
                                     Temperature& Tx)
{
    switch(ExtrapolationMode)
    {
        case ExtrapolationModes::Global:
        {
            for(size_t m = 0; m < Nphases; m++)
            if(fabs(Tx.Tavg - GlobalExtrapolationData[m].eq_temperature) > MaxTemperatureDeviation)
            {
                GlobalExtrapolationData[m].out_of_T_range = true;
            }
            break;
        }
        case ExtrapolationModes::Local:
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
            {
                for(size_t m = 0; m < Nphases; m++)
                if(fabs(Tx.Tx(i,j,k) - LocalExtrapolationData(i,j,k,{m}).eq_temperature) > MaxTemperatureDeviation)
                {
                    LocalExtrapolationData(i,j,k,{m}).out_of_T_range = true;
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        case ExtrapolationModes::None:
        default:
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
            {
                for(size_t m = 0; m < Nphases; m++)
                {
                    LocalExtrapolationData(i,j,k,{m}).out_of_T_range = true;
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
    }
}

void DiffusionProperties::WriteVTK(const Settings& locSettings, const Composition& Cx, const int tStep)
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t phase = 0; phase < Nphases; ++phase)
    for(size_t comp  = 0; comp  < Ncomp;   ++comp)
    for(size_t grad  = 0; grad  < Ncomp;   ++grad)
    {
        const std::string name = "DiffusionCoefficient_" + to_string(phase) + "(" + ElementNames[comp] + "," + ElementNames[grad]+")";
        ListOfFields.push_back((VTK::Field_t) {name, [this, phase, comp, grad](int i,int j,int k){return DiffusionCoefficients(i,j,k,{phase,comp,grad});}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + '_', tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void DiffusionProperties::PrintPointStatistics(int i, int j, int k)
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t comp = 0; comp < Ncomp; comp++)
    for(size_t grad = 0; grad < Ncomp; grad++)
    {
        std::cout << std::fixed
                  << "Diffusion coefficient for phase " << n << " and comp "
                  << ElementNames[comp] << ",  grad "
                  << ElementNames[grad] << "): "
                  << std::scientific
                  << DiffusionCoefficients(i,j,k,{n,comp,grad}) << std::endl;
    }
}

}//namespace openphase
