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

 *   File created :   2016
 *   Main contributors :   Raphael Schiedung
 *
 */

#include "BoundaryConditions.h"
#include "Electrics/ElectricProperties.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Composition.h"
#include "VTK.h"

namespace openphase
{

ElectricProperties::ElectricProperties(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void ElectricProperties::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ElectricProperties";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;
    Nphases = locSettings.Nphases;
    Ncomp   = locSettings.Ncomp;
    size_t Bcells = Grid.Bcells;

    PolarizationDensity.Allocate(Grid, Bcells);
    ElectricField      .Allocate(Grid, Bcells);
    ElectricPotential  .Allocate(Grid, Bcells);
    ChargeDensity      .Allocate(Grid, Bcells);
    WaveVector         .Allocate(Grid, Bcells);

    MolarCharge.Allocate(Ncomp);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ElectricProperties::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Electric properties input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::fstream inpF(InputFileName.c_str(), std::ios::in | std::ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream inp;
    inp << inpF.rdbuf();

    ReadInput(inp);

    inpF.close();
}

void ElectricProperties::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        const std::string name = "MolarCharge_" + comp;
        MolarCharge[comp] = FileInterface::ReadParameterD(inp, moduleLocation, name, false, 0.0);
    }
}

void ElectricProperties::CalculateChargeDensity(const Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChargeDensity,ChargeDensity.Bcells(),)
    {
        ChargeDensity(i,j,k) = 0.0;
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            const double MolarComponentDensity = Cx.MoleFractionsTotal(i,j,k,{comp})*Cx.PT.GetData(Cx.Component[comp].Name).MolarVolume;
            ChargeDensity(i,j,k) += MolarCharge[comp] * MolarComponentDensity;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElectricProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(PolarizationDensity);
    BC.SetY(PolarizationDensity);
    BC.SetZ(PolarizationDensity);

    BC.SetX(ElectricField);
    BC.SetY(ElectricField);
    BC.SetZ(ElectricField);

    BC.SetX(ElectricPotential);
    BC.SetY(ElectricPotential);
    BC.SetZ(ElectricPotential);

    BC.SetX(ChargeDensity);
    BC.SetY(ChargeDensity);
    BC.SetZ(ChargeDensity);
}

void ElectricProperties::WriteVTK(const Settings& locSettings, const int tStep, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Electric Charge Density [A*s/m^3]", [this](int i,int j,int k){return ChargeDensity(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t) {"Electric Potential [V]", [this](int i,int j,int k){return ElectricPotential(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t) {"Wave Vector [1/m]", [this](int i,int j,int k){return WaveVector(i,j,k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ElectricProperties", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

}// end of namespace openphase
