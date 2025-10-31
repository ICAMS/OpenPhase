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

 *   File created :   2017
 *   Main contributors :   Raphael Schiedung
 *
 */

#include "BoundaryConditions.h"
#include "Magnetism/MagneticProperties.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Settings.h"
#include "VTK.h"

using namespace std;

namespace openphase
{

MagneticProperties::MagneticProperties(Settings& locSettings)
{
    Initialize(locSettings);
    ReadInput(DefaultInputFileName);
}

void MagneticProperties::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "MagneticProperties";
    thisobjectname = thisclassname + ObjectNameSuffix;
    //DefaultInputFileName = ProjectInputDir + "MagneticInput.opi";

    const unsigned int BoundaryCells = 1;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;

    Permeability.Allocate(Nphases);
    PermeabilityPhase.Allocate(Nphases);
    RelativePermeability.Allocate(Nphases);
    RelativePermeabilityPhase.Allocate(Nphases);

    MagneticField.Allocate(Grid, BoundaryCells);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void MagneticProperties::ReadInput(std::string InputFileName)
{
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert("Magnetic properties");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened",
                "ReadInput");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream data;
    data << inp.rdbuf();

    ReadInput(data);

    inp.close();

    ConsoleOutput::WriteLine();
}

void MagneticProperties::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    const unsigned int Sublattices =
        FileInterface::ReadParameterI(inp, moduleLocation, std::string("SUB"));

    Magnetisation.Allocate(Grid, {Sublattices}, 1);

    for(unsigned long int alpha = 0; alpha < Nphases; alpha++)
    for(unsigned long int i = 0; i < 3; i++)
    for(unsigned long int j = 0; j < 3; j++)
    {
        stringstream converter;
        converter << alpha << "_" << i+1 << j+1;
        const std::string name = "PMEA" + converter.str();

        const double value = FileInterface::ReadParameterD(inp, moduleLocation, name, false, false);

        RelativePermeability     [alpha](i,j) = value;
        RelativePermeabilityPhase[alpha](i,j) = value;
        Permeability             [alpha](i,j) = value * PhysicalConstants::mu_0;
        PermeabilityPhase        [alpha](i,j) = value * PhysicalConstants::mu_0;
    }

    ConsoleOutput::WriteLine("-");
}

void MagneticProperties::SetGrainsProperties(const PhaseField& Phase)
{
    const unsigned int NPhaseFields = Phase.FieldsProperties.size();
    if(Permeability.size() != NPhaseFields)
    {
        Permeability.Reallocate(NPhaseFields);
        RelativePermeability.Reallocate(NPhaseFields);
    }

    for(unsigned int alpha = 0; alpha < NPhaseFields; alpha++)
    if(Phase.FieldsProperties[alpha].Exist)
    {
        const unsigned int pIndex = Phase.FieldsProperties[alpha].Phase;
        Permeability        [alpha] = PermeabilityPhase        [pIndex];
        RelativePermeability[alpha] = RelativePermeabilityPhase[pIndex];
    }
}

void MagneticProperties::SetGrainsProperties(const PhaseField& Phase,
        const Orientations& OR)
{
    const unsigned int NPhaseFields = Phase.FieldsProperties.size();
    if(Permeability.size() != NPhaseFields)
    {
        Permeability.Reallocate(NPhaseFields);
        RelativePermeability.Reallocate(NPhaseFields);
    }

    /*for(unsigned int alpha = 0; alpha < NPhaseFields; alpha++)
    if(Phase.FieldsProperties[alpha].Exist)
    {
        const unsigned int pIndex = Phase.FieldsProperties[alpha].Phase;
        //Permeability        [alpha] = PermeabilityPhase        [pIndex].rotated(OR.GrainRotations[alpha]);
        //RelativePermeability[alpha] = RelativePermeabilityPhase[pIndex].rotated(OR.GrainRotations[alpha]);
    }*/
}

void MagneticProperties::WriteMagneticFieldVTK(const int tStep, const Settings& locSettings) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "MagneticField_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"MagneticField", [this](int i,int j,int k){return MagneticField(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields);
}

void MagneticProperties::WriteMagnetisationVTK(const int tStep, const Settings& locSettings) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Magnetisation_", tStep, ".vts");
    for (unsigned int n = 0; n < Sublattices; n++)
    {
        ListOfFields.push_back((VTK::Field_t) {"Magnetisation_" + to_string(n), [this,n](int i,int j,int k){return Magnetisation(i,j,k,{n});}});
    }
    VTK::Write(Filename, locSettings, ListOfFields);
}

void MagneticProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetXVector(MagneticField);
    BC.SetYVector(MagneticField);
    BC.SetZVector(MagneticField);

    BC.SetXVector(Magnetisation);
    BC.SetYVector(Magnetisation);
    BC.SetZVector(Magnetisation);
}
}// end of namespace openphase
