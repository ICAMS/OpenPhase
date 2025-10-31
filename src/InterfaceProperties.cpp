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

 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Raphael Schiedung
 *
 */

#include "InterfaceProperties.h"
#include "BoundaryConditions.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "ElasticProperties.h"
#include "VTK.h"

namespace openphase
{
using namespace std;
//Auxiliary function to rotate grain interface normal back to the reference configuration
dVector3 InterfaceOrientation(const PhaseField& Phase, NodePF::citerator it)
{
    dVector3 normal = it->gradient*(-1.0);
    normal.normalize();
    normal = Phase.FieldsProperties[it->index].Orientation.RotationMatrix.transposed()*normal;
    return normal;
}

void InterfaceProperties::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "InterfaceProperties";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    ConsoleOutput::WriteBlankLine();
    Nphases = locSettings.Nphases;

    TripleJunctionEnergy = 0.0;
    TripleJunctionFactor = 0.0;
    RegularizationFactor = 1.0;
    CurvatureFactor      = 1.0;

    maxEnergies.Allocate(Nphases,Nphases);
    maxMobilities.Allocate(Nphases,Nphases);

    NucleiMobilityFactor.Allocate(Nphases, Nphases);
    for (size_t i = 0; i < Nphases; ++i)
    for (size_t j = 0; j < Nphases; ++j)
    {
        NucleiMobilityFactor(i,j) = 1.0;
    }

    FullAnisotropy = false;

    size_t Bcells = Grid.Bcells;
    Properties.Allocate(Grid, Bcells);

    if(Grid.Resolution == Resolutions::Dual)
    {
        GridParameters DoubleDimensions = Grid.DoubleResolution();
        PropertiesDR.Allocate(DoubleDimensions, Bcells*2);
    }
    InterfaceEnergy.Allocate(Nphases, Nphases);
    InterfaceMobility.Allocate(Nphases, Nphases);
    RespectParentBoundaries.Allocate(Nphases, Nphases);

    //PhaseActivationModes = locSettings.PhaseActivationModes;
    //PhaseInteractions = locSettings.PhaseInteractions;

    locSettings.AddForRemeshing(*this);
    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void InterfaceProperties::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName);
    std::string filetype = FileInterface::getFileExtension(InputFileName);
    if (filetype == "opi")
    {
        fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

        if (!inp)
        {
            ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        };

        std::stringstream data;
        data << inp.rdbuf();
        ReadInput(data);

        inp.close();
    }
    else
    if (filetype == "json")
    {
        ReadJSON(InputFileName);
    }
    else
    {
        std::cerr << "Filetype " << filetype << " not recognized. Filetype must be opi or json." << std::endl;
        OP_Exit(1);
    }
}

void InterfaceProperties::ReadInput(stringstream& inp)
{
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert("Interface Properties");

    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

//    // Reading interaction modes for phase pairs.
//    PhaseInteractions.Allocate(Nphases, Nphases);
//    for(size_t alpha =     0; alpha < Nphases; alpha++)
//    for(size_t beta  = alpha; beta  < Nphases; beta++)
//    if(PhaseActivationModes[alpha] == ActivationModes::Enabled and
//       PhaseActivationModes[beta]  == ActivationModes::Enabled)
//    {
//        stringstream converter;
//        converter << string("PhaseInteractions_") << alpha << "_" << beta;
//        PhaseInteractions(alpha,beta) = FileInterface::ReadParameterB(inp, moduleLocation, converter.str(), false, true);
//        PhaseInteractions(beta,alpha) = PhaseInteractions(alpha,beta);
//    }
//    else
//    {
//        PhaseInteractions(alpha,beta) = false;
//        PhaseInteractions(beta,alpha) = false;
//    }

    // Reading interface energies and mobilities
    for (size_t alpha =     0; alpha < Nphases; alpha++)
    for (size_t beta  = alpha;  beta < Nphases; beta++)
    {
        stringstream converter;
        converter << "_" << alpha << "_" << beta;
        string counter = converter.str();

        string locEnergyAnisotropyModel = FileInterface::ReadParameterK(inp, moduleLocation, "EnergyModel" + counter, false, "ISO");

        if (locEnergyAnisotropyModel == "EXT")
        {
            InterfaceEnergy(alpha,beta).Model = InterfaceEnergyModels::Ext;
        }
        if (locEnergyAnisotropyModel == "ISO")
        {
            InterfaceEnergy(alpha,beta).ReadInputIso(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "CUBIC")
        {
            InterfaceEnergy(alpha,beta).ReadInputCubic(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "CUBICFULL")
        {
            InterfaceEnergy(alpha,beta).ReadInputCubicFull(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "HEXBOETTGER")
        {
            InterfaceEnergy(alpha,beta).ReadInputHexBoettger(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "HEXYANG")
        {
            InterfaceEnergy(alpha,beta).ReadInputHexYang(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "HEXSUN")
        {
            InterfaceEnergy(alpha,beta).ReadInputHexSun(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "FACETED")
        {
            InterfaceEnergy(alpha,beta).ReadInputFaceted(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "FACETEDFULL")
        {
            InterfaceEnergy(alpha,beta).ReadInputFacetedFull(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "FACETEDEXP")
        {
            InterfaceEnergy(alpha,beta).ReadInputFacetedEXP(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "DISORIENTATION")
        {
            InterfaceEnergy(alpha,beta).ReadInputDisorientation(inp, moduleLocation, counter);
        }
        string locMobilityAnisotropyModel = FileInterface::ReadParameterK(inp, moduleLocation, "MobilityModel" + counter, false, "ISO");

        if (locMobilityAnisotropyModel == "EXT")
        {
            InterfaceMobility(alpha,beta).Model = InterfaceMobilityModels::Ext;
        }
        if (locMobilityAnisotropyModel == "ISO")
        {
            InterfaceMobility(alpha,beta).ReadInputIso(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "CUBIC")
        {
            InterfaceMobility(alpha,beta).ReadInputCubic(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "HEXBOETTGER")
        {
            InterfaceMobility(alpha,beta).ReadInputHexBoettger(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "HEXYANG")
        {
            InterfaceMobility(alpha,beta).ReadInputHexYang(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "HEXSUN")
        {
            InterfaceMobility(alpha,beta).ReadInputHexSun(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "FACETED")
        {
            InterfaceMobility(alpha,beta).ReadInputFaceted(inp, moduleLocation, counter);
        }

        RespectParentBoundaries(alpha, beta) = FileInterface::ReadParameterB(inp, moduleLocation, string("RPB") + counter, false, false);

        if(alpha != beta)
        {
            InterfaceEnergy(beta,alpha) = InterfaceEnergy(alpha,beta);
            InterfaceMobility(beta,alpha) = InterfaceMobility(alpha,beta);
            RespectParentBoundaries(beta, alpha) = RespectParentBoundaries(alpha, beta);
        }
    }

    // Reading mobility correction factor for nuclei
    for (size_t alpha = 0; alpha < Nphases; alpha++)
    for (size_t beta  = 0;  beta < Nphases; beta++)
    {
        stringstream converter;
        converter << "_" << alpha << "_" << beta;
        string counter = converter.str();

        NucleiMobilityFactor(alpha, beta) = FileInterface::ReadParameterD(inp, moduleLocation, string("NucleiMobilityFactor") + counter, false, 1.0);
    }

    // Setting interface energy model type
    for (size_t i = 0; i < Nphases; i++)
    for (size_t j = 0; j < Nphases; j++)
    {
        if (InterfaceEnergy(i,j).Model != InterfaceEnergyModels::Iso and
            InterfaceEnergy(i,j).Type  == InterfaceEnergyModelType::Energy)
        {
            FullAnisotropy = true;
        }
    }

    if(FullAnisotropy)
    {
        InterfaceStiffnessTMP.Allocate(Grid, Properties.Bcells());
    }

    // Reading triple junction models
    string locTripleJunctionModel = FileInterface::ReadParameterK(inp, moduleLocation, "TripleJunctionModel", false, "MODEL");

    if (locTripleJunctionModel == "MODEL")
    {
        TripleJunctionModel = TripleJunctionModels::Model;
        TripleJunctionFactor = FileInterface::ReadParameterD(inp, moduleLocation, string("TripleJunctionFactor"), false, 0.0);
    }
    if (locTripleJunctionModel == "VALUE")
    {
        TripleJunctionModel = TripleJunctionModels::Value;
        TripleJunctionEnergy = FileInterface::ReadParameterD(inp, moduleLocation, string("TripleJunctionEnergy"), true, 0.0);
    }
    if (locTripleJunctionModel == "NONE")
    {
        TripleJunctionModel = TripleJunctionModels::None;
    }

    RegularizationFactor = FileInterface::ReadParameterD(inp, moduleLocation, string("RegularizationFactor"), false, 1.0);
    CurvatureFactor      = FileInterface::ReadParameterD(inp, moduleLocation, string("CurvatureFactor"),      false, 1.0);

    std::string locExtrapolationMode = FileInterface::ReadParameterK(inp, moduleLocation, std::string("ExtrapolationMode"), false, std::string("NONE"));
    if (locExtrapolationMode == "NONE")
    {
        ExtrapolationMode = ExtrapolationModes::None;
        ConsoleOutput::WriteWarning("Extrapolation mode set to NONE. Curvature will be assumed zero at some points when it is actually not!!! ", thisclassname, "ReadInput()");
    }
    else if (locExtrapolationMode == "ZEROTHORDER")
    {
        ExtrapolationMode = ExtrapolationModes::ZerothOrder;
    }
    else if (locExtrapolationMode == "FIRSTORDER")
    {
        ExtrapolationMode = ExtrapolationModes::FirstOrder;
    }
    else 
    {
        ConsoleOutput::WriteExit("Wrong extrapolation mode selected -> " + locExtrapolationMode, thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }
    if (ExtrapolationMode != ExtrapolationModes::None) 
    {
        if(Grid.Resolution == Resolutions::Dual)
        {
            PropertiesExtrapolations.Allocate(Grid.DoubleResolution(), Grid.Bcells*2);
        }
        else
        {
            PropertiesExtrapolations.Allocate(Grid, Grid.Bcells);
        }
    }
    MaxExtrapolationInterations = FileInterface::ReadParameterI(inp, moduleLocation, std::string("MaxExtrapolationInterations"), false, 6);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void InterfaceProperties::ReadJSON(const string InputFileName)
{
    std::ifstream f(InputFileName);
    json data = json::parse(f);
    if (data.contains(thisclassname))
    {
        json interfaceproperties = data[thisclassname];

        for (size_t alpha =     0; alpha < Nphases; alpha++)
        for (size_t beta  = alpha;  beta < Nphases; beta++)
        {
            std::string locEnergyAnisotropyModel = FileInterface::ReadParameter<std::string>(interfaceproperties, {"EnergyModel", alpha, beta}, "ISO");

            if (locEnergyAnisotropyModel == "EXT")
            {
                InterfaceEnergy(alpha,beta).Model = InterfaceEnergyModels::Ext;
            }
            if (locEnergyAnisotropyModel == "ISO")
            {
                InterfaceEnergy(alpha,beta).ReadJSONIso(interfaceproperties, alpha, beta);
            }
            if (locEnergyAnisotropyModel == "CUBIC")
            {
                InterfaceEnergy(alpha,beta).ReadJSONCubic(interfaceproperties, alpha, beta);
            }
            if (locEnergyAnisotropyModel == "CUBICFULL")
            {
                InterfaceEnergy(alpha,beta).ReadJSONCubicFull(interfaceproperties, alpha, beta);
            }
            if (locEnergyAnisotropyModel == "HEXBOETTGER")
            {
                InterfaceEnergy(alpha,beta).ReadJSONHexBoettger(interfaceproperties, alpha, beta);
            }
            if (locEnergyAnisotropyModel == "HEXYANG")
            {
                InterfaceEnergy(alpha,beta).ReadJSONHexYang(interfaceproperties, alpha, beta);
            }
            if (locEnergyAnisotropyModel == "HEXSUN")
            {
                InterfaceEnergy(alpha,beta).ReadJSONHexSun(interfaceproperties, alpha, beta);
            }

            std::string locMobilityAnisotropyModel = FileInterface::ReadParameter<std::string>(interfaceproperties, {"MobilityModel", alpha, beta}, "ISO");

            if (locMobilityAnisotropyModel == "EXT")
            {
                InterfaceMobility(alpha,beta).Model = InterfaceMobilityModels::Ext;
            }
            if (locMobilityAnisotropyModel == "ISO")
            {
                InterfaceMobility(alpha,beta).ReadJSONIso(interfaceproperties, alpha, beta);
            }
            if (locMobilityAnisotropyModel == "CUBIC")
            {
                InterfaceMobility(alpha,beta).ReadJSONCubic(interfaceproperties, alpha, beta);
            }
            if (locMobilityAnisotropyModel == "HEXBOETTGER")
            {
                InterfaceMobility(alpha,beta).ReadJSONHexBoettger(interfaceproperties, alpha, beta);
            }
            if (locMobilityAnisotropyModel == "HEXYANG")
            {
                InterfaceMobility(alpha,beta).ReadJSONHexYang(interfaceproperties, alpha, beta);
            }
            if (locMobilityAnisotropyModel == "HEXSUN")
            {
                InterfaceMobility(alpha,beta).ReadJSONHexSun(interfaceproperties, alpha, beta);
            }

            RespectParentBoundaries(alpha, beta) = FileInterface::ReadParameter<bool>(interfaceproperties, {"RPB", alpha, beta}, false);

            if(alpha != beta)
            {
                InterfaceEnergy(beta,alpha) = InterfaceEnergy(alpha,beta);
                InterfaceMobility(beta,alpha) = InterfaceMobility(alpha,beta);
                RespectParentBoundaries(beta, alpha) = RespectParentBoundaries(alpha, beta);
            }
        }

        // Reading mobility correction factor for nuclei
        for (size_t alpha = 0; alpha < Nphases; alpha++)
        for (size_t beta  = 0;  beta < Nphases; beta++)
        {
            stringstream converter;
            converter << "_" << alpha << "_" << beta;
            string counter = converter.str();

            NucleiMobilityFactor(alpha, beta) = FileInterface::ReadParameter<double>(interfaceproperties, {"NucleiMobilityFactor", alpha, beta}, 1.0);
        }

        // Setting interface energy model type
        for (size_t i = 0; i < Nphases; i++)
        for (size_t j = 0; j < Nphases; j++)
        {
            if (InterfaceEnergy(i,j).Model != InterfaceEnergyModels::Iso and
                InterfaceEnergy(i,j).Type  == InterfaceEnergyModelType::Energy)
            {
                FullAnisotropy = true;
            }
        }

        if(FullAnisotropy)
        {
            InterfaceStiffnessTMP.Allocate(Grid, Properties.Bcells());
        }

        // Reading triple junction models
        string locTripleJunctionModel = FileInterface::ReadParameter<std::string>(interfaceproperties, {"TripleJunctionModel"}, "MODEL");

        if (locTripleJunctionModel == "MODEL")
        {
            TripleJunctionModel = TripleJunctionModels::Model;
            TripleJunctionFactor = FileInterface::ReadParameter<double>(interfaceproperties, {"TripleJunctionFactor"}, 0.0);
        }
        if (locTripleJunctionModel == "VALUE")
        {
            TripleJunctionModel = TripleJunctionModels::Value;
            TripleJunctionEnergy = FileInterface::ReadParameter<double>(interfaceproperties, {"TripleJunctionEnergy"}, 0.0);
        }
        if (locTripleJunctionModel == "NONE")
        {
            TripleJunctionModel = TripleJunctionModels::None;
        }

        RegularizationFactor = FileInterface::ReadParameter<double>(interfaceproperties, {"RegularizationFactor"}, 1.0);
        CurvatureFactor      = FileInterface::ReadParameter<double>(interfaceproperties, {"CurvatureFactor"}, 1.0);
    }
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void InterfaceProperties::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);

    Properties.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    if(InterfaceStiffnessTMP.IsAllocated())
    {
        InterfaceStiffnessTMP.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    }

    if(Grid.Resolution == Resolutions::Dual)
    {
        PropertiesDR.Reallocate((1+Grid.dNx)*Grid.Nx, (1+Grid.dNy)*Grid.Ny, (1+Grid.dNz)*Grid.Nz);
    }

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void InterfaceProperties::Coarsen(const PhaseField& Phase)
{
    double norm = 1.0/pow(2.0,Grid.Active());
    long int fx = 1 + Grid.dNx;
    long int fy = 1 + Grid.dNy;
    long int fz = 1 + Grid.dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).wide_interface())
        {
            Properties(i,j,k).clear();

            for(int di = -Grid.dNx; di <= Grid.dNx; di+=2)
            for(int dj = -Grid.dNy; dj <= Grid.dNy; dj+=2)
            for(int dk = -Grid.dNz; dk <= Grid.dNz; dk+=2)
            {
                Properties(i,j,k).add_energies_and_mobilities(PropertiesDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2));
            }
            Properties(i,j,k) *= norm;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceProperties::Refine(const PhaseField& Phase)
{
    long int fx = 1 + Grid.dNx;
    long int fy = 1 + Grid.dNy;
    long int fz = 1 + Grid.dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Properties,0,)
    {
        if(Phase.Fields(i,j,k).wide_interface())
        for(int di = -Grid.dNx; di <= Grid.dNx; di+=2)
        for(int dj = -Grid.dNy; dj <= Grid.dNy; dj+=2)
        for(int dk = -Grid.dNz; dk <= Grid.dNz; dk+=2)
        {
            PropertiesDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2) = Properties.at(i+di*0.25,j+dj*0.25,k+dk*0.25);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double InterfaceProperties::ReportMaximumTimeStep(void)
{
    double maxEnergy = 0.0;
    double maxMobility = 0.0;
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    {
        maxEnergy   = std::max(maxEnergy, maxEnergies(n,m));
        maxMobility = std::max(maxMobility, maxMobilities(n,m));
    }
    double max_dt = DBL_MAX;
    if(maxEnergy > DBL_EPSILON and maxMobility > 0.0)
    {
        switch(Grid.Active())
        {
            case 1:
            {
                max_dt = 0.50*Grid.dx*Grid.dx/(maxMobility*maxEnergy*RegularizationFactor);
                break;
            }
            case 2:
            {
                max_dt = 0.25*Grid.dx*Grid.dx/(maxMobility*maxEnergy*RegularizationFactor);
                break;
            }
            case 3:
            {
                max_dt = 0.16*Grid.dx*Grid.dx/(maxMobility*maxEnergy*RegularizationFactor);
                break;
            }
        }
    }
    return max_dt;
}
void InterfaceProperties::Set(const PhaseField& Phase, const BoundaryConditions& BC)
{
    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            SetSR(Phase);
            break;
        }
        case Resolutions::Dual:
        {
            SetDR(Phase);
            SetBoundaryConditionsDR(BC);
            Coarsen(Phase);
            break;
        }
    }
    SetBoundaryConditions(BC);
    Extrapolate(Phase, BC);
}

void InterfaceProperties::Set(const PhaseField& Phase,
                              const Temperature& Tx,
                              const BoundaryConditions& BC)
{
    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            SetSR(Phase);
            SetMobilityThermalEffectSR(Phase, Tx);
            break;
        }
        case Resolutions::Dual:
        {
            SetDR(Phase);
            SetMobilityThermalEffectDR(Phase, Tx);
            SetBoundaryConditionsDR(BC);
            Coarsen(Phase);
            break;
        }
    }
    SetBoundaryConditions(BC);
    Extrapolate(Phase, BC);
}

void InterfaceProperties::SetSR(const PhaseField& Phase)
{
    Matrix<double> locMaxEnergies(Nphases,Nphases);
    Matrix<double> locMaxMobilities(Nphases,Nphases);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Properties, 0, reduction(MatrixDMAX:locMaxEnergies) reduction(MatrixDMAX:locMaxMobilities))
    {
        Properties(i,j,k).clear();

        if(Phase.Fields(i,j,k).wide_interface())
        {
            for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            if(alpha != beta)
            {
                int pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                int pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                // Setting interface energy/stiffness
                double locEnergy = 0.0;
                if(InterfaceEnergy(pIndexA, pIndexB).Model == InterfaceEnergyModels::Iso)
                {
                    locEnergy = InterfaceEnergy(pIndexA, pIndexB).Energy;
                }
                else if(InterfaceEnergy(pIndexA, pIndexB).Model == InterfaceEnergyModels::Disorientation)
                {
                    Quaternion OrientationA = Phase.FieldsProperties[alpha->index].Orientation;
                    Quaternion OrientationB = Phase.FieldsProperties[ beta->index].Orientation;
                    locEnergy = InterfaceEnergy(pIndexA, pIndexB).CalculateDisorientation(OrientationA, OrientationB);
                }
                else if(InterfaceEnergy(pIndexA, pIndexB).Model != InterfaceEnergyModels::Ext)
                {
                    double weightA = 0.0;
                    double weightB = 0.0;

                    if(Phase.FieldsProperties[alpha->index].State == AggregateStates::Solid)
                    {
                        weightA = fabs(alpha->value*(1.0 - alpha->value));
                        if(weightA > DBL_EPSILON)
                        {
                            dVector3 Normal_alpha = InterfaceOrientation(Phase,alpha);
                            locEnergy += InterfaceEnergy(pIndexA, pIndexB).Calculate(Normal_alpha)*weightA;
                        }
                        else
                        {
                            weightA = 0.0;
                        }
                    }

                    if(Phase.FieldsProperties[beta->index].State == AggregateStates::Solid)
                    {
                        weightB = fabs(beta->value*(1.0 - beta->value));

                        if(weightB > DBL_EPSILON)
                        {
                            dVector3 Normal_beta = InterfaceOrientation(Phase,beta);
                            locEnergy += InterfaceEnergy(pIndexA, pIndexB).Calculate(Normal_beta)*weightB;
                        }
                        else
                        {
                            weightB = 0.0;
                        }
                    }

                    if (weightA + weightB > DBL_EPSILON)
                    {
                        locEnergy /= weightA + weightB;
                    }
                    else
                    {
                        locEnergy = InterfaceEnergy(pIndexA, pIndexB).MaxEnergy;
                    }
                }

                // Setting interface mobility
                double locMobility = 0.0;
                if(InterfaceMobility(pIndexA, pIndexB).Model == InterfaceMobilityModels::Iso)
                {
                    locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                }
                else if (InterfaceMobility(pIndexA, pIndexB).Model != InterfaceMobilityModels::Ext)
                {
                    int number_of_solids = 0;
                    if(Phase.FieldsProperties[alpha->index].State == AggregateStates::Solid)
                    if(alpha->value != 0.0 and alpha->value != 1.0)
                    {
                        dVector3 Normal_alpha = InterfaceOrientation(Phase,alpha);
                        locMobility += InterfaceMobility(pIndexA, pIndexB).Calculate(Normal_alpha);
                        number_of_solids ++;
                    }

                    if(Phase.FieldsProperties[beta->index].State == AggregateStates::Solid)
                    if(beta->value != 0.0 and beta->value != 1.0)
                    {
                        dVector3 Normal_beta = InterfaceOrientation(Phase,beta);
                        locMobility += InterfaceMobility(pIndexA, pIndexB).Calculate(Normal_beta);
                        number_of_solids ++;
                    }

                    if(number_of_solids == 0)
                    {
                        locMobility = InterfaceMobility(pIndexA, pIndexB).Mobility;
                    }
                }

                if(RespectParentBoundaries(pIndexA, pIndexB))
                {
                    if(pIndexA == pIndexB and
                       Phase.FieldsProperties[alpha->index].Parent != Phase.FieldsProperties[ beta->index].Parent)
                    {
                        locMobility = 0.0;
                    }

                    if(pIndexA != pIndexB and
                       Phase.FieldsProperties[alpha->index].Parent !=  beta->index and
                       Phase.FieldsProperties[ beta->index].Parent != alpha->index)
                    {
                        locMobility = 0.0;
                    }
                }

                // Modifying mobility for nuclei (if enabled)
                if(NucleiMobilityFactor(pIndexA, pIndexB) != 1.0 and
                   (Phase.FieldsProperties[alpha->index].IsNucleus() or
                    Phase.FieldsProperties[ beta->index].IsNucleus()))
                {
                    double volume_ratio = min(Phase.FieldsProperties[alpha->index].VolumeRatio,
                                              Phase.FieldsProperties[ beta->index].VolumeRatio);

                    locMobility *= NucleiMobilityFactor(pIndexA, pIndexB)*(1.0 - volume_ratio) + 1.0*volume_ratio;
                }

                // Modifying mobility if growth constraints are violated
                if(Phase.FieldsProperties[alpha->index].GrowthConstraintsViolation != GrowthConstraintsViolations::None and
                   Phase.FieldsProperties[ beta->index].GrowthConstraintsViolation != GrowthConstraintsViolations::None)
                {
                    locMobility *= Phase.PairwiseGrowthFactors.get_sym1(alpha->index,beta->index);
                }

                // Putting interface energy and mobility values into the properties storage
                if(locEnergy != 0.0 or locMobility != 0.0)
                {
                    Properties(i,j,k).set_energy_and_mobility(alpha->index, beta->index, locEnergy, locMobility);
                }
                locMaxEnergies(pIndexA,pIndexB) = max(locMaxEnergies(pIndexA,pIndexB),locEnergy);
                locMaxMobilities(pIndexA,pIndexB) = max(locMaxMobilities(pIndexA,pIndexB),locMobility);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxEnergies = locMaxEnergies;
    maxMobilities = locMaxMobilities;

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, maxEnergies.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, maxMobilities.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
}

dVector3 InterfaceProperties::dEnergy_dGradientAlpha(const PhaseField& Phase, NodePF::citerator alpha, NodePF::citerator beta) const
{
    size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
    size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

    dVector3 value = {0,0,0};

    if(InterfaceEnergy(pIndexA, pIndexB).Model != InterfaceEnergyModels::Iso and
       InterfaceEnergy(pIndexA, pIndexB).Model != InterfaceEnergyModels::Ext)
    {
        double number_of_solids = 0.0;
        if(Phase.FieldsProperties[alpha->index].State == AggregateStates::Solid)
        if(alpha->value != 0.0 and alpha->value != 1.0)
        {
            const dMatrix3x3 R = Phase.FieldsProperties[alpha->index].Orientation.RotationMatrix.transposed();
            const dVector3 normal = alpha->gradient.normalized()*(-1);
            const dVector3 locInterfaceOrientationAlpha = R*normal;
            const dMatrix3x3 dlocInterfaceOrientationAlpha_dGradientAlpha = (R*(normal.dyadic(normal) - dMatrix3x3::UnitTensor()))/alpha->gradient.length();
            const dVector3 dsigma_dInterfaceOrientationAlpha = InterfaceEnergy(pIndexA, pIndexB).Derivative(locInterfaceOrientationAlpha);
            value += dlocInterfaceOrientationAlpha_dGradientAlpha*dsigma_dInterfaceOrientationAlpha;
            number_of_solids ++;
        }

        if(Phase.FieldsProperties[beta->index].State == AggregateStates::Solid)
        if(beta->value != 0.0 and beta->value != 1.0)
        {
            number_of_solids ++;
        }

        if(number_of_solids != 0.0)
        {
            value /= number_of_solids;
        }
    }
    return value;
}

void InterfaceProperties::SetDR(const PhaseField& Phase)
{
    Matrix<double> locMaxEnergies(Nphases,Nphases);
    Matrix<double> locMaxMobilities(Nphases,Nphases);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PropertiesDR, 0, reduction(MatrixDMAX:locMaxEnergies) reduction(MatrixDMAX:locMaxMobilities))
    {
        PropertiesDR(i,j,k).clear();

        if(Phase.FieldsDR(i,j,k).wide_interface())
        {
            for(auto alpha  = Phase.FieldsDR(i,j,k).cbegin();
                     alpha != Phase.FieldsDR(i,j,k).cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.FieldsDR(i,j,k).cend(); ++beta)
            if(alpha != beta)
            {
                int pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                int pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                // Setting interface energy/stiffness
                double locEnergy = 0.0;
                if(InterfaceEnergy(pIndexA, pIndexB).Model == InterfaceEnergyModels::Iso)
                {
                    locEnergy = InterfaceEnergy(pIndexA, pIndexB).MaxEnergy;
                }
                else if(InterfaceEnergy(pIndexA, pIndexB).Model == InterfaceEnergyModels::Disorientation)
                {
                    Quaternion OrientationA = Phase.FieldsProperties[alpha->index].Orientation;
                    Quaternion OrientationB = Phase.FieldsProperties[ beta->index].Orientation;
                    locEnergy = InterfaceEnergy(pIndexA, pIndexB).CalculateDisorientation(OrientationA, OrientationB);
                }
                else if(InterfaceEnergy(pIndexA, pIndexB).Model != InterfaceEnergyModels::Ext)
                {
                    double weightA = 0.0;
                    double weightB = 0.0;

                    if(Phase.FieldsProperties[alpha->index].State == AggregateStates::Solid)
                    {
                        weightA = fabs(alpha->value*(1.0 - alpha->value));
                        if(weightA > DBL_EPSILON)
                        {
                            dVector3 Normal_alpha = InterfaceOrientation(Phase,alpha);
                            locEnergy += InterfaceEnergy(pIndexA, pIndexB).Calculate(Normal_alpha)*weightA;
                        }
                        else
                        {
                            weightA = 0.0;
                        }
                    }

                    if(Phase.FieldsProperties[beta->index].State == AggregateStates::Solid)
                    {
                        weightB = fabs(beta->value*(1.0 - beta->value));

                        if(weightB > DBL_EPSILON)
                        {
                            dVector3 Normal_beta = InterfaceOrientation(Phase,beta);
                            locEnergy += InterfaceEnergy(pIndexA, pIndexB).Calculate(Normal_beta)*weightB;
                        }
                        else
                        {
                            weightB = 0.0;
                        }
                    }

                    if (weightA + weightB > DBL_EPSILON)
                    {
                        locEnergy /= weightA + weightB;
                    }
                    else
                    {
                        locEnergy = InterfaceEnergy(pIndexA, pIndexB).MaxEnergy;
                    }
                }

                // Setting interface mobility
                double locMobility = 0.0;
                if(InterfaceMobility(pIndexA, pIndexB).Model == InterfaceMobilityModels::Iso)
                {
                    locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                }
                else if (InterfaceMobility(pIndexA, pIndexB).Model != InterfaceMobilityModels::Ext)
                {
                    int number_of_solids = 0;
                    if(Phase.FieldsProperties[alpha->index].State == AggregateStates::Solid)
                    if(alpha->value != 0.0 and alpha->value != 1.0)
                    {
                        dVector3 Normal_alpha = InterfaceOrientation(Phase,alpha);
                        locMobility += InterfaceMobility(pIndexA, pIndexB).Calculate(Normal_alpha);
                        number_of_solids ++;
                    }

                    if(Phase.FieldsProperties[beta->index].State == AggregateStates::Solid)
                    if(beta->value != 0.0 and beta->value != 1.0)
                    {
                        dVector3 Normal_beta = InterfaceOrientation(Phase,beta);
                        locMobility += InterfaceMobility(pIndexA, pIndexB).Calculate(Normal_beta);
                        number_of_solids ++;
                    }

                    if(number_of_solids == 0)
                    {
                        locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                    }
                }

                if(RespectParentBoundaries(pIndexA, pIndexB))
                {
                    if(pIndexA == pIndexB and
                       Phase.FieldsProperties[alpha->index].Parent != Phase.FieldsProperties[ beta->index].Parent)
                    {
                        locMobility = 0.0;
                    }

                    if(pIndexA != pIndexB and
                       Phase.FieldsProperties[alpha->index].Parent !=  beta->index and
                       Phase.FieldsProperties[ beta->index].Parent != alpha->index)
                    {
                        locMobility = 0.0;
                    }
                }
                // Modifying mobility for nuclei (if enabled)
                if(NucleiMobilityFactor(pIndexA, pIndexB) != 1.0 and
                   (Phase.FieldsProperties[alpha->index].IsNucleus() or
                    Phase.FieldsProperties[ beta->index].IsNucleus()))
                {
                    double volume_ratio = min(Phase.FieldsProperties[alpha->index].VolumeRatio,
                                              Phase.FieldsProperties[ beta->index].VolumeRatio);

                    locMobility *= NucleiMobilityFactor(pIndexA, pIndexB)*(1.0 - volume_ratio) + 1.0*volume_ratio;
                }

                // Modifying mobility if growth constraints are violated
                if(Phase.FieldsProperties[alpha->index].GrowthConstraintsViolation != GrowthConstraintsViolations::None and
                   Phase.FieldsProperties[ beta->index].GrowthConstraintsViolation != GrowthConstraintsViolations::None)
                {
                    locMobility *= Phase.PairwiseGrowthFactors.get_sym1(alpha->index,beta->index);
                }

                // Putting interface energy and mobility values into the properties storage
                if(locEnergy != 0.0 or locMobility != 0.0)
                {
                    PropertiesDR(i,j,k).set_energy_and_mobility(alpha->index, beta->index, locEnergy, locMobility);
                }
                locMaxEnergies(pIndexA,pIndexB) = max(locMaxEnergies(pIndexA,pIndexB),locEnergy);
                locMaxMobilities(pIndexA,pIndexB) = max(locMaxMobilities(pIndexA,pIndexB),locMobility);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxEnergies = locMaxEnergies;
    maxMobilities = locMaxMobilities;

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, maxEnergies.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, maxMobilities.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::SetStoichiometry( std::vector<ThermodynamicPhase> Phase)
{
    for(size_t n = 0; n < Phase[0].Component.size(); n++)
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    for(size_t beta = 0; beta < Nphases; beta++)
    if (alpha != beta)
    {
        if (Phase[alpha].Component[n].isStoichiometric and Phase[beta].Component[n].isStoichiometric)
        {
            InterfaceMobility(alpha, beta).Mobility = 0;
            InterfaceMobility(alpha, beta).MaxMobility = 0;
        }    
    }
}
void InterfaceProperties::SetMobilityThermalEffectSR(const PhaseField& Phase, const Temperature& Tx)
{
    Matrix<double> locMaxMobilities(Nphases, Nphases);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Properties, 0, reduction(MatrixDMAX : locMaxMobilities))
    {
        if(Phase.Fields(i, j, k).wide_interface())
        {
            double invT = 1.0 / (PhysicalConstants::R * Tx(i, j, k));
            for(auto it = Properties(i, j, k).begin();
                     it != Properties(i, j, k).end(); ++it)
            {
                int pIndexA = Phase.FieldsProperties[it->indexA].Phase;
                int pIndexB = Phase.FieldsProperties[it->indexB].Phase;
                it->mobility *= exp(-InterfaceMobility(pIndexA, pIndexB).ActivationEnergy * invT);
                locMaxMobilities(pIndexA, pIndexB) = max(locMaxMobilities(pIndexA, pIndexB), it->mobility);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMobilities = locMaxMobilities;

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, maxMobilities.data(), Nphases * Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::SetMobilityThermalEffectDR(const PhaseField& Phase, const Temperature& Tx)
{
    double locMaxMobility = 0.0;
    Matrix<double> locMaxMobilities(Nphases, Nphases);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PropertiesDR, 0, reduction(max : locMaxMobility) reduction(MatrixDMAX : locMaxMobilities))
    {
        if(Phase.FieldsDR(i, j, k).wide_interface())
        {
            double invT = 1.0 / (PhysicalConstants::R * Tx.at(0.5 * i - Grid.dNx * 0.25, 0.5 * j - Grid.dNy * 0.25, 0.5 * k - Grid.dNz * 0.25));
            for(auto it = PropertiesDR(i, j, k).begin();
                     it != PropertiesDR(i, j, k).end(); ++it)
            {
                int pIndexA = Phase.FieldsProperties[it->indexA].Phase;
                int pIndexB = Phase.FieldsProperties[it->indexB].Phase;
                it->mobility *= exp(-InterfaceMobility(pIndexA, pIndexB).ActivationEnergy * invT);
                locMaxMobility                     = max(locMaxMobility, it->mobility);
                locMaxMobilities(pIndexA, pIndexB) = max(locMaxMobilities(pIndexA, pIndexB), it->mobility);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMobilities = locMaxMobilities;
#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, maxMobilities.data(), Nphases * Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if(Grid.dNx) BC.SetX(Properties);
    if(Grid.dNy) BC.SetY(Properties);
    if(Grid.dNz) BC.SetZ(Properties);
}

void InterfaceProperties::SetBoundaryConditionsDR(const BoundaryConditions& BC)
{
    if(Grid.dNx) BC.SetX(PropertiesDR);
    if(Grid.dNy) BC.SetY(PropertiesDR);
    if(Grid.dNz) BC.SetZ(PropertiesDR);
}

void InterfaceProperties::WriteVTK(const Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    const std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");
    std::string TypeName = "Energy";
    for (size_t alpha =     0; alpha < Nphases; alpha++)
    for (size_t beta  = alpha;  beta < Nphases; beta++)
    {
        if (InterfaceEnergy(alpha,beta).Model != InterfaceEnergyModels::Iso &&
            InterfaceEnergy(alpha,beta).Type == InterfaceEnergyModelType::Stiffness)
        {
            TypeName = "Stiffness";
            break;
        }
    }
    switch (Grid.Resolution)
    {
        case Resolutions::Single:
        {
            ListOfFields.push_back((VTK::Field_t)
            {TypeName, [this](int i,int j,int k)
                {
                    if (Properties(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();

                    double value = 0.0;
                    for (auto it  = Properties(i,j,k).cbegin();
                              it != Properties(i,j,k).cend(); ++it)
                    {
                        value += it->energy;
                    }
                    value /= Properties(i,j,k).size();
                    return value;
                }
            });
            ListOfFields.push_back((VTK::Field_t)
            {"Mobility", [this](int i,int j,int k)
                {
                    if (Properties(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();

                    double value = 0.0;
                    for (auto it  = Properties(i,j,k).cbegin();
                              it != Properties(i,j,k).cend(); ++it)
                    {
                        value += it->mobility;
                    }
                    value /= Properties(i,j,k).size();
                    return value;
                }
            });
            VTK::Write(Filename, locSettings, ListOfFields, precision, 1);
            break;
        }
        case Resolutions::Dual:
        {
            ListOfFields.push_back((VTK::Field_t)
            {TypeName, [this](int i,int j,int k)
                {
                    if (PropertiesDR(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();

                    double value = 0;
                    for (auto it  = PropertiesDR(i,j,k).cbegin();
                              it != PropertiesDR(i,j,k).cend(); ++it)
                    {
                        value += it->energy;
                    }
                    value /= PropertiesDR(i,j,k).size();
                    return value;
                }
            });
            ListOfFields.push_back((VTK::Field_t)
            {"Mobility", [this](int i,int j,int k)
                {
                    if (PropertiesDR(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();

                    double value = 0;
                    for (auto it  = PropertiesDR(i,j,k).cbegin();
                              it != PropertiesDR(i,j,k).cend(); ++it)
                    {
                        value += it->mobility;
                    }
                    value /= PropertiesDR(i,j,k).size();
                    return value;
                }
            });
            VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
            break;
        }
    }
}

std::set<size_t> InterfaceProperties::localPhaseIdx( const int i, const int j, const int k, const PhaseField& Phase) const
{
    auto& PF = (Grid.Resolution != Resolutions::Dual) ? Phase.Fields : Phase.FieldsDR;
    std::set<size_t> locIndex{};
    for (auto ls = Phase.LStencil.cbegin(); ls != Phase.LStencil.cend(); ls++)
    {
        const int ii = ls->di;
        const int jj = ls->dj;
        const int kk = ls->dk;

        for (auto it  = PF(i + ii, j + jj, k + kk).cbegin();
                  it != PF(i + ii, j + jj, k + kk).cend(); ++it)

        if (Phase.FieldsProperties[it->index].Volume > 1.0)
        locIndex.insert(it->index);
    }
    return locIndex;
}

bool InterfaceProperties::FirstOrderExtrapolation(const PhaseField& Phase)
{
    bool Extrapolated = true;
    auto& PF = (Grid.Resolution != Resolutions::Dual) ? Phase.Fields : Phase.FieldsDR;
    auto& IP = (Grid.Resolution != Resolutions::Dual) ? Properties : PropertiesDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,reduction(&&:Extrapolated))
    if (PF(i,j,k).wide_interface())
    {
        std::set<size_t> locIndex = localPhaseIdx(i,j,k,Phase);

        // Extrapolate curvatures for all pairs of phase fields
        for (auto alpha = locIndex.begin(); alpha != locIndex.end(); alpha++)
        for (auto beta = locIndex.begin(); beta != locIndex.end(); beta++)
        if (*beta > *alpha)
        {
            if (IP(i,j,k).get_energy(*alpha, *beta) == 0.0)
            {
                Extrapolated = Extrapolated && false;
                int count = 0;
                double sumEnergies = 0.0;
                double locE1 = 0.0;
                double locE2 = 0.0;

                // First order extrapolation and zeroth order weighting (compare to Riemann problem)
                if (Grid.dNx            ) {locE1 = IP(i+1,j,k)  .get_energy(*alpha, *beta); locE2 = IP(i+2,j,k)  .get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx            ) {locE1 = IP(i-1,j,k)  .get_energy(*alpha, *beta); locE2 = IP(i-2,j,k)  .get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNy            ) {locE1 = IP(i,j+1,k)  .get_energy(*alpha, *beta); locE2 = IP(i,j+2,k)  .get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNy            ) {locE1 = IP(i,j-1,k)  .get_energy(*alpha, *beta); locE2 = IP(i,j-2,k)  .get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNz            ) {locE1 = IP(i,j,k+1)  .get_energy(*alpha, *beta); locE2 = IP(i,j,k+2)  .get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNz            ) {locE1 = IP(i,j,k-1)  .get_energy(*alpha, *beta); locE2 = IP(i,j,k-2)  .get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNy) {locE1 = IP(i-1,j-1,k).get_energy(*alpha, *beta); locE2 = IP(i-2,j-2,k).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNy) {locE1 = IP(i+1,j-1,k).get_energy(*alpha, *beta); locE2 = IP(i+2,j-2,k).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNy) {locE1 = IP(i-1,j+1,k).get_energy(*alpha, *beta); locE2 = IP(i-2,j+2,k).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNy) {locE1 = IP(i+1,j+1,k).get_energy(*alpha, *beta); locE2 = IP(i+2,j+2,k).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNz) {locE1 = IP(i-1,j,k-1).get_energy(*alpha, *beta); locE2 = IP(i-2,j,k-2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNz) {locE1 = IP(i+1,j,k-1).get_energy(*alpha, *beta); locE2 = IP(i+2,j,k-2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNz) {locE1 = IP(i-1,j,k+1).get_energy(*alpha, *beta); locE2 = IP(i-2,j,k+2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNx && Grid.dNz) {locE1 = IP(i+1,j,k+1).get_energy(*alpha, *beta); locE2 = IP(i+2,j,k+2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNy && Grid.dNz) {locE1 = IP(i,j-1,k-1).get_energy(*alpha, *beta); locE2 = IP(i,j-2,k-2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNy && Grid.dNz) {locE1 = IP(i,j+1,k-1).get_energy(*alpha, *beta); locE2 = IP(i,j+2,k-2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNy && Grid.dNz) {locE1 = IP(i,j-1,k+1).get_energy(*alpha, *beta); locE2 = IP(i,j-2,k+2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}
                if (Grid.dNy && Grid.dNz) {locE1 = IP(i,j+1,k+1).get_energy(*alpha, *beta); locE2 = IP(i,j+2,k+2).get_energy(*alpha, *beta); if (locE1 != 0.0 && locE2 != 0.0) { sumEnergies += 2*locE1 - locE2; count++; }}

                if (count > 1)
                {
                    PropertiesExtrapolations(i,j,k).add_energy(*alpha, *beta, sumEnergies/count);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,PF,0)
    if (PF(i,j,k).wide_interface())
    {
        if (PropertiesExtrapolations(i,j,k).size() != 0u)
        {
            for (auto it = PropertiesExtrapolations(i,j,k).begin(); it < PropertiesExtrapolations(i,j,k).end(); it++)
            {
                IP(i,j,k).set_energy(it->indexA, it->indexB, it->energy);
            }
        }
    }
    STORAGE_LOOP_END

    PropertiesExtrapolations.Clear();

    return Extrapolated;
}

void InterfaceProperties::Extrapolate(const PhaseField& Phase, const BoundaryConditions& BC)
{
    bool Complete = false;
    switch (ExtrapolationMode)
    {
        case ExtrapolationModes::None:
            break;
        case ExtrapolationModes::ZerothOrder:
            break;
        case ExtrapolationModes::FirstOrder:
            auto& PF = (Grid.Resolution != Resolutions::Dual) ? Phase.Fields : Phase.FieldsDR;
            auto& IP = (Grid.Resolution != Resolutions::Dual) ? Properties : PropertiesDR;
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,IP,IP.Bcells(),)
            if (PF(i,j,k).interface_halo())
            {
                for (auto it = IP(i,j,k).begin(); it < IP(i,j,k).end(); it++) it->energy = 0.0;
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            for (size_t n = 0; n <= MaxExtrapolationInterations; n++)
            {
                Complete = FirstOrderExtrapolation(Phase);
                if (Complete) break;
            }
            SetBoundaryConditions(BC);
            break;
    }
}

}// namespace openphase

