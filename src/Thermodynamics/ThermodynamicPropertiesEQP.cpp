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
 *  File created :   2011
 *  Main contributors :   Oleg Shchyglo; Marvin Tegeler; Matthias Stratmann;
 *                        Efim Borukhovich;
 *
 */

#include "Thermodynamics/ThermodynamicPropertiesEQP.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Composition.h"
#include "DrivingForce.h"
#include "Temperature.h"
#include "BoundaryConditions.h"
#include "InterfaceProperties.h"
#include "Nucleation.h"
#include "VTK.h"

namespace openphase
{

using namespace std;

void ThermodynamicPropertiesEQP::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ThermodynamicPropertiesEQP";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;
    Ncomp = locSettings.Ncomp;

    PhaseNames = locSettings.PhaseNames;
    ElementNames = locSettings.ElementNames;

    AtStart = true;

    ExtrapolationMode = ExtrapolationModes::None;

    size_t Bcells = Grid.Bcells;

    GlobalExtrapolationData.Allocate({Nphases,Nphases});
    for(size_t m = 0; m < Nphases; m++)
    for(size_t n = 0; n < Nphases; n++)
    {
        GlobalExtrapolationData({m,n}).Allocate(Ncomp);
    }

    LocalExtrapolationData.Allocate(Grid,{Nphases,Nphases}, Bcells);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,LocalExtrapolationData,LocalExtrapolationData.Bcells(),)
    {
        for(size_t m = 0; m < Nphases; m++)
        for(size_t n = 0; n < Nphases; n++)
        {
            LocalExtrapolationData(i,j,k,{m,n}).Allocate(Ncomp);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    IntersectionTemperature.Allocate({Nphases,Nphases});
    IntersectionTemperature.set_to_zero();
    ReferenceTemperature.Allocate({Nphases,Nphases});
    ReferenceTemperature.set_to_zero();

    dMu.Allocate(Grid, {Nphases, Ncomp}, Bcells);

    EnableDrivingForce.Allocate(Nphases, Nphases);
    DrivingForceMap.Allocate(Nphases,Nphases);

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ThermodynamicPropertiesEQP::ReadInput(std::string FileName)
{
    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be opened",
                        thisclassname, "ReadInput");
        OP_Exit(EXIT_FAILURE);
    }

    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert(thisclassname);
    ConsoleOutput::WriteStandard("Source", FileName);

    ReadInput(data);

    ConsoleOutput::WriteLine();
}

void ThermodynamicPropertiesEQP::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    // Reading thermodynamic properties extrapolation modes
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
        MaxCompositionDeviation = FileInterface::ReadParameterD(inp, moduleLocation, std::vector<std::string>{"MaxCompositionDeviation","MDev"}, false, 1E-4);
        MaxTemperatureDeviation = FileInterface::ReadParameterD(inp, moduleLocation, std::vector<std::string>{"MaxTemperatureDeviation","TDev"}, false, 1E-0);
    }

    // Reading driving force models
    std::string dGmodel = FileInterface::ReadParameterK(inp, moduleLocation,std::vector<std::string>{"DrivingForceModel","DrivingForceType"}, false, "STANDARD");
    if (dGmodel == "STANDARD")
    {
        DrivingForceModel = DrivingForceModels::Standard;
    }
    else if (dGmodel == "LOWERSLOPE")
    {
        DrivingForceModel = DrivingForceModels::LowerSlope;
    }
    else if (dGmodel == "WEIGHTED")
    {
        DrivingForceModel = DrivingForceModels::Weighted;
    }
    else if (dGmodel == "USER")
    {
        DrivingForceModel = DrivingForceModels::User;
        for (size_t i = 0; i < Nphases; ++i)
        for (size_t j = 0; j < Nphases; ++j)
        {
            stringstream idx;
            idx << "_" << i << "_" << j;
            DrivingForceMap(i,j) = FileInterface::ReadParameterB(inp, moduleLocation,std::vector<std::string>{"DrivingForceContribution"+idx.str(),"DGContributionMAP"+idx.str()}, false, true);
        }
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong driving force model specified, the standard model used instead!", thisclassname, "ReadInput()");
        DrivingForceModel = DrivingForceModels::Standard;
    }

    // Reading driving force switches
    for (size_t pIndexA = 0; pIndexA < Nphases - 1; pIndexA++)
    for (size_t pIndexB = pIndexA + 1; pIndexB < Nphases; pIndexB++)
    {
        stringstream converter;
        converter << "EnableDrivingForce_" << pIndexA << "_" << pIndexB;
        EnableDrivingForce(pIndexA, pIndexB) = FileInterface::ReadParameterB(inp, moduleLocation, converter.str(), false, true);
        EnableDrivingForce(pIndexB, pIndexA) = EnableDrivingForce(pIndexA, pIndexB);
    }
}

void ThermodynamicPropertiesEQP::PrintStatistics(Composition& Cx, bool matrix,
                             bool warnings, double tStep)
{

}

void ThermodynamicPropertiesEQP::CheckRange(PhaseField& Phi,
                                            Composition& Cx,
                                            Temperature& Tx)
{
    switch(ExtrapolationMode)
    {
        case ExtrapolationModes::Global:
        {
            for(size_t m = 0; m < Nphases; m++)
            for(size_t n = 0; n < Nphases; n++)
            if(m != n and fabs(Tx.Tavg - GlobalExtrapolationData({m,n}).eq_temperature) > MaxTemperatureDeviation)
            {
                GlobalExtrapolationData({m,n}).out_of_T_range = true;
            }
            break;
        }
        case ExtrapolationModes::Local:
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
            {
                if(Phi.Fields(i,j,k).wide_interface())
                {
                    for(size_t m = 0; m < Nphases; m++)
                    for(size_t n = 0; n < Nphases; n++)
                    if(m != n and fabs(Tx.Tx(i,j,k) - LocalExtrapolationData(i,j,k,{m,n}).eq_temperature) > MaxTemperatureDeviation)
                    {
                        LocalExtrapolationData(i,j,k,{m,n}).out_of_T_range = true;
                    }
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
                if(Phi.Fields(i,j,k).wide_interface())
                {
                    for(size_t m = 0; m < Nphases; m++)
                    for(size_t n = 0; n < Nphases; n++)
                    if(m != n)
                    {
                        LocalExtrapolationData(i,j,k,{m,n}).out_of_T_range = true;
                    }
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
    }
}

void ThermodynamicPropertiesEQP::CalculateDrivingForce(PhaseField& Phi,
                                  InterfaceProperties& IP, Composition& Cx,
                                  Temperature& Tx, DrivingForce& dGab)
{
    /** This function calculates the driving force for a phase pair in each
    point in the interface. */
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
    {
        if(Phi.Fields(i,j,k).interface())
        {
            const Tensor<EquilibriumData,2>& locExtrapolationData = ExtrapolationData(i,j,k);

            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phi.Fields(i,j,k).cend(); ++beta)
            {
                size_t pIndexA = Phi.FieldsProperties[alpha->index].Phase;
                size_t pIndexB = Phi.FieldsProperties[ beta->index].Phase;

                if(EnableDrivingForce(pIndexA, pIndexB))
                {
                    if((pIndexA != pIndexB) and
                       (locExtrapolationData({pIndexA,pIndexB}).is_set and
                        locExtrapolationData({pIndexB,pIndexA}).is_set))
                    {
                        dVectorN vCx_A(Cx.Ncomp);
                        dVectorN vCx_B(Cx.Ncomp);
                        dVectorN vEC_A(Cx.Ncomp);
                        dVectorN vEC_B(Cx.Ncomp);
                        dVectorN vmAB(Cx.Ncomp);
                        dVectorN vmBA(Cx.Ncomp);

                        for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                        {
                            vEC_A[comp] = locExtrapolationData({pIndexA,pIndexB}).eq_composition_extra(comp, Tx.Tx(i,j,k));
                            vEC_B[comp] = locExtrapolationData({pIndexB,pIndexA}).eq_composition_extra(comp, Tx.Tx(i,j,k));

                            vCx_A[comp] = Cx.MoleFractions(i,j,k,{pIndexA,comp});
                            vCx_B[comp] = Cx.MoleFractions(i,j,k,{pIndexB,comp});

                            if (locExtrapolationData({pIndexA,pIndexB}).eq_slope({comp}) != 0.0)
                            {
                                vmAB[comp] = 1.0/locExtrapolationData({pIndexA,pIndexB}).eq_slope({comp});
                            }
                            else
                            {
                                vmAB[comp] = 0.0;
                            }

                            if (locExtrapolationData({pIndexB,pIndexA}).eq_slope({comp}) != 0.0)
                            {
                                vmBA[comp] = 1.0/locExtrapolationData({pIndexB,pIndexA}).eq_slope({comp});
                            }
                            else
                            {
                                vmBA[comp] = 0.0;
                            }
                        }

                        dVectorN vdCA = vCx_A - vEC_A;                                   // Actual composition minus equilibrium composition (from last equilibrium)
                        dVectorN vdCB = vCx_B - vEC_B;                                   // Actual composition minus equilibrium composition (from last equilibrium)
                        dVectorN TIELINE = vEC_A - vEC_B;

                        TIELINE.normalize();

                        double projA = vdCA.dot(TIELINE);
                        double projB = vdCB.dot(TIELINE);

                        double lambdaAB = vmAB.dot(TIELINE);
                        double lambdaBA = vmBA.dot(TIELINE);

                        double dE = (locExtrapolationData({pIndexB,pIndexA}).eq_entropy -
                                     locExtrapolationData({pIndexA,pIndexB}).eq_entropy)
                                     /Cx.TotalMolarVolume;

                        double locdG_AB = 0.0;
                        int counter = 0;
                        if (lambdaAB != 0.0)
                        {
                            locdG_AB = projA*std::fabs(dE/lambdaAB);
                            counter++;
                        }
                        double locdG_BA = 0.0;
                        if (lambdaBA != 0.0)
                        {
                            locdG_BA = projB*std::fabs(dE/lambdaBA);
                            counter++;
                        }

                        double dG_AB = 0.0;

                        switch(DrivingForceModel)
                        {
                            case DrivingForceModels::Standard:
                            {
                                dG_AB = (locdG_AB + locdG_BA)/counter;
                                break;
                            }
                            case DrivingForceModels::LowerSlope:
                            {
                                if (std::fabs(lambdaAB) > std::fabs(lambdaBA))
                                {
                                    dG_AB = locdG_AB;
                                }
                                else
                                {
                                    dG_AB = locdG_BA;
                                }
                                break;
                            }
                            case DrivingForceModels::Weighted:
                            {
                                double invlambdaAB = (lambdaAB != 0.0) ? 1.0/lambdaAB : 0.0;
                                double invlambdaBA = (lambdaBA != 0.0) ? 1.0/lambdaBA : 0.0;

                                double total = std::fabs(invlambdaAB) +
                                               std::fabs(invlambdaBA);

                                dG_AB = (std::fabs(invlambdaAB)*locdG_AB +
                                         std::fabs(invlambdaBA)*locdG_BA)/total;
                                break;
                            }
                            case DrivingForceModels::User:
                            {
                                counter = 0;
                                if (DrivingForceMap(pIndexA,pIndexB))
                                {
                                    dG_AB += locdG_AB;
                                    counter++;
                                }
                                if (DrivingForceMap(pIndexB,pIndexA))
                                {
                                    dG_AB += locdG_BA;
                                    counter++;
                                }
                                if(counter)
                                {
                                    dG_AB /= counter;
                                }
                                break;
                            }
                        }

                        dGab.Force(i,j,k).add_raw(alpha->index, beta->index, dG_AB);
                    }
                    else if(pIndexA != pIndexB)
                    {
                        double Prefactor = Pi/Phi.Grid.Eta;
                        double maxdG = dGab.Limit(pIndexA,pIndexB)
                                      *IP.InterfaceEnergy(pIndexA,pIndexB).MaxEnergy
                                      *Prefactor;                                       ///< Maximum driving force allowed with the given settings.

                        if(locExtrapolationData({pIndexA,pIndexB}).eq_fraction >
                           locExtrapolationData({pIndexB,pIndexA}).eq_fraction)
                        {
                            dGab.Force(i,j,k).add_raw(alpha->index,beta->index, maxdG);
                        }
                        else
                        {
                            dGab.Force(i,j,k).add_raw(alpha->index,beta->index,-maxdG);
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ThermodynamicPropertiesEQP::SetBoundaryConditions(const BoundaryConditions& BC)
{

}

void ThermodynamicPropertiesEQP::Remesh(int newNx, int newNy, int newNz,
                                        const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);

    LocalExtrapolationData.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    dMu.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}
void ThermodynamicPropertiesEQP::ClearChemicalPotential(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,dMu,dMu.Bcells(),)
    {
        for(size_t alpha = 0; alpha != Nphases; ++alpha)
        for(size_t  comp = 0;  comp != Ncomp;   ++comp)

        {
            dMu(i,j,k,{alpha, comp}) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ThermodynamicPropertiesEQP::WriteVTK(const Settings& locSettings,
                                              const int tStep,
                                              const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;

    for(size_t alpha = 0; alpha < Nphases; alpha++)
    for(size_t  beta = 0;  beta < Nphases; beta++)
    if (alpha != beta)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            const std::string index1 = "_" + locSettings.PhaseNames[alpha]
                                     + "_" + locSettings.PhaseNames[beta]
                                     + "_" + locSettings.ElementNames[comp];
            const std::string name1 = "EquilibriumSlope" + index1;
            ListOfFields.push_back((VTK::Field_t) {name1, [this, alpha, beta, comp](int i,int j,int k){return ExtrapolationData(i,j,k)({alpha,beta}).eq_slope({comp});}});
            const std::string name2 = "EquilibriumComp" + index1;
            ListOfFields.push_back((VTK::Field_t) {name1, [this, alpha, beta, comp](int i,int j,int k){return ExtrapolationData(i,j,k)({alpha,beta}).eq_composition({comp});}});
        }
        const std::string index2 = "_" + locSettings.PhaseNames[alpha]
                                 + "_" + locSettings.PhaseNames[beta];
        const std::string name1 = "EquilibriumSlope" + index2;
        ListOfFields.push_back((VTK::Field_t) {name1, [this, alpha, beta](int i,int j,int k){return ExtrapolationData(i,j,k)({alpha,beta}).eq_entropy - ExtrapolationData(i,j,k)({beta,alpha}).eq_entropy;}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + '_', tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void ThermodynamicPropertiesEQP::PrintPointStatistics(int i, int j, int k)
{
    for(size_t alpha = 0; alpha < Nphases-1; alpha++)
    for(size_t beta = alpha+1; beta < Nphases; beta++)
    {
        cout << "Entropy difference between phases "
             << alpha << " and " << beta << "     : "
             << ExtrapolationData(i,j,k)({alpha,beta}).eq_entropy -
                ExtrapolationData(i,j,k)({beta,alpha}).eq_entropy << endl;
        cout << endl;
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            cout << "Equilibrium phase composition "
                 << ElementNames[comp] << " (in phase pair "
                 << alpha << "_" << beta << ")    : "
                 << ExtrapolationData(i,j,k)({alpha,beta}).eq_composition({comp}) << endl;
        }
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            cout << "Equilibrium phase composition "
                 << ElementNames[comp] << " (in phase pair "
                 << beta << "_" << alpha << ")    : "
                 << ExtrapolationData(i,j,k)({beta,alpha}).eq_composition({comp}) << endl;
        }
        cout << endl;
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            cout << "Equilibrium phase-diagram slope "
                 << ElementNames[comp] << " (in phase pair "
                 << alpha << "_" << beta << ")    : "
                 << ExtrapolationData(i,j,k)({alpha,beta}).eq_slope({comp}) << endl;
        }
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            cout << "Equilibrium phase-diagram slope "
                 << ElementNames[comp] << " (in phase pair "
                 << beta << "_" << alpha << ")    : "
                 << ExtrapolationData(i,j,k)({beta,alpha}).eq_slope({comp}) << endl;
        }
    }
}

}//namespace openphase
