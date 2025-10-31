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

 *   File created :   2025
 *   Main contributors :   Oleg Shchyglo; Reza Namdar
 *
 */ 
#include "ReactiveFlows/SpeciesTransport.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <limits>

#include <vector>
#include <random>
#include <iostream>

using namespace std;
using namespace openphase;

void SpeciesTransport::ReadInput(string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("SpeciesTransport");
    ConsoleOutput::WriteStandard("Source", InputFile);
    std::fstream inp(InputFile, std::ios::in);
    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << InputFile << "\" could not be opened";
        throw std::runtime_error(message.str());
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();
    int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "SpeciesTransport");
    BCOrder    	         = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("BCOrder"), false, 0);
    TempBurntGas		 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TempBurntGas"), false, 0.0);
    Species		    	 = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("Species"), false, false);
    FuelIndex 	         = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("FuelIndex"), false, 0);
    AdvectionUpwind    	 = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("AdvectionUpwind"), false, false);
    AdvectionCentral     = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("AdvectionCentral"), false, false);
    AdvectionVanLeer     = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("AdvectionVanLeer"), false, true);
    X0BurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("X0BurntZone"), false, -10);
    XNBurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("XNBurntZone"), false, -10);
    Y0BurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Y0BurntZone"), false, -10);
    YNBurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("YNBurntZone"), false, -10);
    Z0BurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Z0BurntZone"), false, -10);
    ZNBurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ZNBurntZone"), false, -10);
    X0FuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("X0FuelZone"), false, -10);
    XNFuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("XNFuelZone"), false, -10);
    Y0FuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Y0FuelZone"), false, -10);
    YNFuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("YNFuelZone"), false, -10);
    Z0FuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Z0FuelZone"), false, -10);
    ZNFuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ZNFuelZone"), false, -10);
    IF_REACTION          = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("IF_REACTION"), false, false);
    nCoeffs    	         = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("nCoeffs"), false, 15);
}

void SpeciesTransport::Initialize(Settings& locSettings)
{
    Grid = locSettings.Grid;
    size_t Bcells = locSettings.Grid.Bcells;

    nSpecies = locSettings.Ncomp;
    HRR.Allocate (locSettings.Grid, Bcells);
    Cp_Species.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    MassDiff_Species.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    W_Species.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    h_Species.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    MassFractions.Allocate(locSettings.Grid,{nSpecies}, Bcells);
    MassFractionsOld.Allocate(locSettings.Grid,{nSpecies}, Bcells);    
    MolecularWeight.Allocate({nSpecies});
    Y0X.Allocate({nSpecies});
    PolyNomCoeffs.Allocate({nSpecies, nCoeffs});
}

void SpeciesTransport::UpdateFields(EnergyTransport& ET)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN (i,j,k,MassFractions,MassFractions.Bcells(),)
    {
        ET.TxOld(i,j,k)         = ET.Tx(i,j,k);
        MassFractionsOld(i,j,k) = MassFractions(i,j,k);
        if(ET.Conjugate)
        {
            ET.TsOld(i,j,k) = ET.Ts(i,j,k);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SpeciesTransport::CalculateSpeciesHeatCapacitiesAndEnthalpies(EnergyTransport& ET, FlowSolverLBM& FL)
{
	double R = 8.314510;   // J/mol.K
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cp_Species,Cp_Species.Bcells(),)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			double T = ET.Tx(i,j,k);
			for (size_t ic = 0 ; ic < nSpecies; ic++)
    		{
				if(T < PolyNomCoeffs({ic,0}))
				{
					Cp_Species(i,j,k,{ic}) = (PolyNomCoeffs({ic,8}) + PolyNomCoeffs({ic,9})*T + PolyNomCoeffs({ic,10})*pow(T,2.0) 
                                            + PolyNomCoeffs({ic,11})*pow(T,3.0) + PolyNomCoeffs({ic,12})*pow(T,4.0) ) 
                                            * R / MolecularWeight({ic});
					h_Species(i,j,k,{ic})  = (PolyNomCoeffs({ic,8}) + PolyNomCoeffs({ic,9})/2.0*T + PolyNomCoeffs({ic,10})/3.0*pow(T,2.0)
                                            + PolyNomCoeffs({ic,11})/4.0*pow(T,3.0) + PolyNomCoeffs({ic,12})/5.0*pow(T,4.0) + PolyNomCoeffs({ic,13})/T) 
                                            * R * T / MolecularWeight({ic});
				}
				else
				{
					Cp_Species(i,j,k,{ic}) = (PolyNomCoeffs({ic,1}) + PolyNomCoeffs({ic,2})*T + PolyNomCoeffs({ic,3})*pow(T,2.0) 
                                            + PolyNomCoeffs({ic,4})*pow(T,3.0) + PolyNomCoeffs({ic,5})*pow(T,4.0) ) 
                                            * R / MolecularWeight({ic});
					h_Species(i,j,k,{ic})  = (PolyNomCoeffs({ic,1}) + PolyNomCoeffs({ic,2})/2.0*T + PolyNomCoeffs({ic,3})/3.0*pow(T,2.0) 
                                            + PolyNomCoeffs({ic,4})/4.0*pow(T,3.0) + PolyNomCoeffs({ic,5})/5.0*pow(T,4.0) + PolyNomCoeffs({ic,6})/T)
                                            * R * T / MolecularWeight({ic});
				}
				
    		}
		}
	}
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SpeciesTransport::SetFreeBCNX(PhaseField& Phase, EnergyTransport& ET)
{
    int iEnd    = Grid.Nx;
    double Temp = ET.T0; 		// Temporary variable for extrapolated temperature
    double f4[nSpecies];    // Temporary array for extrapolated compositions
    double f3[nSpecies];
    double f2[nSpecies];
    double f1[nSpecies];
    // Loop over the extrapolated boundary nodes
    if(Grid.OffsetX+Grid.Nx==Grid.TotalNx)
    for (int i = iEnd; i < Grid.Nx + MassFractions.Bcells(); ++i)
    {
        for (int j = 0; j< Grid.Ny ; ++j)
        {
            for (int k = 0; k< Grid.Nz ; ++k)
            {
                for (size_t ic = 0; ic < nSpecies; ic++)
                {
                    // Load values from previous points for extrapolation
                    f4[ic] = MassFractions(iEnd - 1, j, k,{ic});
                    f3[ic] = MassFractions(iEnd - 2, j, k,{ic});
                    f2[ic] = MassFractions(iEnd - 3, j, k,{ic});
                    f1[ic] = MassFractions(iEnd - 4, j, k,{ic});
                }

                if (BCOrder==0)
                {
                    Temp = ET.Tx(iEnd - 1, j, k);
                    for (size_t ic = 0; ic < nSpecies; ic++)
                    {
                        MassFractions(i, j, k,{ic}) = f4[ic];
                    }
                }
                else if (BCOrder==1)
                {
                    Temp = 2.0 * ET.Tx(iEnd - 1, 0, k) - ET.Tx(iEnd - 2, 0, k);
                    for (size_t ic = 0; ic < nSpecies; ic++)
                    {
                        MassFractions(i, j, k,{ic}) = 2.0 * f4[ic] - f3[ic];
                    }
                }
                else if (BCOrder==2)
                {
                    Temp = ET.Tx(iEnd - 1, j, k) + (ET.Tx(iEnd - 1, j, k) - ET.Tx(iEnd - 2, j, k)) +
                           0.5 * (ET.Tx(iEnd - 1, j, k) - 2.0 * ET.Tx(iEnd - 2, j, k) + ET.Tx(iEnd - 3, j, k));
                    for (size_t ic = 0; ic < nSpecies; ic++)
                    {
                        MassFractions(i, j, k,{ic}) = f4[ic] + (f4[ic] - f3[ic]) +
                                                              0.5 * (f4[ic] - 2.0 * f3[ic] + f2[ic]);
                    }
                }
                else if (BCOrder==3)
                {
                    Temp = ET.Lagrange_polynomial4(0.0, 1.0, 2.0, 3.0, 4.0,
                                                ET.Tx(iEnd - 4, j, k), ET.Tx(iEnd - 3, j, k),
                                                ET.Tx(iEnd - 2, j, k), ET.Tx(iEnd - 1, j, k));
                    for (size_t ic = 0; ic < nSpecies; ic++)
                    {
                        MassFractions(i, j, k,{ic}) = ET.Lagrange_polynomial4(0.0, 1.0, 2.0, 3.0, 4.0, f1[ic], f2[ic], f3[ic], f4[ic]);
                    }
                }
                ET.Tx(i, j, k) = Temp;
            }
        }
    }
}

void SpeciesTransport::SetBoundaryConditions(EnergyTransport& ET, BoundaryConditions& BC)
{
    ET.SetBoundaryConditions(BC);
    if (Grid.dNx > 0) BC.SetX(MassFractions);
    if (Grid.dNy > 0) BC.SetY(MassFractions);
    if (Grid.dNz > 0) BC.SetZ(MassFractions);
}

double SpeciesTransport::CalculateMeanMolarMass(const int i, const int j, const int k, string whichtime)
{
    double MMWi=0.0;
    for(size_t n = 0; n < nSpecies; n++)
    {
        if (whichtime == "old")
        {
            MMWi += MassFractionsOld(i,j,k,{n}) / MolecularWeight({n});
        }
        if (whichtime == "new")
        {
            MMWi += MassFractions(i,j,k,{n}) / MolecularWeight({n});
        }
    }
    double MMW=1.0/MMWi;
    return MMW;
}

double SpeciesTransport::MoleFraction( const int i, const int j, const int k, double MMW, size_t kc)
{
    return ( MMW/MolecularWeight({kc}) *  MassFractionsOld(i,j,k,{kc}) );
}

void SpeciesTransport::UpdateGhostPoints(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, SolidBody& SB)
{
    if(ET.Conjugate) CalculateGhostPointsConjugateHeatTransfer(Phase,ET,FL,SB);
    if(!ET.Conjugate) CalculateGhostPoints(Phase,ET,FL,SB);
}

void SpeciesTransport::CalculateGhostPointsConjugateHeatTransfer(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, SolidBody& SB)
{
    size_t Gasidx = 0;
    size_t Solididx = 1;
    double D = SB.nDist;
    double ep=1e-8;
    
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].State == AggregateStates::Gas)
    {
        Gasidx = idx;
    }
    int iBC=Grid.dNx*2.0;
    int jBC=Grid.dNy*2.0;
    int kBC=Grid.dNz*2.0;

    int idi=Grid.dNx;
    int jdi=Grid.dNy;
    int kdi=Grid.dNz;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ET.Tx,0,)
    {
        bool Set_TsOut=false;
        if(!FL.Obstacle(i,j,k))
        for(int ii = -iBC; ii <= iBC; ii++)
        for(int jj = -jBC; jj <= jBC; jj++)
        for(int kk = -kBC; kk <= kBC; kk++)
        {
            if(FL.Obstacle(i+ii,j+jj,k+kk))
            {
                if(!Set_TsOut)
                {
                    dVector3 Norm{0.0,0.0,0.0};
                    for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
                    {
                        if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
                        {
                            Solididx=it->index;
                            NodeAB<dVector3,dVector3> Normals=Phase.Normals(i,j,k);
                            Norm=Normals.get_asym1(Solididx, Gasidx);
                        }
                    }
                    double Phif = Phase.Fractions(i,j,k,{Gasidx});
                    double X    = fabs(Phase.Grid.iWidth/Pi * asin(1.0-2.0*Phif));
                    dVector3 Xf = {double(i), double(j), double(k)};
                    dVector3 Xb = Xf - Norm*X;
                    dVector3 Xs = Xb - Norm*D;
                    double ic   = Xs.getX();
                    double jc   = Xs.getY();
                    double kc   = Xs.getZ();
                    double i0   = floor(ic);
                    double i1   = i0 + 1.0;
                    double j0   = floor(jc);
                    double j1   = j0 + 1.0;
                    double k0   = floor(kc);
                    double k1   = k0 + 1.0;

                    double Tsi  =  SB.TrilinearInterPolation(ic,jc,kc,
                                                          ET.Ts(i0*idi,j0*jdi,k0*kdi),ET.Ts(i1*idi,j0*jdi,k0*kdi),
                                                          ET.Ts(i0*idi,j1*jdi,k0*kdi),ET.Ts(i1*idi,j1*jdi,k0*kdi),
                                                          ET.Ts(i0*idi,j0*jdi,k1*kdi),ET.Ts(i1*idi,j0*jdi,k1*kdi),
                                                          ET.Ts(i0*idi,j1*jdi,k1*kdi),ET.Ts(i1*idi,j1*jdi,k1*kdi));
                    double Ksi  =  SB.TrilinearInterPolation(ic,jc,kc,
                                                          ET.K_Solid(i0*idi,j0*jdi,k0*kdi),ET.K_Solid(i1*idi,j0*jdi,k0*kdi),
                                                          ET.K_Solid(i0*idi,j1*jdi,k0*kdi),ET.K_Solid(i1*idi,j1*jdi,k0*kdi),
                                                          ET.K_Solid(i0*idi,j0*jdi,k1*kdi),ET.K_Solid(i1*idi,j0*jdi,k1*kdi),
                                                          ET.K_Solid(i0*idi,j1*jdi,k1*kdi),ET.K_Solid(i1*idi,j1*jdi,k1*kdi));
                    double Tf  = ET.Tx(i,j,k);
                    double Kf  = ET.K_Mixture(i,j,k);

                    double Tb  = (Kf/X*Tf+Ksi/D*Tsi)/(Kf/X+Ksi/D);

                    ET.Ts (i,j,k) = Tb * (1.0+X/D) - X/D*Tsi;
                    Set_TsOut=true;
                }
            }
        }

        bool Set_TfOut=false;
        if(FL.Obstacle(i,j,k))
        for(int ii = -iBC; ii <= iBC; ii++)
        for(int jj = -jBC; jj <= jBC; jj++)
        for(int kk = -kBC; kk <= kBC; kk++)
        {
            if(!FL.Obstacle(i+ii,j+jj,k+kk))
            {
                if(!Set_TfOut)
                {
                    dVector3 Norm{0.0,0.0,0.0};
                    for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
                    {
                        if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
                        {
                            Solididx=it->index;
                            NodeAB<dVector3,dVector3> Normals=Phase.Normals(i,j,k);
                            Norm=Normals.get_asym1(Solididx, Gasidx);
                        }
                    }
                    double Phif = Phase.Fractions(i,j,k,{Gasidx});
                    double X    = fabs(Phase.Grid.iWidth/Pi * asin(1.0-2.0*Phif));
                    dVector3 Xs = {double(i), double(j), double(k)};
                    dVector3 Xb = Xs + Norm*X;
                    dVector3 Xf = Xb + Norm*D;
                    double ic   = Xf.getX();
                    double jc   = Xf.getY();
                    double kc   = Xf.getZ();
                    double i0   = floor(ic);
                    double i1   = i0 + 1.0;
                    double j0   = floor(jc);
                    double j1   = j0 + 1.0;
                    double k0   = floor(kc);
                    double k1   = k0 + 1.0;

                    double Tfi  = SB.TrilinearInterPolation(ic,jc,kc,
                                                        ET.Tx(i0*idi,j0*jdi,k0*kdi),ET.Tx(i1*idi,j0*jdi,k0*kdi),
                                                        ET.Tx(i0*idi,j1*jdi,k0*kdi),ET.Tx(i1*idi,j1*jdi,k0*kdi),
                                                        ET.Tx(i0*idi,j0*jdi,k1*kdi),ET.Tx(i1*idi,j0*jdi,k1*kdi),
                                                        ET.Tx(i0*idi,j1*jdi,k1*kdi),ET.Tx(i1*idi,j1*jdi,k1*kdi));
                    double Kfi  = SB.TrilinearInterPolation(ic,jc,kc,
                                                        ET.K_Mixture(i0*idi,j0*jdi,k0*kdi),ET.K_Mixture(i1*idi,j0*jdi,k0*kdi),
                                                        ET.K_Mixture(i0*idi,j1*jdi,k0*kdi),ET.K_Mixture(i1*idi,j1*jdi,k0*kdi),
                                                        ET.K_Mixture(i0*idi,j0*jdi,k1*kdi),ET.K_Mixture(i1*idi,j0*jdi,k1*kdi),
                                                        ET.K_Mixture(i0*idi,j1*jdi,k1*kdi),ET.K_Mixture(i1*idi,j1*jdi,k1*kdi));
                    double Tsi  = ET.Ts(i,j,k);
                    double Ksi  = ET.K_Solid(i,j,k);

                    if(X<ep)
                    {
                        ET.Tx (i,j,k) = Tsi;
                    }
                    else
                    {
                        double Tb   = (Kfi/D*Tfi+Ksi/X*Tsi)/(Kfi/D+Ksi/X);
                        ET.Tx(i,j,k) =  Tb * (1.0+X/D) - X/D * Tfi;
                    }

                    double Cf[nSpecies];
                    for (size_t isp = 0; isp < nSpecies; isp++)
                    {
                        Cf[isp]= SB.TrilinearInterPolation(ic,jc,kc,
                                                        MassFractions(i0*idi,j0*jdi,k0*kdi,{isp}),MassFractions(i1*idi,j0*jdi,k0*kdi,{isp}),
                                                        MassFractions(i0*idi,j1*jdi,k0*kdi,{isp}),MassFractions(i1*idi,j1*jdi,k0*kdi,{isp}),
                                                        MassFractions(i0*idi,j0*jdi,k1*kdi,{isp}),MassFractions(i1*idi,j0*jdi,k1*kdi,{isp}),
                                                        MassFractions(i0*idi,j1*jdi,k1*kdi,{isp}),MassFractions(i1*idi,j1*jdi,k1*kdi,{isp}));
                    }
                    for (size_t isp=0; isp<nSpecies; isp++)
                    {
                        MassFractions(i,j,k,{isp}) = Cf[isp] ;
                    }
                    Set_TfOut=true;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SpeciesTransport::CalculateGhostPoints(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, SolidBody& SB)
{
    size_t Gasidx = 0;
    size_t Solididx = 1;
    double ep=1e-8;
    double Tb=ET.TempSolid;
    double Qb=ET.HeatFlux;
    double D = SB.nDist;
    
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].State == AggregateStates::Gas)
    {
        Gasidx = idx;
    }
    int iBC=Grid.dNx*2.0;
    int jBC=Grid.dNy*2.0;
    int kBC=Grid.dNz*2.0;

    int idi=Grid.dNx;
    int jdi=Grid.dNy;
    int kdi=Grid.dNz;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ET.Tx,0,)
    {
        bool Set_Tfout=false;
        if(FL.Obstacle(i,j,k))
        for(int ii = -iBC; ii <= iBC; ii++)
        for(int jj = -jBC; jj <= jBC; jj++)
        for(int kk = -kBC; kk <= kBC; kk++)
        {
            if(!FL.Obstacle(i+ii,j+jj,k+kk))
            {
                if(!Set_Tfout)
                {
                    dVector3 Norm{0.0,0.0,0.0};
                    for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
                    {
                        if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
                        {
                            Solididx=it->index;
                            NodeAB<dVector3,dVector3> Normals=Phase.Normals(i,j,k);
                            Norm=Normals.get_asym1(Solididx, Gasidx);
                        }
                    }
                    double Phif = Phase.Fractions(i,j,k,{Gasidx});
                    double X    = fabs(Phase.Grid.iWidth/Pi * asin(1.0-2.0*Phif));
                    dVector3 Xs = {double(i), double(j), double(k)};
                    dVector3 Xb = Xs + Norm*X;
                    dVector3 Xf = Xb + Norm*D;
                    double ic   = Xf.getX();
                    double jc   = Xf.getY();
                    double kc   = Xf.getZ();
                    double i0   = floor(ic);
                    double i1   = i0 + 1.0;
                    double j0   = floor(jc);
                    double j1   = j0 + 1.0;
                    double k0   = floor(kc);
                    double k1   = k0 + 1.0;
                    
                    double Tf   =  SB.TrilinearInterPolation(ic,jc,kc,
                                                          ET.Tx(i0*idi,j0*jdi,k0*kdi),ET.Tx(i1*idi,j0*jdi,k0*kdi),
                                                          ET.Tx(i0*idi,j1*jdi,k0*kdi),ET.Tx(i1*idi,j1*jdi,k0*kdi),
                                                          ET.Tx(i0*idi,j0*jdi,k1*kdi),ET.Tx(i1*idi,j0*jdi,k1*kdi),
                                                          ET.Tx(i0*idi,j1*jdi,k1*kdi),ET.Tx(i1*idi,j1*jdi,k1*kdi));
                    double Kf   =  SB.TrilinearInterPolation(ic,jc,kc,
                                                          ET.K_Mixture(i0*idi,j0*jdi,k0*kdi),ET.K_Mixture(i1*idi,j0*jdi,k0*kdi),
                                                          ET.K_Mixture(i0*idi,j1*jdi,k0*kdi),ET.K_Mixture(i1*idi,j1*jdi,k0*kdi),
                                                          ET.K_Mixture(i0*idi,j0*jdi,k1*kdi),ET.K_Mixture(i1*idi,j0*jdi,k1*kdi),
                                                          ET.K_Mixture(i0*idi,j1*jdi,k1*kdi),ET.K_Mixture(i1*idi,j1*jdi,k1*kdi));
                    double Cf[nSpecies];
                    for (size_t isp = 0; isp < nSpecies; isp++)
                    {
                        Cf[isp]= SB.TrilinearInterPolation(ic,jc,kc,
                                                        MassFractions(i0*idi,j0*jdi,k0*kdi,{isp}),MassFractions(i1*idi,j0*jdi,k0*kdi,{isp}),
                                                        MassFractions(i0*idi,j1*jdi,k0*kdi,{isp}),MassFractions(i1*idi,j1*jdi,k0*kdi,{isp}),
                                                        MassFractions(i0*idi,j0*jdi,k1*kdi,{isp}),MassFractions(i1*idi,j0*jdi,k1*kdi,{isp}),
                                                        MassFractions(i0*idi,j1*jdi,k1*kdi,{isp}),MassFractions(i1*idi,j1*jdi,k1*kdi,{isp}));
                    }

                    //dVector3 Xf_2 = Xb + Norm*D*2.0;
                    //double ic_2   = Xf_2.getX();
                    //double jc_2   = Xf_2.getY();
                    //double kc_2   = Xf_2.getZ();
                    //double i0_2   = floor(ic_2);
                    //double i1_2   = i0_2 + 1.0;
                    //double j0_2   = floor(jc_2);
                    //double j1_2   = j0_2 + 1.0;
                    //double k0_2   = floor(kc_2);
                    //double k1_2   = k0_2 + 1.0;

//                    double Tf_2   =  SB.TrilinearInterPolation(ic_2,jc_2,kc_2,
//                                                          ET.Tx(i0_2*idi,j0_2*jdi,k0_2*kdi),ET.Tx(i1_2*idi,j0_2*jdi,k0_2*kdi),
//                                                          ET.Tx(i0_2*idi,j1_2*jdi,k0_2*kdi),ET.Tx(i1_2*idi,j1_2*jdi,k0_2*kdi),
//                                                          ET.Tx(i0_2*idi,j0_2*jdi,k1_2*kdi),ET.Tx(i1_2*idi,j0_2*jdi,k1_2*kdi),
//                                                          ET.Tx(i0_2*idi,j1_2*jdi,k1_2*kdi),ET.Tx(i1_2*idi,j1_2*jdi,k1_2*kdi));
//                    double Kf_2   =  SB.TrilinearInterPolation(ic_2,jc_2,kc_2,
//                                                          ET.K_Mixture(i0_2*idi,j0_2*jdi,k0_2*kdi),ET.K_Mixture(i1_2*idi,j0_2*jdi,k0_2*kdi),
//                                                          ET.K_Mixture(i0_2*idi,j1_2*jdi,k0_2*kdi),ET.K_Mixture(i1_2*idi,j1_2*jdi,k0_2*kdi),
//                                                          ET.K_Mixture(i0_2*idi,j0_2*jdi,k1_2*kdi),ET.K_Mixture(i1_2*idi,j0_2*jdi,k1_2*kdi),
//                                                          ET.K_Mixture(i0_2*idi,j1_2*jdi,k1_2*kdi),ET.K_Mixture(i1_2*idi,j1_2*jdi,k1_2*kdi));
//                    double Cf_2[nSpecies];
//                    for (size_t isp = 0; isp < nSpecies; isp++)
//                    {
//                        Cf_2[isp]= SB.TrilinearInterPolation(ic_2,jc_2,kc_2,
//                                                        MassFractions(i0_2*idi,j0_2*jdi,k0_2*kdi,{isp}),MassFractions(i1_2*idi,j0_2*jdi,k0_2*kdi,{isp}),
//                                                        MassFractions(i0_2*idi,j1_2*jdi,k0_2*kdi,{isp}),MassFractions(i1_2*idi,j1_2*jdi,k0_2*kdi,{isp}),
//                                                        MassFractions(i0_2*idi,j0_2*jdi,k1_2*kdi,{isp}),MassFractions(i1_2*idi,j0_2*jdi,k1_2*kdi,{isp}),
//                                                        MassFractions(i0_2*idi,j1_2*jdi,k1_2*kdi,{isp}),MassFractions(i1_2*idi,j1_2*jdi,k1_2*kdi,{isp}));
//                    }

                    if(ET.IF_ConstTemp[Solididx])
                    {
                        Tb = ET.SurfaceTemp[Solididx];
                        if(X<ep)
                        {
                            ET.Tx(i,j,k) = Tb;
                        }
                        else
                        {
                            ET.Tx(i,j,k) =  Tb * (1.0+X/D) - X/D * Tf;
                        }

                    }
                    else if(ET.IF_ConstFlux[Solididx])
                    {
                        Qb = ET.SurfaceFlux[Solididx];
                        ET.Tx(i,j,k) = Tf-(X+D)*Qb/Kf;
                    }
                    for (size_t isp=0; isp<nSpecies; isp++)
                    {
                        MassFractions(i,j,k,{isp}) = Cf[isp] ;
                    }
                    //HRR(i,j,k) = WTi;
                    Set_Tfout = true;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double SpeciesTransport::CalculatePhiVanLeer(double r)
{
    return (r+fabs(r))/(1.0+fabs(r));
}

void  SpeciesTransport::CalculateDiffusion(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, double dt)
{
    size_t MixtureComp=0;
    std::vector<int> dir(3);
    double dx = Grid.dx;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MassFractions,0,)
    {
        if(!FL.Obstacle(i,j,k))
        {
            double Tc     = ET.TxOld(i,j,k);
            double Rhoc   = FL.DensityWetting(i,j,k,{MixtureComp});
            double cpc    = ET.Cp_Mixture(i,j,k);
            double ktc    = ET.K_Mixture(i,j,k);
            for (int direction=0; direction<3; ++direction)
            {
                dir[0] = (direction == 0) ? 1*Grid.dNx : 0;
                dir[1] = (direction == 1) ? 1*Grid.dNy : 0;
                dir[2] = (direction == 2) ? 1*Grid.dNz : 0;
                int ip=i+dir[0];
                int jp=j+dir[1];
                int kp=k+dir[2];
                int im=i-dir[0];
                int jm=j-dir[1];
                int km=k-dir[2];
                double Tp   = ET.TxOld(ip,jp,kp);
                double Tm   = ET.TxOld(im,jm,km);
                double ktp  = ET.K_Mixture(ip,jp,kp);
                double ktm  = ET.K_Mixture(im,jm,km);

                ET.Tx(i, j, k) += dt/(Rhoc*cpc) * (ktc * (Tp+Tm-2.0*Tc)/dx/dx + (ktp-ktm)/dx/2.0 * (Tp-Tm)/dx/2.0);
                double Rhop    = FL.DensityWetting(ip,jp,kp,{MixtureComp});
                double Rhom    = FL.DensityWetting(im,jm,km,{MixtureComp});
                int ipp=i+2.0*dir[0];
                int jpp=j+2.0*dir[1];
                int kpp=k+2.0*dir[2];
                int imm=i-2.0*dir[0];
                int jmm=j-2.0*dir[1];
                int kmm=k-2.0*dir[2];
                /*
                double YkVkc1[nSpecies];
                double YkVkp1[nSpecies];
                double YkVkm1[nSpecies];

                double MMWc  = CalculateMeanMolarMass(i, j, k, "old");
                double MMWp  = CalculateMeanMolarMass(ip ,jp ,kp  ,"old");
                double MMWpp = CalculateMeanMolarMass(ipp,jpp,kpp ,"old");
                double MMWm  = CalculateMeanMolarMass(im ,jm ,km  ,"old");
                double MMWmm = CalculateMeanMolarMass(imm,jmm,kmm ,"old");

                double YkVcc=0.0;
                double YkVcp=0.0;
                double YkVcm=0.0;

                for(size_t ic =0; ic < nSpecies; ic++)
                {
                    double Mk    = MolecularWeight({ic});
                    double Xc    = MoleFraction(i,j,k, MMWc ,ic);
                    double Xp    = MoleFraction(ip ,jp ,kp, MMWp,ic);
                    double Xpp   = MoleFraction(ipp,jpp,kpp, MMWpp,ic);
                    double Xm    = MoleFraction(im ,jm ,km, MMWm,ic);
                    double Xmm   = MoleFraction(imm,jmm,kmm, MMWmm,ic);

                    YkVkc1[ic]   = - MassDiff_Species(i,j,k,{ic})   *Mk/MMWc*(Xp-Xm) /dx/2.0;
                    YkVkp1[ic]   = - MassDiff_Species(ip,jp,kp,{ic})*Mk/MMWp*(Xpp-Xc)/dx/2.0;
                    YkVkm1[ic]   = - MassDiff_Species(im,jm,km,{ic})*Mk/MMWm*(Xc-Xmm)/dx/2.0;

                    YkVcc       +=  MassFractionsOld(i,j,k)   ({ic})*MassDiff_Species(i,j,k,{ic})   * Mk/MMWc*(Xp-Xm) /dx/2.0;
                    YkVcp       +=  MassFractionsOld(ip,jp,kp,{ic})*MassDiff_Species(ip,jp,kp,{ic})* Mk/MMWp*(Xpp-Xc)/dx/2.0;
                    YkVcm       +=  MassFractionsOld(im,jm,km,{ic})*MassDiff_Species(im,jm,km,{ic})* Mk/MMWm*(Xc-Xmm)/dx/2.0;
                }

                for(size_t ic =0; ic < nSpecies; ic++)
                {
                    ET.Tx(i, j, k) -= dt/cpc * (Cp_Species(i,j,k,{ic}) * (YkVkc1[ic]+ YkVcc)) * (Tp-Tm)/dx/2.0;
                    MassFractionsOld(i, j, k,{ic}) -= dt/Rhoc * (Rhop * (YkVkp1[ic]+ YkVcp) - Rhom*(YkVkm1[ic]+ YkVcm))/dx/2.0;
                }
                */
                double YkVkc[nSpecies];
                double YkVkp[nSpecies];
                double YkVkm[nSpecies];

                for(size_t ic =0; ic < nSpecies; ic++)
                {
                    double MFc  = MassFractionsOld(i,j,k,{ic});
                    double MFp  = MassFractionsOld(ip,jp,kp,{ic});
                    double MFpp = MassFractionsOld(ipp,jpp,kpp,{ic});
                    double MFm  = MassFractionsOld(im,jm,km,{ic});
                    double MFmm = MassFractionsOld(imm,jmm,kmm,{ic});

                    YkVkc[ic]   = - MassDiff_Species(i,j,k,{ic})   *(MFp-MFm)/dx/2.0;
                    YkVkp[ic]   = - MassDiff_Species(ip,jp,kp,{ic})*(MFpp-MFc)/dx/2.0;
                    YkVkm[ic]   = - MassDiff_Species(im,jm,km,{ic})*(MFc-MFmm)/dx/2.0;
                }

                for(size_t ic =0; ic < nSpecies; ic++)
                {
                    ET.Tx(i, j, k) -= dt/cpc * (Cp_Species(i,j,k,{ic}) * (YkVkc[ic])) * (Tp-Tm)/dx/2.0;
                    MassFractions(i, j, k,{ic}) -= dt/Rhoc * (Rhop * (YkVkp[ic]) - Rhom*(YkVkm[ic]))/dx/2.0;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SpeciesTransport::CalculateReaction(EnergyTransport& ET, FlowSolverLBM& FL, int nDim,  double dt)
{
	size_t MixtureComp=0;
	FuelConsumptionRate=0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MassFractions,0,)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			ET.Tx(i,j,k) += dt/(FL.DensityWetting(i,j,k,{MixtureComp}) * ET.Cp_Mixture(i,j,k)) * HRR(i,j,k);
			for (size_t ic=0; ic<nSpecies; ic++)
    		{
        		MassFractions(i, j, k,{ic}) += dt / FL.DensityWetting(i, j, k,{MixtureComp}) * (W_Species(i, j, k,{ic}));
    		}
    		// Compute Fuel Consumption only if species is fuel
    		FuelConsumptionRate += abs(W_Species(i, j, k,{FuelIndex}) * pow(Grid.dx, nDim));
		}
	}
    OMP_PARALLEL_STORAGE_LOOP_END
	#ifdef MPI_PARALLEL
    	OP_MPI_Allreduce(OP_MPI_IN_PLACE, &FuelConsumptionRate, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
	#endif
}

void SpeciesTransport::CalculateAdvection(EnergyTransport& ET, FlowSolverLBM& FL, FlowMixture& FM, double dt)
{
    std::vector<int> dir(3);
    double eps=1e-10;
    double dx = Grid.dx;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MassFractions,0,)
    {	
        if(!FL.Obstacle(i,j,k))
        {
            double Tc = ET.TxOld(i,j,k);
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

                double Tp = ET.TxOld(ip,jp,kp);
                double Tm = ET.TxOld(im,jm,km);
                double Vc = FM.V_Mixture(i,j,k)[direction];

                if(AdvectionUpwind)
                {
                    if(Vc>=0)
                    {
                        ET.Tx(i,j,k) += -dt/dx * Vc * (Tc-Tm);
                        for(size_t n =0; n < nSpecies; n++)
                        {
                            double C   = MassFractionsOld(i,j,k,{n});
                            double Cm  = MassFractionsOld(im,jm,km,{n});
                            MassFractions(i, j, k,{n}) += ( -dt/dx * Vc * (C-Cm));
                        }
                    }
                    else
                    {
                        ET.Tx(i,j,k) += -dt/dx * Vc * (Tp-Tc);
                        for(size_t n =0; n < nSpecies; n++)
                        {
                            double C   = MassFractionsOld(i,j,k,{n});
                            double Cp  = MassFractionsOld(ip,jp,kp,{n});
                            MassFractions(i, j, k,{n}) += ( -dt/dx * Vc * (Cp-C));
                        }
                    }                       
                }
                else if(AdvectionVanLeer)
			    {
                    double Thp = (Tc+Tp)/2.0;
			    	double Thm = (Tc+Tm)/2.0;
			    	if(Vc>=0)
                	{
                        int imm=i-2.0*dir[0];
    		            int jmm=j-2.0*dir[1];
    		            int kmm=k-2.0*dir[2];
                        double Tmm  = ET.TxOld(imm,jmm,kmm);
                        double r    = (Tc-Tm)/(Tp-Tc+eps);
			            double phi  = (r+fabs(r))/(1.0+fabs(r));
			            Thp         = Tc+0.5*phi*(Tp-Tc);
			            double rm  = (Tm-Tmm)/(Tc-Tm+eps);
			            double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			            Thm=Tm+0.5*phim*(Tc-Tm);
			            ET.Tx(i,j,k) +=  -dt*Vc*(Thp-Thm)/dx;
                        for(size_t n =0; n < nSpecies; n++)
                        {
                            double C    = MassFractionsOld(i,j,k,{n});
                            double Cp    = MassFractionsOld(ip,jp,kp,{n});
                            double Cm   = MassFractionsOld(im,jm,km,{n});
                            double Cmm  = MassFractionsOld(imm,jmm,kmm,{n});
                            double r    = (C-Cm)/(Cp-C+eps);
			                double phi  = (r+fabs(r))/(1.0+fabs(r));
			                double Chp         = C+0.5*phi*(Cp-C);
			                double rm   = (Cm-Cmm)/(C-Cm+eps);
			                double phim = (rm+fabs(rm))/(1.0+fabs(rm));
			                double Chm         = Cm+0.5*phim*(C-Cm);
			                MassFractions(i, j, k,{n}) +=  -dt * Vc * (Chp-Chm)/dx;
                        }
    	        	}
    	        	else
                	{
                        int ipp=i+2.0*dir[0];
    		            int jpp=j+2.0*dir[1];
    		            int kpp=k+2.0*dir[2];
                        double Tpp  = ET.TxOld(ipp,jpp,kpp);
			            double r  = (Tpp-Tp)/(Tp-Tc+eps);
			            double phi=(r+fabs(r))/(1.0+fabs(r));
			            Thp=Tp-0.5*phi*(Tp-Tc);
			            double rm  = (Tp-Tc)/(Tc-Tm+eps);
			            double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			            Thm=Tc-0.5*phim*(Tc-Tm);
			            ET.Tx(i,j,k)  +=  -dt*Vc*(Thp-Thm)/dx;

                        for(size_t n =0; n < nSpecies; n++)
                        {
                            double C    = MassFractionsOld(i,j,k,{n});
                            double Cm   = MassFractionsOld(im,jm,km,{n});
                            double Cp    = MassFractionsOld(ip,jp,kp,{n});
                            double Cpp  = MassFractionsOld(ipp,jpp,kpp,{n});
                            double r    = (Cpp-Cp)/(Cp-C+eps);
			                double phi  = (r+fabs(r))/(1.0+fabs(r));
			                double Chp  = Cp-0.5*phi*(Cp-C);
			                double rm   = (Cp-C)/(C-Cm+eps);
			                double phim = (rm+fabs(rm))/(1.0+fabs(rm));
			                double Chm  = C-0.5*phim*(C-Cm);
			                MassFractions(i, j, k,{n}) +=  -dt * Vc * (Chp-Chm)/dx;
                        }
                	}
			    }
                else if (AdvectionCentral)
                {
                    ET.Tx(i,j,k) += -dt/dx * Vc * (Tp-Tm)/2.0;
                    for(size_t n =0; n < nSpecies; n++)
                    {
                        double Cm    = MassFractionsOld(im,jm,km,{n});
                        double Cp    = MassFractionsOld(ip,jp,kp,{n});
                        MassFractions(i, j, k,{n}) += ( -dt/dx * Vc * (Cp-Cm)/2.0);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
} 

void SpeciesTransport::WriteVTKMassFractions(Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < nSpecies; comp++)
    {
        const std::string nameComp1 = "MassFraction_" + locSettings.ElementNames[comp];
        ListOfFields.push_back((VTK::Field_t) {nameComp1, [this, comp](int i,int j,int k){return MassFractions(i,j,k,{comp});}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"MassFractions_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void SpeciesTransport::WriteVTKHHR(Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"omega_T (W/m^3)", [this](int i,int j,int k){return HRR(i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"HHR_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void SpeciesTransport::SetInitial(PhaseField& Phase, FlowSolverLBM& FL, EnergyTransport& ET,
                            const vector<double>& BM, const vector<double>& FM , const vector<double>& AM, 
                            double TB, double TF, double TA)
{
    //Simple composition set up, uniform component fractions per phase.
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, MassFractions, MassFractions.Bcells(),)
    {
        HRR(i,j,k) = 0.0;
        if( (i+Phase.Grid.OffsetX)>=X0BurntZone and (i+Phase.Grid.OffsetX)<XNBurntZone and
            (j+Phase.Grid.OffsetY)>=Y0BurntZone and (j+Phase.Grid.OffsetY)<YNBurntZone and
            (k+Phase.Grid.OffsetZ)>=Z0BurntZone and (k+Phase.Grid.OffsetZ)<ZNBurntZone and !FL.Obstacle(i,j,k))
        {
            ET.Tx(i,j,k) = TB;
            for (size_t ic=0; ic<nSpecies; ic++)
            {
                MassFractions(i,j,k,{ic}) = BM[ic];
            }
        }
        else if( (i+Phase.Grid.OffsetX)>=X0FuelZone and (i+Phase.Grid.OffsetX)<XNFuelZone and
                 (j+Phase.Grid.OffsetY)>=Y0FuelZone and (j+Phase.Grid.OffsetY)<YNFuelZone and
                 (k+Phase.Grid.OffsetZ)>=Z0FuelZone and (k+Phase.Grid.OffsetZ)<ZNFuelZone and !FL.Obstacle(i,j,k))
        {
            ET.Tx(i,j,k) = TF;
            for (size_t ic=0; ic<nSpecies; ic++)
            {
                MassFractions(i,j,k,{ic}) = FM[ic];
            }
        }
        else 
        {
            if(!FL.Obstacle(i,j,k))
            {
                ET.Tx(i,j,k) = TA;
                for (size_t ic=0; ic<nSpecies; ic++)
                {
                    MassFractions(i,j,k,{ic}) = AM[ic];
                }
            }
            else
            {
                for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
                {
                    if(Phase.FieldsProperties[it->index].State ==AggregateStates::Solid)
                    {
                        int PFIdx=it->index;
                        ET.Tx(i,j,k) = ET.SurfaceTemp[PFIdx];
                    }
                }
                for (size_t ic=0; ic<nSpecies; ic++)
                {
                    MassFractions(i,j,k,{ic}) = AM[ic];
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(ET.Conjugate)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, ET.Tx, ET.Tx.Bcells(),)
        {
            ET.Ts(i,j,k)=ET.Tx(i,j,k);
            ET.TsOld(i,j,k)=ET.Tx(i,j,k);
            ET.Cp_Solid(i,j,k)=ET.Cps;
            ET.K_Solid(i,j,k)=ET.Ks;
            ET.Density_Solid(i,j,k)=ET.Rhos;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}
