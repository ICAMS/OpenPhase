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
#include "ReactiveFlows/EnergyTransport.h"
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

void EnergyTransport::ReadInput(string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("EnergyTransport");
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
    int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "EnergyTransport");
    Pr					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Pr"), false, 0.71);
    Cp					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Cp"), false, 1005.0);
    Mw					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Mw"), false, 0.02896);
    AdvectionUpwind    	 = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("AdvectionUpwind"), false, false);
    AdvectionCentral     = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("AdvectionCentral"), false, false);
    AdvectionVanLeer     = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("AdvectionVanLeer"), false, true);
    Conjugate		     = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("Conjugate"), false, false);
    Ts0					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Ts0"), false, 300.0);
    Cps					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Cps"), false, 900.0);
    Rhos				 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Rhos"), false, 5000.0);
    Ks					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Ks"), false, 237);
    T0					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("T0"), false, 300.0);
    ATstar  			 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ATstar"), false, 100.0);
    TempSolid	    	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TempSolid"), false, 300.0);
    TempSolidCold     	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TempSolidCold"), false, 300.0);
    TMu0              	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TMu0"), false, 273.0);
    SMu0     	         = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("SMu0"), false, 110.5);
    Mu0     	         = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Mu0"), false, 1.68e-05);
    IF_ENERGY            = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("IF_ENERGY"), false, true);
    BCOrder              = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("BCOrder"), false, 0);
}

void EnergyTransport::Initialize(Settings& locSettings)
{
    Grid = locSettings.Grid;
    size_t Bcells = locSettings.Grid.Bcells;
    Tx.Allocate         (locSettings.Grid, Bcells);
    TxOld.Allocate      (locSettings.Grid, Bcells);
    Cp_Mixture.Allocate (locSettings.Grid, Bcells);
    K_Mixture.Allocate  (locSettings.Grid, Bcells);

    if(Conjugate)
    {
        Ts.Allocate            (locSettings.Grid, Bcells);
        TsOld.Allocate         (locSettings.Grid, Bcells);
        Cp_Solid.Allocate      (locSettings.Grid, Bcells);
        K_Solid.Allocate       (locSettings.Grid, Bcells);
        Density_Solid.Allocate (locSettings.Grid, Bcells);
    }
}

void EnergyTransport::UpdateFields()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN (i,j,k,Tx,Tx.Bcells(),)
    {
        TxOld(i,j,k)=Tx(i,j,k);
        if(Conjugate)
        {
            TsOld(i,j,k)=Ts(i,j,k);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EnergyTransport::CalculateMixtureSpecificHeatCapacityAndThermalConductivity(FlowSolverLBM& FL)
{
    size_t MixtureComp=0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cp_Mixture, Cp_Mixture.Bcells(),)
    {
        Cp_Mixture(i,j,k) = Cp;
        K_Mixture(i,j,k)=FL.DensityWetting(i,j,k,{MixtureComp})*FL.nut(i,j,k,{MixtureComp})*Cp_Mixture(i,j,k)/Pr;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EnergyTransport::SetFreeBCNX(PhaseField& Phase)
{
    int iEnd = Grid.Nx;
    double Temp=T0; 		// Temporary variable for extrapolated temperature
    // Loop over the extrapolated boundary nodes
    if(Grid.OffsetX+Grid.Nx==Grid.TotalNx)
    for (int i = iEnd; i < Grid.Nx + Tx.Bcells(); ++i)
    {
        for (int j = 0; j< Grid.Ny ; ++j)
        {
            for (int k = 0; k< Grid.Nz ; ++k)
            {	
                if (BCOrder==0)
                {
                    Temp = Tx(iEnd - 1, j, k);
                }
                else if (BCOrder==1)
                {
                    Temp = 2.0 * Tx(iEnd - 1, 0, k) - Tx(iEnd - 2, 0, k);
                }
                else if (BCOrder==2)
                {
                    Temp = Tx(iEnd - 1, j, k) + (Tx(iEnd - 1, j, k) - Tx(iEnd - 2, j, k)) +
                           0.5 * (Tx(iEnd - 1, j, k) - 2.0 * Tx(iEnd - 2, j, k) + Tx(iEnd - 3, j, k));
                }
                else if (BCOrder==3)
                {
                    Temp = Lagrange_polynomial4(0.0, 1.0, 2.0, 3.0, 4.0,
                                                Tx(iEnd - 4, j, k), Tx(iEnd - 3, j, k),
                                                Tx(iEnd - 2, j, k), Tx(iEnd - 1, j, k));
                }
                Tx(i, j, k) = Temp;
            }
        }
    }
}

void EnergyTransport::SetBoundaryConditions(BoundaryConditions& BC)
{
    if (Grid.dNx > 0) BC.SetX(Tx);
    if (Grid.dNy > 0) BC.SetY(Tx);
    if (Grid.dNz > 0) BC.SetZ(Tx);

    if(Conjugate)
    {
        if (Grid.dNx > 0) BC.SetX(Ts);
        if (Grid.dNy > 0) BC.SetY(Ts);
        if (Grid.dNz > 0) BC.SetZ(Ts);
    }
}

void EnergyTransport::UpdateGhostPoints(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB)
{
    if(Conjugate) CalculateGhostPointsConjugateHeatTransfer(Phase, FL, SB);
    if(!Conjugate) CalculateGhostPoints(Phase, FL, SB);
}

void EnergyTransport::CalculateGhostPointsConjugateHeatTransfer(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB)
{
    size_t Gasidx   = 0;
    size_t Solididx = 1;
    double ep=1e-8;
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

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
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
                                                          Ts(i0*idi,j0*jdi,k0*kdi),Ts(i1*idi,j0*jdi,k0*kdi),
                                                          Ts(i0*idi,j1*jdi,k0*kdi),Ts(i1*idi,j1*jdi,k0*kdi),
                                                          Ts(i0*idi,j0*jdi,k1*kdi),Ts(i1*idi,j0*jdi,k1*kdi),
                                                          Ts(i0*idi,j1*jdi,k1*kdi),Ts(i1*idi,j1*jdi,k1*kdi));
                    double Ksi  =  SB.TrilinearInterPolation(ic,jc,kc,
                                                          K_Solid(i0*idi,j0*jdi,k0*kdi),K_Solid(i1*idi,j0*jdi,k0*kdi),
                                                          K_Solid(i0*idi,j1*jdi,k0*kdi),K_Solid(i1*idi,j1*jdi,k0*kdi),
                                                          K_Solid(i0*idi,j0*jdi,k1*kdi),K_Solid(i1*idi,j0*jdi,k1*kdi),
                                                          K_Solid(i0*idi,j1*jdi,k1*kdi),K_Solid(i1*idi,j1*jdi,k1*kdi));
                    double Tf  = Tx(i,j,k);
                    double Kf  = K_Mixture(i,j,k);

                    double Tb  = (Kf/X*Tf+Ksi/D*Tsi)/(Kf/X+Ksi/D);
                    Ts (i,j,k) = Tb * (1.0+X/D) - X/D*Tsi;
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
                                                        Tx(i0*idi,j0*jdi,k0*kdi),Tx(i1*idi,j0*jdi,k0*kdi),
                                                        Tx(i0*idi,j1*jdi,k0*kdi),Tx(i1*idi,j1*jdi,k0*kdi),
                                                        Tx(i0*idi,j0*jdi,k1*kdi),Tx(i1*idi,j0*jdi,k1*kdi),
                                                        Tx(i0*idi,j1*jdi,k1*kdi),Tx(i1*idi,j1*jdi,k1*kdi));
                    double Kfi  = SB.TrilinearInterPolation(ic,jc,kc,
                                                        K_Mixture(i0*idi,j0*jdi,k0*kdi),K_Mixture(i1*idi,j0*jdi,k0*kdi),
                                                        K_Mixture(i0*idi,j1*jdi,k0*kdi),K_Mixture(i1*idi,j1*jdi,k0*kdi),
                                                        K_Mixture(i0*idi,j0*jdi,k1*kdi),K_Mixture(i1*idi,j0*jdi,k1*kdi),
                                                        K_Mixture(i0*idi,j1*jdi,k1*kdi),K_Mixture(i1*idi,j1*jdi,k1*kdi));
                    double Tsi  = Ts(i,j,k);
                    double Ksi  = K_Solid(i,j,k);
                    if(X<ep)
                    {
                        Tx (i,j,k) = Tsi;
                    }
                    else
                    {
                        double Tb  = (Kfi/D*Tfi+Ksi/X*Tsi)/(Kfi/D+Ksi/X);
                        Tx(i,j,k) =  Tb * (1.0+X/D) - X/D * Tfi;
                    }
                    Set_TfOut=true;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EnergyTransport::CalculateGhostPoints(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB)
{
    size_t Gasidx = 0;
    size_t Solididx = 1;
    double Tb=TempSolid;
    double Qb=HeatFlux;
    double ep=1e-8;
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

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    {
        bool Set_Tfout=false;
        if(FL.Obstacle(i,j,k))
        for(int ii = -iBC; ii <= iBC; ii++)
        for(int jj = -jBC; jj <= jBC; jj++)
        for(int kk = -kBC; kk <= kBC; kk++)
        {
            if(!Set_Tfout)
            {
                if(!FL.Obstacle(i+ii,j+jj,k+kk))
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
                    dVector3 Xf = Xb + Norm*SB.nDist;

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
                                                          Tx(i0*idi,j0*jdi,k0*kdi),Tx(i1*idi,j0*jdi,k0*kdi),
                                                          Tx(i0*idi,j1*jdi,k0*kdi),Tx(i1*idi,j1*jdi,k0*kdi),
                                                          Tx(i0*idi,j0*jdi,k1*kdi),Tx(i1*idi,j0*jdi,k1*kdi),
                                                          Tx(i0*idi,j1*jdi,k1*kdi),Tx(i1*idi,j1*jdi,k1*kdi));
                    double Kf   =  SB.TrilinearInterPolation(ic,jc,kc,
                                                          K_Mixture(i0*idi,j0*jdi,k0*kdi),K_Mixture(i1*idi,j0*jdi,k0*kdi),
                                                          K_Mixture(i0*idi,j1*jdi,k0*kdi),K_Mixture(i1*idi,j1*jdi,k0*kdi),
                                                          K_Mixture(i0*idi,j0*jdi,k1*kdi),K_Mixture(i1*idi,j0*jdi,k1*kdi),
                                                          K_Mixture(i0*idi,j1*jdi,k1*kdi),K_Mixture(i1*idi,j1*jdi,k1*kdi));
                    if(IF_ConstTemp[Solididx])
                    {
                        Tb = SurfaceTemp[Solididx];
                        if(X<ep)
                        {
                            Tx(i,j,k) = Tb;
                        }
                        else
                        {
                            Tx(i,j,k) =  Tb * (1.0+X/D) - X/D * Tf;
                        }
                    }
                    else if(IF_ConstFlux[Solididx])
                    {
                        Qb = SurfaceFlux[Solididx];
                        Tx(i,j,k) = Tf-(X+D)*Qb/Kf;
                    }                    
                    Set_Tfout = true;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double EnergyTransport::CalculatePhiVanLeer(double r)
{
    return (r+fabs(r))/(1.0+fabs(r));
}

void EnergyTransport::CalculateThermodynamicPressure(FlowSolverLBM& FL, double P0, double T0)
{
    double SumTinv=0.0;
    double fcells=0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    {
        if(!FL.Obstacle(i,j,k))
        {
            SumTinv += 1.0/Tx(i,j,k);
            fcells++; 
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    #ifdef MPI_PARALLEL
    	OP_MPI_Allreduce(OP_MPI_IN_PLACE, &SumTinv, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    	OP_MPI_Allreduce(OP_MPI_IN_PLACE, &fcells, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
	#endif
    FL.Pth = P0/T0*fcells/SumTinv;    
}

void EnergyTransport::CalculateDiffusion(PhaseField& Phase, FlowSolverLBM& FL, double dt)
{
    size_t MixtureComp=0;
    std::vector<int> dir(3);
    double dx	  = Grid.dx;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    {
        if(!FL.Obstacle(i,j,k))
        {
            double Tc     = TxOld(i,j,k);
            double Rhoc   = FL.DensityWetting(i,j,k,{MixtureComp});
            double cpc    = Cp_Mixture(i,j,k);
            double ktc    = K_Mixture(i,j,k);
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

                double Tp     = TxOld(ip,jp,kp);
                double Tm     = TxOld(im,jm,km);
                double ktp    = K_Mixture(ip,jp,kp);
                double ktm    = K_Mixture(im,jm,km);

                Tx(i, j, k) += dt/(Rhoc*cpc) * (ktc * (Tp+Tm-2.0*Tc)/dx/dx + (ktp-ktm)/dx/2.0 * (Tp-Tm)/dx/2.0);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EnergyTransport::CalculatePhaseFieldCoupling(PhaseField& Phase, FlowSolverLBM& FL, double hT, double dt)
{
    double Ts=TempSolid;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    if(!FL.Obstacle(i,j,k))
    {
        double SF=0;
        for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            size_t PhaseIdx =  Phase.FieldsProperties[it->index].Phase;
            if(Phase.FieldsProperties[PhaseIdx].State == AggregateStates::Solid)
            {
                int PFIdx=it->index;
                Ts = SurfaceTemp[PFIdx];
                SF += it->value;
            }
        }
        double ST  = hT * K_Mixture(i,j,k)* pow(SF,2)*(1.0-SF) / pow(Phase.Grid.iWidth*Grid.dx,2.0) * (TxOld(i,j,k)-Ts);
        Tx(i,j,k) -=   dt / FL.DensityWetting(i,j,k,{0}) / Cp_Mixture(i,j,k) * ST; 
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EnergyTransport::CalculateAdvection(FlowSolverLBM& FL, FlowMixture& FM, double dt)
{
    std::vector<int> dir(3);
    double eps=1e-10;
    double dx = Grid.dx;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    {	
        if(!FL.Obstacle(i,j,k))
        {
            double Tc = TxOld(i,j,k);
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

                double Tp = TxOld(ip,jp,kp);
                double Tm = TxOld(im,jm,km);
                double Vc = FM.V_Mixture(i,j,k)[direction];

                if(AdvectionUpwind)
                {
                    if(Vc>=0)
                    {
                        Tx(i,j,k) += -dt/dx * Vc * (Tc-Tm);
                    }
                    else
                    {
                        Tx(i,j,k) += -dt/dx * Vc * (Tp-Tc);
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
                        double Tmm  = TxOld(imm,jmm,kmm);
                        double r    = (Tc-Tm)/(Tp-Tc+eps);
			            double phi  = (r+fabs(r))/(1.0+fabs(r));
			            Thp         = Tc+0.5*phi*(Tp-Tc);
			            double rm  = (Tm-Tmm)/(Tc-Tm+eps);
			            double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			            Thm=Tm+0.5*phim*(Tc-Tm);
			            Tx(i,j,k) +=  -dt*Vc*(Thp-Thm)/dx;
    	        	}
    	        	else
                	{
                        int ipp=i+2.0*dir[0];
    		            int jpp=j+2.0*dir[1];
    		            int kpp=k+2.0*dir[2];
                        double Tpp  = TxOld(ipp,jpp,kpp);
			            double r  = (Tpp-Tp)/(Tp-Tc+eps);
			            double phi=(r+fabs(r))/(1.0+fabs(r));
			            Thp=Tp-0.5*phi*(Tp-Tc);
			            double rm  = (Tp-Tc)/(Tc-Tm+eps);
			            double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			            Thm=Tc-0.5*phim*(Tc-Tm);
			            Tx(i,j,k)  +=  -dt*Vc*(Thp-Thm)/dx;
                	}
			    }
                else if (AdvectionCentral)
                {
                    Tx(i,j,k) += -dt/dx * Vc * (Tp-Tm)/2.0;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
} 

void EnergyTransport::CalculateSolidDiffusion(PhaseField& Phase, FlowSolverLBM& FL, double dt)
{
    std::vector<int> dir(3);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Ts,0,)
    {
        if(FL.Obstacle(i,j,k))
        {
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
                double dx	  = Grid.dx;
                double Tc     = TsOld(i,j,k);
                double Tp     = TsOld(ip,jp,kp);
                double Tm     = TsOld(im,jm,km);
                double Rhoc   = Density_Solid(i,j,k);
                double cpc    = Cp_Solid(i,j,k);
                double ktc    = K_Solid(i,j,k);
                Ts(i, j, k) += dt/(Rhoc*cpc) * (ktc * (Tp+Tm-2.0*Tc)/dx/dx);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EnergyTransport::WriteVTKTemperature(Settings& locSettings, FlowSolverLBM& FL, const int tStep, const int precision)
{
    if(Conjugate) IntegrateTemperature(FL);
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Temperature", [this](int i,int j,int k){return Tx(i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"Temperature_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void EnergyTransport::IntegrateTemperature(FlowSolverLBM& FL)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Ts,0,)
    {
        if (FL.Obstacle(i,j,k))
        {
            Tx(i,j,k)=Ts(i,j,k);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double EnergyTransport::Lagrange_polynomial4(double T1, double T2, double T3, double T4, double T5,
                                          double f1, double f2, double f3, double f4)
{
    return ((T5 - T2) * (T5 - T3) * (T5 - T4) * f1) / (-6.0) +
           ((T5 - T1) * (T5 - T3) * (T5 - T4) * f2) / (2.0) +
           ((T5 - T1) * (T5 - T2) * (T5 - T4) * f3) / (-2.0) +
           ((T5 - T1) * (T5 - T2) * (T5 - T3) * f4) / (6.0);
}

void EnergyTransport::SetSolidPhaseTemp(PhaseField& Phase, bool DI)
{
    double CRI = (DI)?(1.0-DBL_EPSILON):(0.5-DBL_EPSILON);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    {
        for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            size_t PhaseIdx =  Phase.FieldsProperties[it->index].Phase;
            if(Phase.FieldsProperties[PhaseIdx].State == AggregateStates::Solid)
            {
                if(Phase.Fractions(i,j,k)[PhaseIdx]>=CRI)
                {
                    int PFIdx = it->index;
                    Tx(i, j, k) = SurfaceTemp[PFIdx];
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EnergyTransport::SetInitial(PhaseField& Phase, FlowSolverLBM& FL)
{
    
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, Tx, Tx.Bcells(),)
    {
        Tx(i,j,k) = T0;
        if(Conjugate)
        {
            Ts(i,j,k) = T0;
            Cp_Solid(i,j,k)=Cps;
            K_Solid(i,j,k)=Ks;
            Density_Solid(i,j,k)=Rhos;
        }
        if(FL.Obstacle(i,j,k))
        {
            Ts(i,j,k) = Ts0;
            Tx(i,j,k) = Ts0;
        }   
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}