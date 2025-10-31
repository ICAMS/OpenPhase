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

#include "ReactiveFlows/FlowMixture.h"
#include "Containers/dVector3.h"

using namespace std;
using namespace openphase;

constexpr double lbcs2 = 1.0/3.0;                                               ///< Speed of sound squared [lattice units]

void FlowMixture::ReadInput(string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("FlowMixture");
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
    int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "FlowMixture");
    LengthScale    	     = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("LengthScale"), false, 1.0);
    Updating_Velocity    = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("Update_Vel"), false, false);
    UniformVel           = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("UniformVel"), false, true);
    PoiseuilleVel        = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("PoiseuilleVel"), false, false);
    Austar               = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Austar"), false, 300.0);
    Species              = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("Species"), false, false);
    DI                   = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("DI"), false, false);
    DO_ZCORNER           = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("ZCORNER"), false, false);
    DO_YCORNER           = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("YCORNER"), false, false);
    SecOrdBB             = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("SecOrdBB"), false, false);
    lbnu                 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("lbnu"), false, 0.005);
    MaxU                 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("MaxU"), false, 1.0);
    Kp                   = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Kp"), false, 1.0);
}
void FlowMixture::Initialize(Settings& locSettings, PhaseField& Phase, FlowSolverLBM& FL, Velocities& Vel, BoundaryConditions& BC)
{
    V_Mixture.Allocate(locSettings.Grid, locSettings.Grid.Bcells);

    SetInitialVelocity(Phase, BC, FL, Vel, FL.U0X);
    InitializingFlowProperties(Phase,FL,Vel, BC);

    InletMassFlowRate=InletDensity*InletArea*FL.U0X;

    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
    {
        cout<<"Inlet area (m^2): "<<InletArea<<endl;
        cout<<"Inlet mass flow (kg/s): "<<InletMassFlowRate<<endl;
    }

    //FL.U0X = InletMassFlowRate/InletDensity/InletArea;
}

void FlowMixture::WriteVTKMixtureVelocity(Settings& locSettings,  const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"Velocity_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t){"Velocity",  [this](int i,int j,int k){return V_Mixture(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields);
}

void FlowMixture::UpdateMixtureVelocity(Velocities& Vel)
{
    size_t Flowidx=0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,V_Mixture, V_Mixture.Bcells(),)
    {
        V_Mixture (i,j,k)  = Vel.Phase(i,j,k,{Flowidx});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowMixture::LBMLimits(PhaseField& Phase, FlowSolverLBM& FL, double Umax, double lbnu, double dt)
{

    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
    {
        double dxt  = Umax/0.1;
        double dx = Phase.Grid.dx;

        double nuphy=InletViscosity;
        double maxdx=nuphy/lbnu/dxt;
        double maxdt=maxdx/dxt;

        cout<<"For lbnu = "<<lbnu<<", Umax = "<<Umax<<", maximum of dx is: "<< maxdx<<", and max timestep is: "<<maxdt<<endl;

        cout<<"lbnu for running simulation is: "<< nuphy * dt/dx/dx<<endl;
        cout<<"Physical kinematic viscosity for running simulation is: "<< nuphy<<endl;
        cout<<"Inlet Density: "<<InletDensity<<endl;
        cout<<"Reynolds number (U Lz/nu) for running simulation: "<<FL.U0X*LengthScale/InletViscosity<<endl;
    }
}

void FlowMixture::SetInitialVelocity(PhaseField& Phase, BoundaryConditions& BC, FlowSolverLBM& FL, Velocities& Vel, double Um)
{
    size_t gasindex=0;
    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
        if(Phase.FieldsProperties[idx].State == AggregateStates::Gas)
        {
            gasindex=idx;
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Vel.Average,0,)
    if (!FL.Obstacle(i,j,k))
    {
        Vel.Phase(i,j,k,{gasindex})={Um,0.0,0.0};
    }
    else
    {
        Vel.Phase(i,j,k,{gasindex})={0.0,0.0,0.0};
    }

    OMP_PARALLEL_STORAGE_LOOP_END
}
void FlowMixture::InitializingFlowProperties(PhaseField& Phase, FlowSolverLBM& FL, Velocities& Vel, BoundaryConditions& BC)
{
    bool getfirst=false;
    InletDensity=0.0;
    InletViscosity=0.0;
    int iref=0;
    cout<<"Calculate inlet density at x node: "<<iref<<endl;
    if(Phase.Grid.OffsetX==0.0 and Phase.Grid.OffsetZ==0.0)
    {
        for (int j = 0; j < FL.Grid.Ny; j++)
        {
            for (int k = 0; k < FL.Grid.Nz; k++)
            {
                if(!FL.Obstacle(iref,j,k))
                {
                    if(!getfirst)
                    {
                        InletDensity=FL.DensityWetting(iref,j,k,{0});
                        InletViscosity=FL.nut(iref,j,k,{0});
                        getfirst=true;
                    }
                }
            }
        }
        cout<<"Inlet density: "<<InletDensity<<endl;
        cout<<"Inlet viscosity: "<<InletViscosity<<endl;
    }
    #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &InletDensity, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &InletViscosity, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    if(InletDensity==0.0)
    {
        iref = FL.Grid.Nx/2.0;
        cout<<"Calculate inlet density at x node: "<<iref<<endl;
        if(Phase.Grid.OffsetX==0.0 and Phase.Grid.OffsetZ==0.0)
        {
            for (int j = 0; j < FL.Grid.Ny; j++)
            {
                for (int k = 0; k < FL.Grid.Nz; k++)
                {
                    if(!FL.Obstacle(iref,j,k))
                    {
                        if(!getfirst)
                        {
                            InletDensity=FL.DensityWetting(iref,j,k,{0});
                            InletViscosity=FL.nut(iref,j,k,{0});
                            getfirst=true;
                        }
                    }
                }
            }
            cout<<"Inlet density: "<<InletDensity<<endl;
            cout<<"Inlet viscosity: "<<InletViscosity<<endl;
        }
        #ifdef MPI_PARALLEL
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &InletDensity, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &InletViscosity, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif
    }

    double nodes=0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,0,)
    if (!FL.Obstacle(i,j,k) and i == 0 and j ==0)
    {
        nodes++;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    InletArea = nodes * FL.Grid.dx;
    FL.SetObstacleNodes(Phase, Vel);
    FL.SetInitialPopulationsTC(BC,Vel);
    FL.SetBoundaryConditions(BC);
}

void FlowMixture::Collision(FlowSolverLBM& FL, Velocities& Vel, BoundaryConditions& BC)
{
    FL.CollisionTC(Vel);
    if(FL.Grid.dNx) BC.SetXVector(FL.lbPopulations);
    if(FL.Grid.dNy) BC.SetYVector(FL.lbPopulations);
    if(FL.Grid.dNz) BC.SetZVector(FL.lbPopulations);
}

void FlowMixture::Propagation(FlowSolverLBM& FL, PhaseField& Phase, BoundaryConditions& BC)
{
    if (SecOrdBB) 
    {
        FL.PropagationSecondOrderBB(Phase,BC);
    } 
    else 
    {
        FL.Propagation(Phase, BC);
    }
    if(FL.Grid.dNx) BC.SetXVector(FL.lbPopulations);
    if(FL.Grid.dNy) BC.SetYVector(FL.lbPopulations);
    if(FL.Grid.dNz) BC.SetZVector(FL.lbPopulations);
}

void FlowMixture::CalculateForceDragbyPF(PhaseField& Phase, FlowSolverLBM& FL, const Velocities& Vel, double hHyd )
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,0,)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        if((i+Phase.Grid.OffsetX) > 0)
        if(!FL.Obstacle(i,j,k))
        {
            const double mu = FL.DensityWetting(i,j,k,{n})*FL.nut(i,j,k,{n});  //Dynamic viscosity [kg/m/s]
            dVector3 vel = Vel.Phase(i,j,k,{n});
            double SF = FL.SolidFraction(i,j,k,Phase);
            double CF = hHyd * mu /pow(Phase.Grid.iWidth*FL.Grid.dx,2.0) * pow(SF,2)*(1.0-SF);
            FL.ForceDensity(i,j,k,{n}) +=  vel *(-CF); 
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double FlowMixture::CalculateSutherlandViscosity(double Temp)
{
    double Ts=273.0;
    double mus=1.68e-05; //it's for air
    double S=110.5;
    double mu=mus*(Temp/Ts)*sqrt(Temp/Ts)*(Ts+S)/(Temp+S);
    return mu;
}
double FlowMixture::CalculateIdealGasDensity(double p, double Rm, double Temp)
{
    double rho= p/(Rm * Temp);
    return rho;
}

void FlowMixture::UpdateFluidProperties(FlowSolverLBM& FL, EnergyTransport& ET, double mw) 
{
    double R = 8.314462618; // J/(mol*K)
    size_t MixtureComp = 0 ;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.nut, FL.nut.Bcells(),)
    {
        FL.DensityWetting(i,j,k,{MixtureComp})=CalculateIdealGasDensity(FL.Pth, R/mw, ET.Tx(i,j,k));
        FL.nut(i,j,k,{MixtureComp})=CalculateSutherlandViscosity(ET.Tx(i,j,k))/FL.DensityWetting(i,j,k,{MixtureComp});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowMixture::UpdatingInletVelocity(FlowSolverLBM& FL, double FuelConsumptionRate, double FuelMF)
{
    FL.U0X  = (FuelConsumptionRate) / (FuelMF * InletDensity * InletArea);
}

void FlowMixture::DetectObstacles(FlowSolverLBM& FL, const PhaseField& Phase, bool DI)
{
    bool locObstaclesChanged = false;
	double CRI = (DI)?(1.0-DBL_EPSILON):(0.5-DBL_EPSILON);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,FL.DensityWetting.Bcells(), reduction(or:locObstaclesChanged))
    {
        if (FL.SolidFraction(i,j,k,Phase) >= CRI)
        {
            if (not FL.Obstacle(i,j,k)) locObstaclesChanged = true;
            FL.Obstacle(i,j,k) = true;
        }
        else
        {
            if (FL.Obstacle(i,j,k)) locObstaclesChanged = true;
            FL.Obstacle(i,j,k) = false;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    FL.ObstaclesChanged = locObstaclesChanged;
}

void FlowMixture::CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL,  EnergyTransport& ET, BoundaryConditions& BC, double dt)
{
    double CRI = (DI)?(1.0-DBL_EPSILON):(0.5-DBL_EPSILON);
    std::vector<int> dir(3);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DivVel,0,)
    {
        for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            size_t phaseIdx =  Phase.FieldsProperties[it->index].Phase;
            if(Phase.PhaseAggregateStates[phaseIdx] == AggregateStates::Gas)
            {
                if(Phase.Fractions(i,j,k)[phaseIdx]>=CRI)
                {
                    double Tc     = ET.Tx(i,j,k);
                    double TempoT = (ET.Tx(i,j,k) - ET.TxOld(i,j,k))/dt/Tc;
                    double SpatiT = 0.0;
                    
                    for (int direction=0; direction<3; ++direction)
                    {
                        dir[0] = (direction == 0) ? 1*ET.Grid.dNx : 0;
                        dir[1] = (direction == 1) ? 1*ET.Grid.dNy : 0;
                        dir[2] = (direction == 2) ? 1*ET.Grid.dNz : 0;
                        double Vc 	= V_Mixture(i,j,k)[direction];
                        double Tp   = ET.Tx(i+dir[0], j+dir[1], k+dir[2]);
                        double Tm   = ET.Tx(i-dir[0], j-dir[1], k-dir[2]);
                        SpatiT      +=   Vc*(Tp-Tm)/(2.0*FL.Grid.dx)/Tc;
                        FL.DivVel(i,j,k,{0}) = (TempoT + SpatiT);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowMixture::CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL,  EnergyTransport& ET,
                                              SpeciesTransport& ST, BoundaryConditions& BC, double dt)
{
    double CRI = (DI)?(1.0-DBL_EPSILON):(0.5-DBL_EPSILON);
    std::vector<int> dir(3);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DivVel,0,)
    {
        for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            size_t phaseIdx =  Phase.FieldsProperties[it->index].Phase;
            if(Phase.PhaseAggregateStates[phaseIdx] == AggregateStates::Gas)
            {
                if(Phase.Fractions(i,j,k)[phaseIdx]>=CRI)
                {
                    double Tc     = ET.Tx(i,j,k);
                    double MMWc;
                    double TempoT = (ET.Tx(i,j,k) - ET.TxOld(i,j,k))/dt/Tc;
                    double SpatiT = 0.0;

                    double TempoC = 0.0;
                    double SpatiC = 0.0;

                    MMWc= ST.CalculateMeanMolarMass(i, j, k, "new");
                    for(size_t n =0; n < ST.nSpecies; n++)
                    {
                        TempoC += (ST.MassFractions(i,j,k,{n})-ST.MassFractionsOld(i,j,k,{n}))/dt*MMWc/ST.MolecularWeight({n});
                    }
                    
                    for (int direction=0; direction<3; ++direction)
                    {
                        dir[0] = (direction == 0) ? 1*ET.Grid.dNx : 0;
                        dir[1] = (direction == 1) ? 1*ET.Grid.dNy : 0;
                        dir[2] = (direction == 2) ? 1*ET.Grid.dNz : 0;
                        double Vc 	= V_Mixture(i,j,k)[direction];
                        double Tp   = ET.Tx(i+dir[0], j+dir[1], k+dir[2]);
                        double Tm   = ET.Tx(i-dir[0], j-dir[1], k-dir[2]);
                        SpatiT      +=   Vc*(Tp-Tm)/(2.0*FL.Grid.dx)/Tc;
                        for(size_t n =0; n < ST.nSpecies; n++)
                        {
                            //double Cc   = ST.MassFractions(i, j, k,{n});
                            double Cp   = ST.MassFractions(i+dir[0], j+dir[1], k+dir[2],{n});
                            double Cm   = ST.MassFractions(i-dir[0], j-dir[1], k-dir[2],{n});
                            SpatiC +=   ( Vc*(Cp-Cm)/(2.0*FL.Grid.dx) ) *MMWc/ST.MolecularWeight({n}) ;
                        }
                        FL.DivVel(i,j,k,{0}) = (TempoC + TempoT + SpatiC + SpatiT);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowMixture::SettingDivergenceVelocityZero(FlowSolverLBM& FL)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DivVel,FL.DivVel.Bcells(),)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        FL.DivVel(i,j,k,{n})=0.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowMixture::SetInletVelocity(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um)
{
    if(UniformVel) SetInletUniformVelocity(FL,BC,Phase,Um);
    if(PoiseuilleVel) SetInletPoiseuilleVelocity(FL,BC,Phase,Um);
}

void FlowMixture::SetInletUniformVelocity(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um)
{
    int i=0;
    if(FL.Grid.OffsetX==0)
    for (int k = 0; k < FL.Grid.Nz; ++k)
    {
        for (int j = 0; j < FL.Grid.Ny; ++j)
        {
            if(!FL.Obstacle(i,j,k))
            {
                for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
                {
                    int ip      = i+1;
                    double rhop = FL.DensityWetting(ip,j,k,{n});
                    double rho  = FL.DensityWetting(i,j,k,{n});

                    double pp = 0.0;
                    double p  = 0.0;

                    dVector3 Mvp = {0.0, 0.0, 0.0};
                    dVector3 u   = {Um,0.0,0.0};
                    for (int ii = -FL.Grid.dNx; ii <= FL.Grid.dNx; ++ii)
                    for (int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                    for (int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                    {
                        double rrp = FL.lbPopulations(ip,j,k,{n})(ii,jj,kk);

                        pp  += rrp*FL.dP;

                        Mvp[0] += rrp*ii*FL.dm;
                        Mvp[1] += rrp*jj*FL.dm;
                        Mvp[2] += rrp*kk*FL.dm;
                        if(FL.Obstacle(i-ii,j-jj,k-kk))
                        {
                            u = {0.0,0.0,0.0};
                        }
                    }
                    p=pp;

                    dVector3 up = (Mvp * 3.0 + FL.ForceDensity(ip,j,k,{n})*FL.dt/2.0)/rhop;

                    pp +=  FL.cs2/2.0 *FL.dt* ( (up[0]*  FL.GradRho(ip,j,k,{n})[0] )
                                         +      (up[1]*  FL.GradRho(ip,j,k,{n})[1] )
                                         +      (up[2]*  FL.GradRho(ip,j,k,{n})[2] ) + rhop * FL.DivVel(ip,j,k,{n}) );
                    p  +=  FL.cs2/2.0 *FL.dt* ( (u [0]*  FL.GradRho(i ,j,k,{n})[0] )
                                         +      (u [1]*  FL.GradRho(i ,j,k,{n})[1] )
                                         +      (u [2]*  FL.GradRho(i ,j,k,{n})[2] ) + rho  * FL.DivVel(i ,j,k,{n}) );

                    dVector3 lbup = up * FL.dt/FL.Grid.dx;
                    dVector3 lbu  = u  * FL.dt/FL.Grid.dx;

                    double lbu2p  = lbup[0]*lbup[0]+lbup[1]*lbup[1]+lbup[2]*lbup[2];
                    double lbu2   = lbu [0]*lbu [0]+lbu [1]*lbu [1]+lbu [2]*lbu [2];

                    double lbtaup  = FL.nut(ip,j,k,{n})/FL.dnu/lbcs2 + 0.5;
                    double lbtau   = FL.nut(i ,j,k,{n})/FL.dnu/lbcs2 + 0.5;
                    
                    for(int ii = -FL.Grid.dNx; ii <= FL.Grid.dNx; ++ii)
                    {
                        for(int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                        {
                            for(int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                            {
                                if(!FL.Obstacle(i-ii,j-jj,k-kk))
                                {
                                    if(ii==1)
                                    {
                                        dVector3 lbFb        = FL.ForceDensity(i ,j,k,{n})/FL.df;
                                        dVector3 lbFbp       = FL.ForceDensity(ip,j,k,{n})/FL.df;

                                        dVector3 lbGradRho   = FL.GradRho(i ,j,k,{n})*FL.Grid.dx;
                                        dVector3 lbGradRhop  = FL.GradRho(ip,j,k,{n})*FL.Grid.dx;

                                        double   lbrho       = FL.DensityWetting(i ,j,k,{n})/FL.dRho;
                                        double   lbrhop      = FL.DensityWetting(ip,j,k,{n})/FL.dRho;

                                        double lbp           = p /FL.dP;
                                        double lbpp          = pp/FL.dP;

                                        double lbDivVel      = FL.DivVel(i ,j,k,{n})*FL.dt;
                                        double lbDivVelp     = FL.DivVel(ip,j,k,{n})*FL.dt;
                                        double factorp       = (1.0+1.0/(2.0*lbtaup));

                                        double CR            = (ii-lbu [0])*lbGradRho [0]+(jj-lbu [1])*lbGradRho [1]+(kk-lbu [2])*lbGradRho [2];
                                        double CRp           = (ii-lbup[0])*lbGradRhop[0]+(jj-lbup[1])*lbGradRhop[1]+(kk-lbup[2])*lbGradRhop[2];

                                        double CF            = (ii-lbu [0])*lbFb [0]+(jj-lbu [1])*lbFb [1]+(kk-lbu [2])*lbFb [2];
                                        double CFp           = (ii-lbup[0])*lbFbp[0]+(jj-lbup[1])*lbFbp[1]+(kk-lbup[2])*lbFbp[2];

                                        double lbw           = FL.lbWeights[ii+1][jj+1][kk+1];

                                        double lbcu          = ii*lbu [0] + jj*lbu [1] + kk*lbu [2];
                                        double lbcup         = ii*lbup[0] + jj*lbup[1] + kk*lbup[2];

                                        double Feq           = lbrho  * lbw *(1.0 + lbcu /lbcs2- lbu2 /(2.0*lbcs2) + lbcu *lbcu /(2.0*lbcs2*lbcs2) );
                                        double Feqp          = lbrhop * lbw *(1.0 + lbcup/lbcs2- lbu2p/(2.0*lbcs2) + lbcup*lbcup/(2.0*lbcs2*lbcs2) );
                                        double geq           = Feq *lbcs2 + (lbp  - lbrho  * lbcs2 ) * lbw; 
                                        double geqp          = Feqp*lbcs2 + (lbpp - lbrhop * lbcs2 ) * lbw; 

                                        double A             = 0.5 * (lbw*lbcs2*lbrho *lbDivVel +CR *lbcs2*(Feq /lbrho -lbw)+ CF *Feq /lbrho );
                                        double Ap            = 0.5 * (lbw*lbcs2*lbrhop*lbDivVelp+CRp*lbcs2*(Feqp/lbrhop-lbw)+ CFp*Feqp/lbrhop);

                                        double gbarp         = FL.lbPopulations(ip,j,k,{n})(ii,jj,kk);
                                        double gp            = (gbarp+1/lbtaup/2.0*geqp+Ap)/factorp;

                                        double g             = geq  + gp - geqp;
                                        FL.lbPopulations(i,j,k,{n})(ii,jj,kk) = g - 1/2.0/lbtau * (geq-g) - A;
                                        //FL.lbPopulations(i,j,k,{n})(ii,jj,kk) = geq +gbarp - geqp;
                                        FL.lbPopulations(i,j,k,{n})(ii,jj,kk) = geq  - A;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void FlowMixture::SetVelocityInletUniformBCX0(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um)
{
    int i=0;
    if(FL.Grid.OffsetX==0)
    for (int k = 0; k < FL.Grid.Nz; ++k)
    {
        for (int j = 0; j < FL.Grid.Ny; ++j)
        {
            if(!FL.Obstacle(i,j,k))
            {
                for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
                {
                    int ii=1;
                    for(int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                    {
                        for(int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                        {
                            double lbUx=Um*FL.dt/FL.Grid.dx;
                            double lbcu = ii*lbUx;
                            FL.lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                                        +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k,{n})/FL.dRho * lbcu;
                        }
                    }
                }
            }
        }
    }
}
void FlowMixture::SetCornerAtInlet(FlowSolverLBM& FL, const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.lbPopulations, 0,)
    {
        if((i+Phase.Grid.OffsetX)==0)
        {
        if(!FL.Obstacle(i,j,k) and FL.Obstacle(i,j,k-1)) //left corners
        {
            double lbUx=0.0;
            for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
            {
                for(int ii = -FL.Grid.dNx; ii <= FL.Grid.dNx; ++ii)
                {
                    for(int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                    {
                        for(int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                        {
                            if(ii==1 and kk==0)
                            {
                                double cu = ii*lbUx;
                                FL.lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                                                                    +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k,{n})/FL.dRho * cu;
                            }
                            if(ii==1 and kk==-1)
                            {
                                double cu = ii*lbUx;
                                FL.lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=FL.lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                                        +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k,{n})/FL.dRho * cu;
                            }
                        }
                    }
                }
            }
        }
        if(!FL.Obstacle(i,j,k) and FL.Obstacle(i,j,k+1))   //right corners
        {
            double lbUx=0.0;
            for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
            {
                for(int ii = -FL.Grid.dNx; ii <= FL.Grid.dNx; ++ii)
                {
                    for(int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                    {
                        for(int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                        {
                            if(ii==1 and kk==0)
                            {
                                double cu = ii*lbUx;
                                FL.lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                                                                    +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k,{n})/FL.dRho * cu;
                            }
                            if(ii==1 and kk==1)
                            {
                                double cu = ii*lbUx;
                                FL.lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=FL.lbPopulations(i,j,k,{n})(-ii,-jj,-kk)
                                        +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k,{n})/FL.dRho * cu;
                            }
                        }
                    }
                }
            }
        }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void FlowMixture::CalculatingMassFlowRate(FlowSolverLBM& FL, PhaseField& Phase, double &Sum, int i)
{
    double delY=1.0;
    double delZ=1.0;
    double dx=FL.Grid.dx;
    
    if(FL.Grid.dNy)
    {
        delY=dx;
    }
    if(FL.Grid.dNz)
    {
        delZ=dx;
    } 

    for(int j=0;j<FL.Grid.Ny;j++)
    {
        for(int k=0;k<FL.Grid.Nz;k++)
        {
            if(!FL.Obstacle(i,j,k))
            {
                Sum += V_Mixture(i,j,k)[0]*FL.DensityWetting(i,j,k,{0})*delY*delZ;
            }
        }
    }
}

void FlowMixture::SetInletPoiseuilleVelocity(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um)
{
    double Umax=1.5*Um;
    int i=0;
    double k0=0;
    double Nl=0;
    bool foundk0=false;
    if(FL.Grid.OffsetX==0)
    for (int k = 0; k < FL.Grid.Nz; ++k)
    {
        for (int j = 0; j < FL.Grid.Ny; ++j)
        {
            if(FL.Obstacle(i,j,k))
            {
                if(!foundk0)
                {
                    k0++;
                }
            }
            if(!FL.Obstacle(i,j,k))
            {
                foundk0=true;
                Nl++;
            }
        }
    }
    #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &Nl, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &k0, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif
    
    if(FL.Grid.OffsetX==0)
    for (int k = 0; k < FL.Grid.Nz; ++k)
    {
        for (int j = 0; j < FL.Grid.Ny; ++j)
        {
            if(!FL.Obstacle(i,j,k))
            {
                for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
                {
                    int ip      = i+1;
                    double rhop = FL.DensityWetting(ip,j,k,{n});
                    double rho  = FL.DensityWetting(i,j,k,{n});

                    double pp = 0.0;
                    double p  = 0.0;

                    dVector3 Mvp = {0.0, 0.0, 0.0};

                    double ku=k+FL.Grid.OffsetZ-k0;
                    double ux=4.0*Umax*ku/Nl*(1.0-ku/Nl);
                    dVector3 u   = {ux,0.0,0.0};

                    for (int ii = -FL.Grid.dNx; ii <= FL.Grid.dNx; ++ii)
                    for (int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                    for (int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                    {
                        double rrp = FL.lbPopulations(ip,j,k,{n})(ii,jj,kk);

                        pp  += rrp*FL.dP;

                        Mvp[0] += rrp*ii*FL.dm;
                        Mvp[1] += rrp*jj*FL.dm;
                        Mvp[2] += rrp*kk*FL.dm;
                        if(FL.Obstacle(i-ii,j-jj,k-kk))
                        {
                            u= {0.0,0.0,0.0};
                        }
                    }
                    p=pp;

                    dVector3 up = (Mvp * 3.0 + FL.ForceDensity(ip,j,k,{n})*FL.dt/2.0)/rhop;

                    pp +=  FL.cs2/2.0 *FL.dt* ( (up[0]*  FL.GradRho(ip,j,k,{n})[0] )
                                         +      (up[1]*  FL.GradRho(ip,j,k,{n})[1] )
                                         +      (up[2]*  FL.GradRho(ip,j,k,{n})[2] ) + rhop * FL.DivVel(ip,j,k,{n}) );
                    p  +=  FL.cs2/2.0 *FL.dt* ( (u [0]*  FL.GradRho(i ,j,k,{n})[0] )
                                         +      (u [1]*  FL.GradRho(i ,j,k,{n})[1] )
                                         +      (u [2]*  FL.GradRho(i ,j,k,{n})[2] ) + rho  * FL.DivVel(i ,j,k,{n}) );

                    dVector3 lbup = up * FL.dt/FL.Grid.dx;
                    dVector3 lbu  = u  * FL.dt/FL.Grid.dx;

                    double lbu2p  = lbup[0]*lbup[0]+lbup[1]*lbup[1]+lbup[2]*lbup[2];
                    double lbu2   = lbu [0]*lbu [0]+lbu [1]*lbu [1]+lbu [2]*lbu [2];

                    //double lbtaup  = FL.nut(ip,j,k,{n})/FL.dnu/lbcs2 + 0.5;
                    //double lbtau   = FL.nut(i ,j,k,{n})/FL.dnu/lbcs2 + 0.5;
                    
                    //int ii=1;
                    for(int ii = -FL.Grid.dNx; ii <= FL.Grid.dNx; ++ii)
                    {
                    for(int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                    {
                        for(int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                        {
                            //if(!FL.Obstacle(i-ii,j-jj,k-kk))
                            {
                            //dVector3 lbFb        = FL.ForceDensity(i ,j,k,{n})/FL.df;
                            //dVector3 lbFbp       = FL.ForceDensity(ip,j,k,{n})/FL.df;

                            //dVector3 lbGradRho   = FL.GradRho(i ,j,k,{n})*FL.Grid.dx;
                            //dVector3 lbGradRhop  = FL.GradRho(ip,j,k,{n})*FL.Grid.dx;

                            double   lbrho       = FL.DensityWetting(i ,j,k,{n})/FL.dRho;
                            double   lbrhop      = FL.DensityWetting(ip,j,k,{n})/FL.dRho;

                            double lbp           = p /FL.dP;
                            double lbpp          = pp/FL.dP;

                            //double lbDivVel      = FL.DivVel(i ,j,k,{n})*FL.dt;
                            //double lbDivVelp     = FL.DivVel(ip,j,k,{n})*FL.dt;
                            //double factorp       = (1.0+1.0/(2.0*lbtaup));

                            //double CR            = (ii-lbu [0])*lbGradRho [0]+(jj-lbu [1])*lbGradRho [1]+(kk-lbu [2])*lbGradRho [2];
                            //double CRp           = (ii-lbup[0])*lbGradRhop[0]+(jj-lbup[1])*lbGradRhop[1]+(kk-lbup[2])*lbGradRhop[2];

                            //double CF            = (ii-lbu [0])*lbFb [0]+(jj-lbu [1])*lbFb [1]+(kk-lbu [2])*lbFb [2];
                            //double CFp           = (ii-lbup[0])*lbFbp[0]+(jj-lbup[1])*lbFbp[1]+(kk-lbup[2])*lbFbp[2];

                            double lbw           = FL.lbWeights[ii+1][jj+1][kk+1];

                            double lbcu          = ii*lbu [0] + jj*lbu [1] + kk*lbu [2];
                            double lbcup         = ii*lbup[0] + jj*lbup[1] + kk*lbup[2];

                            double Feq           = lbrho  * lbw *(1.0 + lbcu /lbcs2- lbu2 /(2.0*lbcs2) + lbcu *lbcu /(2.0*lbcs2*lbcs2) );
                            double Feqp          = lbrhop * lbw *(1.0 + lbcup/lbcs2- lbu2p/(2.0*lbcs2) + lbcup*lbcup/(2.0*lbcs2*lbcs2) );
                            double geq           = Feq *lbcs2 + (lbp  - lbrho  * lbcs2 ) * lbw; 
                            double geqp          = Feqp*lbcs2 + (lbpp - lbrhop * lbcs2 ) * lbw; 

                            //double A             = 0.5 * (lbw*lbcs2*lbrho *lbDivVel +CR *lbcs2*(Feq /lbrho -lbw)+ CF *Feq /lbrho );
                            //double Ap            = 0.5 * (lbw*lbcs2*lbrhop*lbDivVelp+CRp*lbcs2*(Feqp/lbrhop-lbw)+ CFp*Feqp/lbrhop);

                            double gbarp         = FL.lbPopulations(ip,j,k,{n})(ii,jj,kk);
                            //double gp            = (gbarp+1/lbtaup/2.0*geqp+Ap)/factorp;

                            //double g             = geq  + gp - geqp;
                            //FL.lbPopulations(i,j,k,{n})(ii,jj,kk) = g - 1/2.0/lbtau * (geq-g) - A;
                            FL.lbPopulations(i,j,k,{n})(ii,jj,kk) = geq + gbarp - geqp;
                            //FL.lbPopulations(i,j,k,{n})(ii,jj,kk) = geq  - A;
                            }
                        }
                    }
                    }
                }
            }
        }
    }
}

void FlowMixture::SetVelocityInletPoiseuilleFlow(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um)
{
    double Umax=1.5*Um;
    int i=0;
    double k0=0;
    double Nl=0;
    bool foundk0=false;
    if(FL.Grid.OffsetX==0)
    for (int k = 0; k < FL.Grid.Nz; ++k)
    {
        for (int j = 0; j < FL.Grid.Ny; ++j)
        {
            if(FL.Obstacle(i,j,k))
            {
                if(!foundk0)
                {
                    k0++;
                }
            }
            if(!FL.Obstacle(i,j,k))
            {
                foundk0=true;
                Nl++;
            }
        }
    }
    #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &Nl, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &k0, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    if(FL.Grid.OffsetX==0)
    for (int k = 0; k < FL.Grid.Nz; ++k)
    {
        for (int j = 0; j < FL.Grid.Ny; ++j)
        {
            if(!FL.Obstacle(i,j,k))
            {
                for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
                {
                    int ii=1;
                    for(int jj = -FL.Grid.dNy; jj <= FL.Grid.dNy; ++jj)
                    {
                        for(int kk = -FL.Grid.dNz; kk <= FL.Grid.dNz; ++kk)
                        {
                            double ku=k+FL.Grid.OffsetZ-k0;
                            double lbUx=Umax/FL.dv;
                            double lbcu=4.0*lbUx*ku/Nl*(1.0-ku/Nl);
                            FL.lbPopulations(i-ii,j-jj,k-kk,{n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk,{n})(-ii,-jj,-kk)
                                        +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k,{n})/FL.dRho * lbcu;
                        }
                    }
                }
            }
        }
    }
}
