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
 *  Main contributors :   Reza Namdar
 *
 */

#include "ReactiveFlows/MixtureFlow.h"
#include "Containers/dVector3.h"

using namespace std;
using namespace openphase;

constexpr double lbcs2 = 1.0/3.0;                                               ///< Speed of sound squared [lattice units]

void MixtureFlow::ReadInput(string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("MixtureFlow");
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
    int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "MixtureFlow");
    LengthScale    	     = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("LengthScale"), false, 1.0);
    ReynoldsNumber    	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ReynoldsNumber"), false, 1.0);
    UPDATE_VELOCITY      = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("UPDATE_VELOCITY"), false, false);
    POISEFLOW            = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("POISEFLOW"), false, false);
    Austar               = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Austar"), false, 300.0);
    DI                   = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("DI"), false, false);
    DO_ZCORNER           = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("ZCORNER"), false, false);
    DO_YCORNER           = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("YCORNER"), false, false);
    SecOrdBB             = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("SecOrdBB"), false, false);
    lbnu                 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("lbnu"), false, 0.005);
    MaxU                 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("MaxU"), false, 1.0);
    Kp                   = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Kp"), false, 1.0);
}
void MixtureFlow::Initialize(Settings& locSettings)
{
    Grid = locSettings.Grid;
    VelocityMix. Allocate(Grid, Grid.Bcells);
    DensityMix.  Allocate(Grid, Grid.Bcells);
    GradDensity. Allocate(Grid, Grid.Bcells);
    ViscosityMix.Allocate(Grid, Grid.Bcells);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,VelocityMix,0,)
    {
        VelocityMix(i,j,k).set_to_zero();
        DensityMix(i,j,k) = 1;
        ViscosityMix(i,j,k) = 1e-5;
        GradDensity(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::SetInitial(PhaseField &Phase, FlowSolverLBM& FL, BoundaryConditions& BC)
{
    SetInitialVelocity(Phase, BC, FL);
    SetFlowProperties(Phase,FL,BC);
    LBMLimits(Phase, FL, MaxU, lbnu, FL.dt);
}

void MixtureFlow::SetInitialVelocity(PhaseField& Phase, BoundaryConditions& BC, FlowSolverLBM& FL)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,VelocityMix,VelocityMix.Bcells(),)
    if (!FL.Obstacle(i,j,k))
    {
        VelocityMix(i,j,k)={FL.U0X,FL.U0Y,FL.U0Z};
    }
    else
    {
        VelocityMix(i,j,k)={0.0,0.0,0.0};
    }

    OMP_PARALLEL_STORAGE_LOOP_END
    if(Grid.dNx) BC.SetXVector(VelocityMix);
    if(Grid.dNy) BC.SetYVector(VelocityMix);
    if(Grid.dNz) BC.SetZVector(VelocityMix);
}

void MixtureFlow::SetFlowProperties(PhaseField& Phase, FlowSolverLBM& FL, BoundaryConditions& BC)
{
    bool getfirst=false;
    InletDensity=0.0;
    InletViscosity=0.0;
    int iref=0;
    cout<<"Calculate inlet density at x node: "<<iref<<endl;
    if(Phase.Grid.OffsetX==0.0 and Phase.Grid.OffsetY==0.0 and Phase.Grid.OffsetZ==0.0)
    {
        for (int j = 0; j < Grid.Ny; j++)
        {
            for (int k = 0; k < Grid.Nz; k++)
            {
                if(!FL.Obstacle(iref,j,k))
                {
                    if(!getfirst)
                    {
                        InletDensity=DensityMix(iref,j,k);
                        InletViscosity=ViscosityMix(iref,j,k);
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
        iref = Grid.Nx/2.0;
        cout<<"Calculate inlet density at x node: "<<iref<<endl;
        if(Phase.Grid.OffsetX==0.0 and Phase.Grid.OffsetY==0.0 and Phase.Grid.OffsetZ==0.0)
        {
            for (int j = 0; j < Grid.Ny; j++)
            {
                for (int k = 0; k < Grid.Nz; k++)
                {
                    if(!FL.Obstacle(iref,j,k))
                    {
                        if(!getfirst)
                        {
                            InletDensity=DensityMix(iref,j,k);
                            InletViscosity=ViscosityMix(iref,j,k);
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
    if(Phase.Grid.OffsetX==0.0 and Phase.Grid.OffsetY==0.0 and Phase.Grid.OffsetZ==0.0)
    {
        size_t Ni=Grid.Nx;
        size_t Nj=Grid.Ny;
        size_t Nk=Grid.Nz;
        if(FL.U0X!=0) Ni=1; 
        if(FL.U0Y!=0) Nj=1; 
        if(FL.U0Z!=0) Nk=1; 
        for (size_t i = 0; i < Ni; i++)
        {
            for (size_t j = 0; j < Nj; j++)
            {
                for (size_t k = 0; k < Nk; k++)
                {
                    if(!FL.Obstacle(i,j,k))
                        nodes++;
                }
            }
        }
    }
    #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &nodes, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif
    InletArea = nodes * Grid.dx;
    SetInitialPopulations(FL,BC);

    InletMassFlowRate=InletDensity*InletArea*FL.U0X;

    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetY==0 and Phase.Grid.OffsetZ==0)
    {
        cout<<"Inlet area (m^2): "<<InletArea<<endl;
        cout<<"Inlet mass flow (kg/s): "<<InletMassFlowRate<<endl;
    }

    double L=1;
    if(FL.U0X!=0) L = Grid.TotalNx * Grid.dx;
    if(FL.U0Y!=0) L = Grid.TotalNy * Grid.dx;
    if(FL.U0Z!=0) L = Grid.TotalNz * Grid.dx;
    Kp = 0.2 * sqrt(FL.cs2) / L;
    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetY==0 and Phase.Grid.OffsetZ==0)
    {
        cout<<"Non-reflecting pressure constant Kp: "<<Kp<<endl;
    }
    //FL.U0X = InletMassFlowRate/InletDensity/InletArea;
}

void MixtureFlow::SetInitialPopulations(FlowSolverLBM& FL, BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.HydroPressure,FL.HydroPressure.Bcells(),)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        FL.HydroPressure (i,j,k,{n}) = FL.Poutlet;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityMix,DensityMix.Bcells(),)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        dVector3 lbvel      = VelocityMix(i,j,k)/FL.dv;
        dVector3 lbFb       = FL.ForceDensity(i,j,k,{n})/FL.df;
        dVector3 lbGradRho  = GradDensity(i,j,k)*Grid.dx;
        double lbnu         = ViscosityMix(i,j,k)/FL.dnu;
        double lbrho        = DensityMix(i,j,k)/FL.dRho;
        double lbph         = FL.HydroPressure(i,j,k,{n})/FL.dP;
        double lbDivVel     = FL.DivVel(i,j,k,{n})*FL.dt;
        FL.lbPopulations(i,j,k,{n}) = FL.EquilibriumDistributionTC( lbvel, lbph, lbnu,  lbDivVel,  lbFb, FL.lbWeights, lbrho, lbGradRho, n);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    BC.SetXVector(FL.lbPopulations);
    BC.SetYVector(FL.lbPopulations);
    BC.SetZVector(FL.lbPopulations);
}

void MixtureFlow::WriteVTK(Settings& locSettings, FlowSolverLBM& FL, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"Flow_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t){"Velocity",  [this](int i,int j,int k){return VelocityMix(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t){"Pressure",  [this, FL](int i,int j,int k){return FL.HydroPressure(i,j,k,{0});}});
    VTK::Write(Filename, locSettings, ListOfFields);
}

void MixtureFlow::WriteVTKDivergenceVelocity(Settings& locSettings, FlowSolverLBM& FL, const int tStep, const int precision)
{
    const size_t MixtureFlow = 0 ;
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"DivVel_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t){"DivU",  [this, &FL](int i,int j,int k){return FL.DivVel(i,j,k,{MixtureFlow});}});
    VTK::Write(Filename, locSettings, ListOfFields);
}

void MixtureFlow::LBMLimits(PhaseField& Phase, FlowSolverLBM& FL, double Umax, double lbnu, double dt)
{

    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetY==0 and Phase.Grid.OffsetZ==0)
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
    }
}

void MixtureFlow::Collision(FlowSolverLBM& FL, BoundaryConditions &BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityMix,0,)
    if (!FL.Obstacle(i,j,k))
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        dVector3 lbvel      = VelocityMix(i,j,k)/FL.dv;
        dVector3 lbFb       = FL.ForceDensity(i,j,k,{n})/FL.df;
        dVector3 lbGradRho  = GradDensity(i,j,k)*Grid.dx;
        double lbu2         = lbvel[0]*lbvel[0]+lbvel[1]*lbvel[1]+lbvel[2]*lbvel[2];
        double lbnu         = ViscosityMix(i,j,k)/FL.dnu;
        double lbrho        = DensityMix(i,j,k)/FL.dRho;
        double lbph         = FL.HydroPressure(i,j,k,{n})/FL.dP;
        double lbDivVel     = FL.DivVel(i,j,k,{n})*FL.dt;
        double lbtau        = lbnu/lbcs2 + 0.5;
        double factor       = (1.0-1.0/(2.0*lbtau));

        for(int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
        for(int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
        for(int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
        {
            double CR = (ii-lbvel[0])*lbGradRho[0]+(jj-lbvel[1])*lbGradRho[1]+(kk-lbvel[2])*lbGradRho[2];
            double CF = (ii-lbvel[0])*lbFb[0]+(jj-lbvel[1])*lbFb[1]+(kk-lbvel[2])*lbFb[2];
            double lbw  = FL.lbWeights[ii+1][jj+1][kk+1];
            double lbcu = ii*lbvel[0] + jj*lbvel[1] + kk*lbvel[2];
            double Feq  = lbrho * lbw *(1.0 + lbcu/lbcs2- lbu2/(2.0*lbcs2) + lbcu*lbcu/(2.0*lbcs2*lbcs2) );
            double geq  = Feq*lbcs2 + (lbph - lbrho * lbcs2 ) * lbw; 

            double geqb = geq - 0.5 * (lbrho * lbcs2 * lbw * lbDivVel * factor  + CR*(Feq/lbrho-lbw)*lbcs2*factor + CF*Feq/lbrho);
            double gb  = FL.lbPopulations(i,j,k,{n})(ii,jj,kk);
            FL.lbPopulations(i,j,k,{n})(ii,jj,kk) = gb + 1.0/lbtau * (geqb - gb) + lbcs2 * lbw * lbrho * lbDivVel * factor + CR*(Feq/lbrho-lbw)*lbcs2*factor + CF*Feq/lbrho;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if(Grid.dNx) BC.SetXVector(FL.lbPopulations);
    if(Grid.dNy) BC.SetYVector(FL.lbPopulations);
    if(Grid.dNz) BC.SetZVector(FL.lbPopulations);
}

void MixtureFlow::CalculateVelocityAndPressure(FlowSolverLBM& FL)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,VelocityMix,0,)
    //if(!FL.Obstacle(i,j,k) or FL.ObstacleAppeared(i,j,k))
    if(!FL.Obstacle(i,j,k))
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        FL.HydroPressure(i,j,k,{n}) = 0.0;
        VelocityMix(i,j,k).set_to_zero();

        for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
        for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
        for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
        {
            double rr = FL.lbPopulations(i,j,k,{n})(ii,jj,kk);
            FL.HydroPressure  (i,j,k,{n})    += rr*FL.dP;
            VelocityMix(i,j,k)[0] += rr*ii;
            VelocityMix(i,j,k)[1] += rr*jj;
            VelocityMix(i,j,k)[2] += rr*kk;
        }
        
        VelocityMix(i,j,k) = (VelocityMix(i,j,k)/FL.cs2*pow(Grid.dx,3.0)/pow(FL.dt,3.0) + FL.ForceDensity(i,j,k,{n})*FL.dt/2.0)/DensityMix(i,j,k);
        FL.HydroPressure (i,j,k,{n}) +=  FL.cs2*FL.dt/2.0 * ( (VelocityMix(i,j,k)[0]*GradDensity(i,j,k)[0] +
                                                      VelocityMix(i,j,k)[1]*GradDensity(i,j,k)[1] + 
                                                      VelocityMix(i,j,k)[2]*GradDensity(i,j,k)[2]) +
                                                      DensityMix(i,j,k) * FL.DivVel(i,j,k,{n}));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::Propagation(FlowSolverLBM& FL, PhaseField& Phase, SolidBody& SB, BoundaryConditions& BC)
{
    if (SecOrdBB and !DI) 
    {
        PropagationSecondOrderBB(Phase, FL, SB, BC);
    } 
    else 
    {
        FL.Propagation(Phase, BC);
    }
    if(Grid.dNx) BC.SetXVector(FL.lbPopulations);
    if(Grid.dNy) BC.SetYVector(FL.lbPopulations);
    if(Grid.dNz) BC.SetZVector(FL.lbPopulations);
}

void MixtureFlow::PropagationSecondOrderBB(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.lbPopulationsTMP,0,)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        FL.lbPopulationsTMP(i,j,k,{n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.lbPopulations,0,)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        double lbDensityChange = 0.0;
        for(int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
        for(int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
        for(int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
        {
            if (FL.Obstacle(i, j, k))
            {
                FL.lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) = FL.lbPopulations(i,j,k,{n})(ii,jj,kk);
            }
            else if (FL.Obstacle(i-ii, j-jj, k-kk))
            {
                FL.lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) = SecondOrderBounceBack(i, j, k, ii, jj, kk, n, Phase, FL, SB, BC, lbDensityChange);
            }
            else
            {
                FL.lbPopulationsTMP(i,j,k,{n})(ii,jj,kk) = FL.lbPopulations(i-ii, j-jj, k-kk,{n})(ii,jj,kk);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.lbPopulations,0,)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        FL.lbPopulations(i,j,k,{n}) = FL.lbPopulationsTMP(i,j,k,{n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double MixtureFlow::SecondOrderBounceBack(const int i, const int j, const int k,
        const int ii, const int jj, const int kk, const size_t n,
        PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB,
        const BoundaryConditions& BC, double& lbDensityChange)
{
    double NewPopulation = 0.0;

    // ------------------------------------------------------------
    // 1) Geometry setup
    // ------------------------------------------------------------
    // Lattice node positions
    dVector3 Xf = { double(i),      double(j),      double(k)     };      // fluid node
    dVector3 Xs = { double(i - ii), double(j - jj), double(k - kk) };     // neighbor (solid) node

    dVector3 ql = SB.DistanceToInterface(Xf,Xs);
    // Periodic link vector from Xf to Xs
    dVector3 Xsf;
    SB.CalculateDistancePeriodic(Xs, Xf, Xsf, Grid.Nx, Grid.Ny, Grid.Nz);
    const double l = Xsf.length();        // physical link length

    // Normalize q to [0,1] for bounce-back rules
    double q = (l > 1e-30) ? (ql.length() / l) : 0.0;

    // ------------------------------------------------------------
    // 2) Second-order bounce-back using the normalized q
    // ------------------------------------------------------------
    if (q < 0.5)
    {
        NewPopulation =
            2.0 * q * FL.lbPopulations(i,j,k,{n})(-ii, -jj, -kk) +
            (1.0 - 2.0 * q) * FL.lbPopulations(i+ii,j+jj,k+kk,{n})(-ii, -jj, -kk);
    }
    else // q >= 0.5
    {
        NewPopulation =
            1.0 / (2.0 * q) * FL.lbPopulations(i, j, k, {n})(-ii, -jj, -kk) +
            (2.0 * q - 1.0) / (2.0 * q) * FL.lbPopulations(i, j, k, {n})(ii, jj, kk);
    }

    return NewPopulation;
}

void MixtureFlow::ApplyForces(PhaseField& Phase, FlowSolverLBM& FL)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.ForceDensity,0,)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        FL.ForceDensity(i,j,k,{n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(DI) CalculateForceDragbyPF(Phase,FL,Austar);
    if (FL.Do_Gravity) CalculateForceGravity(Phase, FL);
}

void MixtureFlow::CalculateForceGravity(PhaseField& Phase, FlowSolverLBM& FL)
{
    // Calculate gravity on solid bodies
    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if (Phase.FieldsProperties[idx].State == AggregateStates::Solid and
        Phase.FieldsProperties[idx].Mobile)
    {
        Phase.FieldsProperties[idx].Acm += FL.GA;
    }

    // Calculate gravity on fluid nodes
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.ForceDensity,0,)
    if (not FL.Obstacle(i,j,k))
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        FL.ForceDensity(i,j,k,{n}) += FL.GA * DensityMix(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::CalculateForceDragbyPF(PhaseField& Phase, FlowSolverLBM& FL, double hHyd )
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityMix,0,)
    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
    {
        if (Phase.Fields(i,j,k).interface())
        {
            const double mu = DensityMix(i,j,k)*ViscosityMix(i,j,k);  //Dynamic viscosity [kg/m/s]
            dVector3 vel = VelocityMix(i,j,k);
            double SF = FL.SolidFraction(i,j,k,Phase);
            double CF = hHyd * mu /pow(Phase.Grid.iWidth*Grid.dx,2.0) * pow(SF,3)*(1.0-SF);
            FL.ForceDensity(i,j,k,{n}) +=  vel *(-CF);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::CalculateDensityGradient(FlowSolverLBM& FL, Energy& EN)
{
    constexpr double eps = 1e-10;
    const double invDx  = 1.0 / Grid.dx;
    const double inv2Dx = 0.5 * invDx;

    // Active directions only
    int activeDirs[3]; int nActive = 0;
    if (Grid.dNx) activeDirs[nActive++] = 0;
    if (Grid.dNy) activeDirs[nActive++] = 1;
    if (Grid.dNz) activeDirs[nActive++] = 2;

    auto vanleer_phi = [](double r) -> double 
    {
        const double ar = fabs(r);
        return (r + ar) / (1.0 + ar);
    };

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityMix,0,)
    {
        if (FL.Obstacle(i,j,k)) continue;

        // Avoid repeated enum comparisons in inner loops
        const auto scheme = EN.AdvScheme;

        for (int a = 0; a < nActive; ++a)
        {
            const int direction = activeDirs[a];

            int di=0, dj=0, dk=0;
            if (direction == 0) di = Grid.dNx;
            else if (direction == 1) dj = Grid.dNy;
            else dk = Grid.dNz;

            const int ip = i + di, jp = j + dj, kp = k + dk;
            const int im = i - di, jm = j - dj, km = k - dk;

            const double v = VelocityMix(i,j,k)[direction];

            for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
            {
                const double Rho  = DensityMix(i,j,k);
                const double Rhop = DensityMix(ip,jp,kp);
                const double Rhom = DensityMix(im,jm,km);

                if (scheme == AdvectionSchemes::Upwind)
                {
                    GradDensity(i,j,k)[direction] = (v >= 0.0) ? ((Rho - Rhom) * invDx)
                                                                 : ((Rhop - Rho) * invDx);
                }
                else if (scheme == AdvectionSchemes::Central)
                {
                    GradDensity(i,j,k)[direction] = (Rhop - Rhom) * inv2Dx;
                }
                else if (scheme == AdvectionSchemes::VanLeer)
                {
                    if (v >= 0.0)
                    {
                        const int imm = i - 2*di, jmm = j - 2*dj, kmm = k - 2*dk;
                        const double Rhomm = DensityMix(imm,jmm,kmm);

                        const double r   = (Rho  - Rhom)  / (Rhop - Rho  + eps);
                        const double rm  = (Rhom - Rhomm) / (Rho  - Rhom + eps);

                        const double Rhohp = Rho  + 0.5 * vanleer_phi(r)  * (Rhop - Rho);
                        const double Rhohm = Rhom + 0.5 * vanleer_phi(rm) * (Rho  - Rhom);

                        GradDensity(i,j,k)[direction] = (Rhohp - Rhohm) * invDx;
                    }
                    else
                    {
                        const int ipp = i + 2*di, jpp = j + 2*dj, kpp = k + 2*dk;
                        const double Rhopp = DensityMix(ipp,jpp,kpp);

                        const double r   = (Rhopp - Rhop) / (Rhop - Rho  + eps);
                        const double rm  = (Rhop  - Rho)  / (Rho  - Rhom + eps);

                        const double Rhohp = Rhop - 0.5 * vanleer_phi(r)  * (Rhop - Rho);
                        const double Rhohm = Rho  - 0.5 * vanleer_phi(rm) * (Rho  - Rhom);

                        GradDensity(i,j,k)[direction] = (Rhohp - Rhohm) * invDx;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::UpdateIncomingVelocityValue(FlowSolverLBM& FL, double FuelConsumptionRate, double FuelMF)
{
    FL.U0X  = (FuelConsumptionRate) / (FuelMF * InletDensity * InletArea);
}

void MixtureFlow::SetVanishedObstacleNodes(PhaseField &Phase, FlowSolverLBM& FL, SolidBody &SB, Energy& EN, Species& SP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.ObstacleVanished,FL.ObstacleVanished.Bcells(),)
    {
        if(FL.ObstacleVanished(i,j,k))
        {
            for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
            {
                //constexpr double kTiny = 1e-30;
                const double D  = SB.nDist;
                //const double dx = Grid.dx;
                const dVector3 n0 = SB.NormField(i,j,k);
                const dVector3 Xf = {double(i),double(j),double(k)};
                const dVector3 Xb  = Xf - n0 * SB.DistanceField(i,j,k);
                const dVector3 Xf1 = Xb + n0 * D;               // fluid-side normal sample

                const dVector3 Xfb = Xb + n0 * SB.nFallBack ;    
                const dVector3 ufb  = VelocityMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
                const dVector3 Ffb  = FL.ForceDensity(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ()),{n}); 
                const dVector3 drhofb  = GradDensity(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
                const double nufb  = ViscosityMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
                const double rhofb  = DensityMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
                const double pfb  = FL.HydroPressure(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ()),{n}); 
                const double dvelfb  = FL.DivVel(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ()),{n}); 

                auto keepFluidCorner = [&](int ii, int jj, int kk) -> bool 
                {
                    return !FL.Obstacle(ii, jj, kk);
                };

                EN.T(i,j,k) = EN.InterfaceTemperature(Phase,FL,SB,i,j,k);
                SP.MassFrac(i,j,k) = SP.InterfaceSpecies(FL,SB,i,j,k);

                const dVector3 vel = SB.MaskedTrilinearVector3(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk) { return VelocityMix(ii, jj, kk); }, keepFluidCorner,ufb); 

                const dVector3 Fb = SB.MaskedTrilinearVector3(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk) { return FL.ForceDensity(ii, jj, kk,{n}); }, keepFluidCorner,Ffb);  

                const dVector3 GradRho  = SB.MaskedTrilinearVector3(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk) { return GradDensity(ii, jj, kk); }, keepFluidCorner,drhofb);

                const double nu = SB.MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk){ return ViscosityMix(ii, jj, kk); }, keepFluidCorner,nufb);

                const double rho = SB.MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk){ return DensityMix(ii, jj, kk); }, keepFluidCorner,rhofb);

                const double ph = SB.MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk){ return FL.HydroPressure(ii, jj, kk,{n}); }, keepFluidCorner,pfb);

                const double DivVel = SB.MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk){ return FL.DivVel(ii, jj, kk,{n}); }, keepFluidCorner,dvelfb);

                const dVector3 lbvel      = vel/FL.dv;
                const dVector3 lbFb       = Fb/FL.df;
                const dVector3 lbGradRho  = GradRho*Grid.dx;
                const double lbnu         = nu/FL.dnu;
                const double lbrho        = rho/FL.dRho;
                const double lbph         = ph/FL.dP;
                const double lbDivVel     = DivVel*FL.dt;
                FL.lbPopulations(i,j,k,{n}) = FL.EquilibriumDistributionTC( lbvel, lbph, lbnu,  lbDivVel,  lbFb, FL.lbWeights, lbrho, lbGradRho, n);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::SetAppeardObstacleNodes(FlowSolverLBM& FL, SolidBody &SB)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.ObstacleAppeared,FL.ObstacleAppeared.Bcells(),)
    {
        if(FL.ObstacleAppeared(i,j,k))
        {
            for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
            {
                const double nu = ViscosityMix(i, j, k);
                const double rho = DensityMix(i, j, k);

                const dVector3 lbvel = dVector3::ZeroVector();
                const dVector3 lbFb  = dVector3::ZeroVector();
                const dVector3 lbGradRho = dVector3::ZeroVector();
                const double lbnu         = nu/FL.dnu;
                const double lbrho        = rho/FL.dRho;
                const double lbph         = 0.0;
                const double lbDivVel     = 0.0;
                FL.lbPopulations(i,j,k,{n}) = FL.EquilibriumDistributionTC( lbvel, lbph, lbnu,  lbDivVel,  lbFb, FL.lbWeights, lbrho, lbGradRho, n);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::DetectObstacles(FlowSolverLBM& FL, const PhaseField& Phase, BoundaryConditions &BC, bool DI)
{
    double CRI = (DI)?(1.0):(0.5);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityMix,DensityMix.Bcells(),)
    {
        if (FL.SolidFraction(i,j,k,Phase) >= CRI)
        {
            if(!FL.Obstacle(i,j,k))
            {
                FL.ObstacleAppeared(i,j,k) = true;
            }
            else
            {
                FL.ObstacleAppeared(i,j,k) = false;
            }
            FL.Obstacle(i,j,k) = true;
        }
        else
        {
            if(FL.Obstacle(i,j,k)) 
            {
                FL.ObstacleVanished(i,j,k) = true;
            }
            else
            {
                FL.ObstacleVanished(i,j,k) = false;
            }
            FL.Obstacle(i,j,k) = false;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (Grid.dNx) BC.SetX(FL.Obstacle);
    if (Grid.dNy) BC.SetY(FL.Obstacle);
    if (Grid.dNz) BC.SetZ(FL.Obstacle);
}

void MixtureFlow::CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL, Energy& EN, SolidBody& SB, double dt)
{
    constexpr double alpha_min = 0.1;
    const size_t FlowMix = 0;

    //const double CRI = (DI) ? 1.0 : 0.5;
    const double inv2Dx = 0.5 / Grid.dx;

    // Active directions only
    int activeDirs[3]; int nActive = 0;
    if (Grid.dNx) activeDirs[nActive++] = 0;
    if (Grid.dNy) activeDirs[nActive++] = 1;
    if (Grid.dNz) activeDirs[nActive++] = 2;

    auto unit_or_fallback = [](dVector3 n, const dVector3& fallback) -> dVector3 {
        const double nl = n.length();
        if (nl > 1e-30) n *= (1.0 / nl);
        else n = fallback;
        return n;
    };

        // Neighbor sampling (OLD) with obstacle reconstruction
    auto SampleNeighborNew = [&](int i,int j,int k, int in,int jn,int kn,
                                 const double Tc,
                                 double& Tn)
    {
        if (!FL.Obstacle(in,jn,kn))
        {
            Tn = EN.TOld(in,jn,kn);
            return;
        }

        dVector3 Xf{double(i),double(j),double(k)};
        dVector3 Xs{double(in),double(jn),double(kn)};
        dVector3 dist = SB.DistanceToInterface(Xf, Xs);

        const double alpha = std::clamp(dist.length(), alpha_min, 1.0);

        dVector3 nf = SB.NormField(i,j,k);
        dVector3 ns = SB.NormField(in,jn,kn);
        dVector3 n  = unit_or_fallback(nf + ns, nf);

        const double Tw = EN.SurfaceTemperature(Phase, FL, SB, Xf + dist, n, "new");
        Tn = (Tw - Tc) / alpha + Tc;
    };

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DivVel,0,)
    {
        FL.DivVel(i,j,k,{FlowMix}) = 0.0;
        if (FL.Obstacle(i,j,k)) continue;

        const double Tc = EN.T(i,j,k);
        if (Tc <= 0.0) continue;

        // Temporal temperature term
        const double TempoT = (Tc - EN.TOld(i,j,k)) / (dt * Tc);

        double SpatiT = 0.0;

        for (int a = 0; a < nActive; ++a)
        {
            const int direction = activeDirs[a];
            int di=0, dj=0, dk=0;
            if (direction == 0) di = Grid.dNx;
            else if (direction == 1) dj = Grid.dNy;
            else dk = Grid.dNz;

            const int ip=i+di, jp=j+dj, kp=k+dk;
            const int im=i-di, jm=j-dj, km=k-dk;

            const double Vc = VelocityMix(i,j,k)[direction];

            double Tp=Tc, Tm=Tc;

            SampleNeighborNew(i,j,k, ip,jp,kp, Tc, Tp);
            SampleNeighborNew(i,j,k, im,jm,km, Tc, Tm);

            SpatiT += Vc * (Tp - Tm) * inv2Dx / Tc;
        }

        FL.DivVel(i,j,k,{FlowMix}) = (TempoT + SpatiT);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MixtureFlow::CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL, Energy& EN,
                                                        Species& SP, SolidBody& SB, double dt)
{
    constexpr double alpha_min = 0.1;
    const size_t FlowMix = 0;

    //const double CRI = (DI) ? 1.0 : 0.5;
    const double inv2Dx = 0.5 / Grid.dx;

    // Active directions only
    int activeDirs[3]; int nActive = 0;
    if (Grid.dNx) activeDirs[nActive++] = 0;
    if (Grid.dNy) activeDirs[nActive++] = 1;
    if (Grid.dNz) activeDirs[nActive++] = 2;

    auto unit_or_fallback = [](dVector3 n, const dVector3& fallback) -> dVector3 {
        const double nl = n.length();
        if (nl > 1e-30) n *= (1.0 / nl);
        else n = fallback;
        return n;
    };

        // Neighbor sampling (OLD) with obstacle reconstruction
    auto SampleNeighborNew = [&](int i,int j,int k, int in,int jn,int kn,
                                 const double Tc,
                                 const Tensor<double,1>& Yc,
                                 double& Tn,
                                 Tensor<double,1>& Yn)
    {
        if (!FL.Obstacle(in,jn,kn))
        {
            Tn = EN.TOld(in,jn,kn);
            Yn = SP.MassFracOld(in,jn,kn);
            return;
        }

        dVector3 Xf{double(i),double(j),double(k)};
        dVector3 Xs{double(in),double(jn),double(kn)};
        dVector3 dist = SB.DistanceToInterface(Xf, Xs);

        const double alpha = std::clamp(dist.length(), alpha_min, 1.0);

        dVector3 nf = SB.NormField(i,j,k);
        dVector3 ns = SB.NormField(in,jn,kn);
        dVector3 n  = unit_or_fallback(nf + ns, nf);

        const double Tw = EN.SurfaceTemperature(Phase, FL, SB, Xf + dist, n, "new");
        Tn = (Tw - Tc) / alpha + Tc;

        Tensor<double,1> Yw = SP.SurfaceSpecies(FL, SB, Xf + dist, n, "new");
        for (size_t s=0; s<SP.nSpecies; ++s)
            Yn[s] = Yc[s] + (Yw[s] - Yc[s]) / alpha;
    };

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DivVel,0,)
    {
        FL.DivVel(i,j,k,{FlowMix}) = 0.0;
        if (FL.Obstacle(i,j,k)) continue;

        const double Tc = EN.T(i,j,k);
        if (Tc <= 0.0) continue;

        // Temporal temperature term
        const double TempoT = (Tc - EN.TOld(i,j,k)) / (dt * Tc);

        // Mean molar mass from NEW composition
        const double MMWc = SP.CalculateMeanMolarMass(i, j, k, "new");

        double TempoC = 0.0;
        for (size_t s = 0; s < SP.nSpecies; ++s)
        {
            TempoC += (SP.MassFrac(i,j,k,{s}) - SP.MassFracOld(i,j,k,{s})) / dt
                      * MMWc / SP.MwSp({s});
        }

        double SpatiT = 0.0;
        double SpatiC = 0.0;

        const Tensor<double,1> Yc = SP.MassFrac(i,j,k);

        for (int a = 0; a < nActive; ++a)
        {
            const int direction = activeDirs[a];
            int di=0, dj=0, dk=0;
            if (direction == 0) di = Grid.dNx;
            else if (direction == 1) dj = Grid.dNy;
            else dk = Grid.dNz;

            const int ip=i+di, jp=j+dj, kp=k+dk;
            const int im=i-di, jm=j-dj, km=k-dk;

            const double Vc = VelocityMix(i,j,k)[direction];

            Tensor<double,1> Yp({SP.nSpecies}), Ym({SP.nSpecies});
            double Tp=Tc, Tm=Tc;

            SampleNeighborNew(i,j,k, ip,jp,kp, Tc, Yc, Tp, Yp);
            SampleNeighborNew(i,j,k, im,jm,km, Tc, Yc, Tm, Ym);

            SpatiT += Vc * (Tp - Tm) * inv2Dx / Tc;

            for (size_t s = 0; s < SP.nSpecies; ++s)
            {
                SpatiC += Vc * (Yp[s] - Ym[s]) * inv2Dx * MMWc / SP.MwSp({s});
            }
        }

        FL.DivVel(i,j,k,{FlowMix}) = (TempoC + TempoT + SpatiC + SpatiT);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


void MixtureFlow::SetPressureBoundary(FlowSolverLBM& FL, double Pout, PressureOutletMode mode, const dVector3& A,  const dVector3& B, int n)
{
    if (n != 1 && n != -1) return;

    const int ia = A.getX(), ja = A.getY(), ka = A.getZ();
    const int ib = B.getX(), jb = B.getY(), kb = B.getZ();

    const bool xConst = (ia == ib);
    const bool yConst = (ja == jb);
    const bool zConst = (ka == kb);

    // Must be exactly one axis-aligned boundary face
    if ((int)xConst + (int)yConst + (int)zConst != 1) return;

    const int Nx = Grid.Nx, Ny = Grid.Ny, Nz = Grid.Nz;
    const double dx = Grid.dx;
    const double dt = FL.dt;
    const double cs2 = lbcs2*dx*dx/dt/dt;
    const double cs = std::sqrt(cs2);

    auto clamp = [](int v, int lo, int hi){ return std::max(lo, std::min(hi, v)); };

    auto owns_face = [&](int axis)->bool
    { 
        if (axis == 0) return ( (n > 0 && Grid.OffsetX == 0) || (n < 0 && Grid.OffsetX + Nx == Grid.TotalNx) );
        if (axis == 1) return ( (n > 0 && Grid.OffsetY == 0) || (n < 0 && Grid.OffsetY + Ny == Grid.TotalNy) );
        return                ( (n > 0 && Grid.OffsetZ == 0) || (n < 0 && Grid.OffsetZ + Nz == Grid.TotalNz) );
    };

    auto tangential_indices = [&](int axis, int& t1, int& t2)
    {
        if (axis == 0) { t1 = 1; t2 = 2; return; }
        if (axis == 1) { t1 = 0; t2 = 2; return; }
        /* axis==2 */   t1 = 0; t2 = 1;
    };

    // -------------------------
    // Kernel at a single boundary node
    // -------------------------
    auto apply_at_node = [&](int i, int j, int k, int axis)
    {
        if (FL.Obstacle(i,j,k)) return;

        int im=i, jm=j, km=k;
        int imm=i, jmm=j, kmm=k;

        if (axis == 0) { im = i + n;  imm = i + 2*n; }
        if (axis == 1) { jm = j + n;  jmm = j + 2*n; }
        if (axis == 2) { km = k + n;  kmm = k + 2*n; }

        // Need at least first interior for both modes.
        if (im < 0 || im >= Nx || jm < 0 || jm >= Ny || km < 0 || km >= Nz) return;
        if (FL.Obstacle(im, jm, km)) return;

        int t1, t2;
        tangential_indices(axis, t1, t2);

        for (size_t nf = 0; nf < FL.N_Fluid_Comp; ++nf)
        {
            // Always needed fields (boundary + first interior)
            const double Pbold  = FL.HydroPressure(i,  j,  k,  {nf});
            const double Pbmold = FL.HydroPressure(im, jm, km, {nf});

            const double uN_b  = VelocityMix(i,  j,  k )[axis];
            const double uN_m  = VelocityMix(im, jm, km)[axis];

            const double uT1_b = VelocityMix(i,  j,  k )[t1];
            const double uT1_m = VelocityMix(im, jm, km)[t1];

            const double uT2_b = VelocityMix(i,  j,  k )[t2];
            const double uT2_m = VelocityMix(im, jm, km)[t2];

            const double rhob  = DensityMix(i,  j,  k);
            const double rhobm = DensityMix(im, jm, km);
    
            // --------------------------------------------------
            // 1) Compute pm, um in the first interior cell (im)
            //    (shared for both modes)
            // --------------------------------------------------
            double pm = 0.0;
            dVector3 Mvm = {0.0, 0.0, 0.0};

            for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
            {
                const double rrm = FL.lbPopulations(im, jm, km, {nf})(ii,jj,kk);
                pm     += rrm * FL.dP;
                Mvm[0] += rrm * ii;
                Mvm[1] += rrm * jj; 
                Mvm[2] += rrm * kk;
            }
        
            dVector3 um = (Mvm/cs2*pow(dx,3.0)/pow(dt,3.0) + FL.ForceDensity(im,jm,km,{nf})*dt/2.0)/rhobm;
            pm +=  cs2*dt/2.0 * 
                ( (um[0]*GradDensity(im,jm,km)[0] 
                +  um[1]*GradDensity(im,jm,km)[1] 
                +  um[2]*GradDensity(im,jm,km)[2]) + rhobm * FL.DivVel(im,jm,km,{nf}));

            // --------------------------------------------------
            // 2) Choose boundary target (p_bn, u_bn) by mode
            // --------------------------------------------------
            double p_bn = Pbold;
            double uN_bn = uN_b, uT1_bn = uT1_b, uT2_bn = uT2_b;

            if (mode == PressureOutletMode::NonReflecting)
            {
                // Need second interior only if fallback triggers, but still must be in bounds for safe access.
                if (imm < 0 || imm >= Nx || jmm < 0 || jmm >= Ny || kmm < 0 || kmm >= Nz) return;

                double L1 = Kp * (Pbold - Pout);
                //double L1 = (uN_b - cs) * ( (Pbold - Pbmold)/dx - rhob * cs * (uN_b - uN_m)/dx );
                double L3 = uN_b * (uT1_b - uT1_m) / dx;
                double L4 = uN_b * (uT2_b - uT2_m) / dx;
                double L5 = (uN_b + cs) * ( (Pbold - Pbmold)/dx + rhob * cs * (uN_b - uN_m)/dx );

                //uN_bn  = uN_b  - (L5 - L1) * dt / (2.0 * cs * rhob);
                uT1_bn = uT1_b - L3 * dt;
                uT2_bn = uT2_b - L4 * dt;
                //p_bn   = Pbold - (L5 + L1) * dt / 2.0;
                p_bn   = pm + (L5/(uN_b+cs) + L1/(uN_b-cs)) * dx / 2.0;
                uN_bn   = um[0] + (L5/(uN_b+cs) - L1/(uN_b-cs)) * dx / (2.0 * cs * rhob);
                //cin.get();
                
                if (uN_m < 0.0)
                {
                    const double Pbmmold = FL.HydroPressure(imm, jmm, kmm, {nf});
                    const double uN_mm  = VelocityMix(imm, jmm, kmm)[axis];
                    const double uT1_mm = VelocityMix(imm, jmm, kmm)[t1];
                    const double uT2_mm = VelocityMix(imm, jmm, kmm)[t2];

                    L1 = Kp * (Pbmold - Pout);
                    L3 = uN_m * (uT1_m - uT1_mm) / dx;
                    L4 = uN_m * (uT2_m - uT2_mm) / dx;
                    L5 = (uN_m + cs) * ( (Pbmold - Pbmmold)/dx + rhobm * cs * (uN_m - uN_mm)/dx );

                    p_bn   = Pbmold - (L5 + L1) * dt / 2.0;
                    uN_bn  = uN_m   - (L5 - L1) * dt / (2.0 * cs * rhobm);
                    uT1_bn = uT1_m  - L3 * dt;
                    uT2_bn = uT2_m  - L4 * dt;
                }
                
            }
            else // FixedPressure
            {
                // Fixed pressure: set p, take velocity from the first interior cell for stability.
                // (This avoids extra characteristic work and is the common “do-nothing velocity” choice.)
                p_bn  = Pout;

                uN_bn  = um[0];
                uT1_bn = um[1];
                uT2_bn = um[2];

                // Optional safety: prevent inflow at an outlet (uncomment if desired)
                // if (uN_bn < 0.0) uN_bn = 0.0;
            }
            dVector3 lbum = um*dt/dx;
            dVector3 lbub = {uN_bn*dt/dx, uT1_bn*dt/dx, uT2_bn*dt/dx};

            // --------------------------------------------------
            // 3) Hoist per-node constants for incoming loop
            // --------------------------------------------------
            const double lbrhob = rhob  / FL.dRho;
            const double lbrhom = rhobm / FL.dRho;
            const double lbphb  = p_bn  / FL.dP;
            const double lbphm  = pm    / FL.dP;

            const double lbu2b = lbub[0]*lbub[0] + lbub[1]*lbub[1] + lbub[2]*lbub[2];
            const double lbu2m = lbum[0]*lbum[0] + lbum[1]*lbum[1] + lbum[2]*lbum[2];

            // --------------------------------------------------
            // 4) Update only incoming populations (normal component == n)
            // --------------------------------------------------
            //cout<<"axis= "<<axis<<endl;
            
            if (axis == 0)
            {
                const int ii = n;
                //for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                {
                    const double lbw   = FL.lbWeights[ii+1][jj+1][kk+1];

                    const double lbcub = ii*lbub[0] + jj*lbub[1] + kk*lbub[2];
                    const double lbcum = ii*lbum[0] + jj*lbum[1] + kk*lbum[2];

                    const double Feqb  = lbrhob * lbw * (1.0 + lbcub/lbcs2 - lbu2b/(2.0*lbcs2) + (lbcub*lbcub)/(2.0*lbcs2*lbcs2));
                    const double Feqm  = lbrhom * lbw * (1.0 + lbcum/lbcs2 - lbu2m/(2.0*lbcs2) + (lbcum*lbcum)/(2.0*lbcs2*lbcs2));
                    const double geqb  = Feqb*lbcs2 + (lbphb - lbrhob*lbcs2) * lbw;
                    const double geqm  = Feqm*lbcs2 + (lbphm - lbrhom*lbcs2) * lbw;

                    const double gbarm = FL.lbPopulations(im, jm, km, {nf})(ii,jj,kk);
                    FL.lbPopulations(i,j,k,{nf})(ii,jj,kk) = geqb + gbarm - geqm;
                    FL.lbPopulations(i,j,k,{nf})(ii,jj,kk) = geqb;
                }
            }
            else if (axis == 1)
            {
                const int jj = n;
                for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                {
                    const double lbw   = FL.lbWeights[ii+1][jj+1][kk+1];

                    const double lbcub = ii*lbub[0] + jj*lbub[1] + kk*lbub[2];
                    const double lbcum = ii*lbum[0] + jj*lbum[1] + kk*lbum[2];

                    const double Feqb  = lbrhob * lbw * (1.0 + lbcub/lbcs2 - lbu2b/(2.0*lbcs2) + (lbcub*lbcub)/(2.0*lbcs2*lbcs2));
                    const double Feqm  = lbrhom * lbw * (1.0 + lbcum/lbcs2 - lbu2m/(2.0*lbcs2) + (lbcum*lbcum)/(2.0*lbcs2*lbcs2));
                    const double geqb  = Feqb*lbcs2 + (lbphb - lbrhob*lbcs2) * lbw;
                    const double geqm  = Feqm*lbcs2 + (lbphm - lbrhom*lbcs2) * lbw;

                    const double gbarm = FL.lbPopulations(im, jm, km, {nf})(ii,jj,kk);
                    FL.lbPopulations(i,j,k,{nf})(ii,jj,kk) = geqb + gbarm - geqm;
                }
            }
            else // axis == 2
            {
                const int kk = n;
                for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                {
                    const double lbw   = FL.lbWeights[ii+1][jj+1][kk+1];

                    const double lbcub = ii*lbub[0] + jj*lbub[1] + kk*lbub[2];
                    const double lbcum = ii*lbum[0] + jj*lbum[1] + kk*lbum[2];

                    const double Feqb  = lbrhob * lbw * (1.0 + lbcub/lbcs2 - lbu2b/(2.0*lbcs2) + (lbcub*lbcub)/(2.0*lbcs2*lbcs2));
                    const double Feqm  = lbrhom * lbw * (1.0 + lbcum/lbcs2 - lbu2m/(2.0*lbcs2) + (lbcum*lbcum)/(2.0*lbcs2*lbcs2));
                    const double geqb  = Feqb*lbcs2 + (lbphb - lbrhob*lbcs2) * lbw;
                    const double geqm  = Feqm*lbcs2 + (lbphm - lbrhom*lbcs2) * lbw;

                    const double gbarm = FL.lbPopulations(im, jm, km, {nf})(ii,jj,kk);
                    FL.lbPopulations(i,j,k,{nf})(ii,jj,kk) = geqb + gbarm - geqm;
                }
            }
        }
    };

    // -------------------------
    // Select axis, check ownership, loop patch
    // -------------------------
    const int axis = xConst ? 0 : (yConst ? 1 : 2);
    if (!owns_face(axis)) return;

    if (axis == 0)
    {
        const int iFace = ia - Grid.OffsetX;
        if ((n > 0 && iFace != 0) || (n < 0 && iFace != Nx-1)) return;

        const int jMin = clamp(std::min(ja, jb), 0, Ny-1);
        const int jMax = clamp(std::max(ja, jb), 0, Ny-1);
        const int kMin = clamp(std::min(ka, kb), 0, Nz-1);
        const int kMax = clamp(std::max(ka, kb), 0, Nz-1);

        for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
            apply_at_node(iFace, j, k, 0);
        return;
    }

    if (axis == 1)
    {
        const int jFace = ja - Grid.OffsetY;
        if ((n > 0 && jFace != 0) || (n < 0 && jFace != Ny-1)) return;

        const int iMin = clamp(std::min(ia, ib), 0, Nx-1);
        const int iMax = clamp(std::max(ia, ib), 0, Nx-1);
        const int kMin = clamp(std::min(ka, kb), 0, Nz-1);
        const int kMax = clamp(std::max(ka, kb), 0, Nz-1);

        for (int k = kMin; k <= kMax; ++k)
        for (int i = iMin; i <= iMax; ++i)
            apply_at_node(i, jFace, k, 1);

        return;
    }

    if(axis == 2)
    {
        const int kFace = ka - Grid.OffsetZ;
        if ((n > 0 && kFace != 0) || (n < 0 && kFace != Nz-1)) return;

        const int iMin = clamp(std::min(ia, ib), 0, Nx-1);
        const int iMax = clamp(std::max(ia, ib), 0, Nx-1);
        const int jMin = clamp(std::min(ja, jb), 0, Ny-1);
        const int jMax = clamp(std::max(ja, jb), 0, Ny-1);

        for (int j = jMin; j <= jMax; ++j)
        for (int i = iMin; i <= iMax; ++i)
            apply_at_node(i, j, kFace, 2);

        return;
    }
}

double MixtureFlow::CalculateMassFlowRate(const FlowSolverLBM& FL, const dVector3& A, const dVector3& B, int nSign)
{
    if (nSign != 1 && nSign != -1) return 0.0;

    const int ia = A.getX(), ja = A.getY(), ka = A.getZ();
    const int ib = B.getX(), jb = B.getY(), kb = B.getZ();

    const bool xConst = (ia == ib);
    const bool yConst = (ja == jb);
    const bool zConst = (ka == kb);

    // Must be exactly one axis-aligned face (one constant coordinate)
    if ((int)xConst + (int)yConst + (int)zConst != 1) return 0.0;

    const int Nx = Grid.Nx, Ny = Grid.Ny, Nz = Grid.Nz;

    auto clamp = [](int v, int lo, int hi){ return std::max(lo, std::min(hi, v)); };

    // Determine axis normal to the plane
    const int axis = xConst ? 0 : (yConst ? 1 : 2);

    // Convert GLOBAL fixed coordinate to LOCAL
    const int faceLocal =
        (axis == 0) ? (ia - Grid.OffsetX) :
        (axis == 1) ? (ja - Grid.OffsetY) :
                      (ka - Grid.OffsetZ);

    // If this MPI subdomain doesn't contain the plane, local contribution is zero
    if (axis == 0) { if (faceLocal < 0 || faceLocal >= Nx) return 0.0; }
    if (axis == 1) { if (faceLocal < 0 || faceLocal >= Ny) return 0.0; }
    if (axis == 2) { if (faceLocal < 0 || faceLocal >= Nz) return 0.0; }

    // Tangential ranges (GLOBAL -> LOCAL)
    const int la_i = ia - Grid.OffsetX, lb_i = ib - Grid.OffsetX;
    const int la_j = ja - Grid.OffsetY, lb_j = jb - Grid.OffsetY;
    const int la_k = ka - Grid.OffsetZ, lb_k = kb - Grid.OffsetZ;

    int iMin=0,iMax=0,jMin=0,jMax=0,kMin=0,kMax=0;

    if (axis == 0)
    {
        jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);
        kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);
    }
    else if (axis == 1)
    {
        iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);
    }
    else // axis == 2
    {
        iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);
    }

    // Differential "area" element:
    // 3D: face area dx^2
    // 2D: line length dx (per unit depth)
    // Use dN? == 0 to infer inactive dimension(s)
    const int activeDims =
        (Grid.dNx > 0 ? 1 : 0) +
        (Grid.dNy > 0 ? 1 : 0) +
        (Grid.dNz > 0 ? 1 : 0);

    double dA = 1.0;
    if (activeDims >= 3) dA = Grid.dx * Grid.dx;
    else if (activeDims == 2) dA = Grid.dx;
    else dA = 1.0; // 1D fallback

    double mdot = 0.0;

    auto add_node_flux = [&](int i, int j, int k)
    {
        if (FL.Obstacle(i,j,k)) return;

        double massFluxNormal = 0.0; // sum_c rho_c * u_{c,n}
        for (size_t c = 0; c < FL.N_Fluid_Comp; ++c)
        {
            const double rho_c = DensityMix(i,j,k);
            const double uN_c  = VelocityMix(i,j,k)[axis];
            massFluxNormal += rho_c * uN_c;
        }
        mdot += (double)nSign * massFluxNormal * dA;
    };

    // Loop the patch on this subdomain
    if (axis == 0)
    {
        const int i = faceLocal;
        for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
            add_node_flux(i, j, k);
        return mdot;
    }

    if (axis == 1)
    {
        const int j = faceLocal;
        for (int k = kMin; k <= kMax; ++k)
        for (int i = iMin; i <= iMax; ++i)
            add_node_flux(i, j, k);
        return mdot;
    }

    // axis == 2
    {
        const int k = faceLocal;
        for (int j = jMin; j <= jMax; ++j)
        for (int i = iMin; i <= iMax; ++i)
            add_node_flux(i, j, k);
        return mdot;
    }
}

void MixtureFlow::SetVelocityBoundary(FlowSolverLBM& FL, const dVector3& UinMean, const dVector3& A, const dVector3& B, int n, bool usePoiseuille)
{
    if (n != 1 && n != -1) return;

    const int ia = A.getX(), ja = A.getY(), ka = A.getZ();
    const int ib = B.getX(), jb = B.getY(), kb = B.getZ();

    const bool xConst = (ia == ib);
    const bool yConst = (ja == jb);
    const bool zConst = (ka == kb);

    // Must be exactly one axis-aligned face
    if ((int)xConst + (int)yConst + (int)zConst != 1) return;

    const int Nx = Grid.Nx, Ny = Grid.Ny, Nz = Grid.Nz;

    auto clamp = [](int v, int lo, int hi){ return std::max(lo, std::min(hi, v)); };

    auto owns_face = [&](int axis)->bool
    {
        if (axis == 0) return ( (n > 0 && Grid.OffsetX == 0) ||
                                (n < 0 && Grid.OffsetX + Nx == Grid.TotalNx) );
        if (axis == 1) return ( (n > 0 && Grid.OffsetY == 0) ||
                                (n < 0 && Grid.OffsetY + Ny == Grid.TotalNy) );
        return            ( (n > 0 && Grid.OffsetZ == 0) ||
                            (n < 0 && Grid.OffsetZ + Nz == Grid.TotalNz) );
    };

    // Determine which face axis we are on
    const int axis = xConst ? 0 : (yConst ? 1 : 2);
    if (!owns_face(axis)) return;

    // Convert face coordinate from global to local
    const int faceLocal =
        (axis == 0) ? (ia - Grid.OffsetX) :
        (axis == 1) ? (ja - Grid.OffsetY) :
                      (ka - Grid.OffsetZ);

    // Must be the local boundary plane
    if (axis == 0)
    {
        if (faceLocal < 0 || faceLocal >= Nx) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Nx-1)) return;
    }
    else if (axis == 1)
    {
        if (faceLocal < 0 || faceLocal >= Ny) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Ny-1)) return;
    }
    else
    {
        if (faceLocal < 0 || faceLocal >= Nz) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Nz-1)) return;
    }

    // Segment AB in GLOBAL index space (used for Poiseuille coordinate)
    const double tX = double(ib - ia);
    const double tY = double(jb - ja);
    const double tZ = double(kb - ka);
    const double t2 = tX*tX + tY*tY + tZ*tZ; // length^2 in index space

    // Helper: compute imposed velocity at a node (global indices gx,gy,gz)
    auto imposed_velocity = [&](int gx, int gy, int gz)->dVector3
    {
        if (!usePoiseuille) return UinMean;
        if (t2 <= 0.0)      return UinMean; // A==B -> fallback to uniform

        // projection s of (P-A) on AB, clamped to [0,1]
        const double dxp = double(gx - ia);
        const double dyp = double(gy - ja);
        const double dzp = double(gz - ka);

        double s = (dxp*tX + dyp*tY + dzp*tZ) / t2;
        if (s < 0.0) s = 0.0;
        if (s > 1.0) s = 1.0;

        // Poiseuille factor with mean=1 and max=1.5
        const double f = 6.0 * s * (1.0 - s);

        return UinMean * f;
    };

    // -------------------------
    // Core kernel (optimized: no dead work). Takes imposed u explicitly.
    // -------------------------
    auto apply_at_node = [&](int i, int j, int k, const dVector3& u, int axisLocal)
    {
        if (FL.Obstacle(i,j,k)) return;

        int ip=i, jp=j, kp=k;
        if (axisLocal == 0) ip = i + n;
        if (axisLocal == 1) jp = j + n;
        if (axisLocal == 2) kp = k + n;

        if (ip < 0 || ip >= Nx || jp < 0 || jp >= Ny || kp < 0 || kp >= Nz) return;
        if (FL.Obstacle(ip,jp,kp)) return;

        for (size_t nf = 0; nf < FL.N_Fluid_Comp; ++nf)
        {
            const double rhop = DensityMix(ip,jp,kp);
            const double rho  = DensityMix(i ,j ,k );

            // Forward-cell reconstruction for pp and Mvp
            double pp = 0.0;
            dVector3 Mvp = {0.0, 0.0, 0.0};

            for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
            for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
            for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
            {
                const double rrp = FL.lbPopulations(ip,jp,kp,{nf})(ii,jj,kk);
                pp += rrp * FL.dP;
                Mvp[0] += rrp * ii * FL.dm;
                Mvp[1] += rrp * jj * FL.dm;
                Mvp[2] += rrp * kk * FL.dm;
            }

            // Boundary pressure uses forward pressure (your original approach)
            double p = pp;

            const dVector3 up = (Mvp * 3.0 + FL.ForceDensity(ip,jp,kp,{nf}) * FL.dt / 2.0) / rhop;

            // Keep your pressure corrections
            pp += FL.cs2/2.0 * FL.dt *
                  ( up[0]*GradDensity(ip,jp,kp)[0]
                  + up[1]*GradDensity(ip,jp,kp)[1]
                  + up[2]*GradDensity(ip,jp,kp)[2]
                  + rhop *FL.DivVel(ip,jp,kp,{nf}) );

            p  += FL.cs2/2.0 * FL.dt *
                  ( u[0]*GradDensity(i,j,k)[0]
                  + u[1]*GradDensity(i,j,k)[1]
                  + u[2]*GradDensity(i,j,k)[2]
                  + rho * FL.DivVel(i,j,k,{nf}) );

            const dVector3 lbup = up * (FL.dt / Grid.dx);
            const dVector3 lbu  = u  * (FL.dt / Grid.dx);

            const double lbu2p = lbup[0]*lbup[0] + lbup[1]*lbup[1] + lbup[2]*lbup[2];
            const double lbu2  = lbu [0]*lbu [0] + lbu [1]*lbu [1] + lbu [2]*lbu [2];

            const double lbrho  = rho  / FL.dRho;
            const double lbrhop = rhop / FL.dRho;

            const double lbp  = p  / FL.dP;
            const double lbpp = pp / FL.dP;

            // Incoming-only update: normal component == n
            if (axisLocal == 0)
            {
                const int ii = n;
                for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                {
                    const double lbw  = FL.lbWeights[ii+1][jj+1][kk+1];

                    const double lbcu  = ii*lbu[0]  + jj*lbu[1]  + kk*lbu[2];
                    const double lbcup = ii*lbup[0] + jj*lbup[1] + kk*lbup[2];

                    const double Feq  = lbrho  * lbw * (1.0 + lbcu /lbcs2 - lbu2 /(2.0*lbcs2) + (lbcu *lbcu)/(2.0*lbcs2*lbcs2));
                    const double Feqp = lbrhop * lbw * (1.0 + lbcup/lbcs2 - lbu2p/(2.0*lbcs2) + (lbcup*lbcup)/(2.0*lbcs2*lbcs2));

                    const double geq  = Feq *lbcs2  + (lbp  - lbrho *lbcs2)  * lbw;
                    const double geqp = Feqp*lbcs2  + (lbpp - lbrhop*lbcs2) * lbw;

                    const double gbarp = FL.lbPopulations(ip,jp,kp,{nf})(ii,jj,kk);
                    FL.lbPopulations(i,j,k,{nf})(ii,jj,kk) = geq + gbarp - geqp;
                }
            }
            else if (axisLocal == 1)
            {
                const int jj = n;
                for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                {
                    const double lbw  = FL.lbWeights[ii+1][jj+1][kk+1];

                    const double lbcu  = ii*lbu[0]  + jj*lbu[1]  + kk*lbu[2];
                    const double lbcup = ii*lbup[0] + jj*lbup[1] + kk*lbup[2];

                    const double Feq  = lbrho  * lbw * (1.0 + lbcu /lbcs2 - lbu2 /(2.0*lbcs2) + (lbcu *lbcu)/(2.0*lbcs2*lbcs2));
                    const double Feqp = lbrhop * lbw * (1.0 + lbcup/lbcs2 - lbu2p/(2.0*lbcs2) + (lbcup*lbcup)/(2.0*lbcs2*lbcs2));

                    const double geq  = Feq *lbcs2  + (lbp  - lbrho *lbcs2)  * lbw;
                    const double geqp = Feqp*lbcs2  + (lbpp - lbrhop*lbcs2) * lbw;

                    const double gbarp = FL.lbPopulations(ip,jp,kp,{nf})(ii,jj,kk);
                    FL.lbPopulations(i,j,k,{nf})(ii,jj,kk) = geq + gbarp - geqp;
                }
            }
            else // axisLocal == 2
            {
                const int kk = n;
                for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                {
                    const double lbw  = FL.lbWeights[ii+1][jj+1][kk+1];

                    const double lbcu  = ii*lbu[0]  + jj*lbu[1]  + kk*lbu[2];
                    const double lbcup = ii*lbup[0] + jj*lbup[1] + kk*lbup[2];

                    const double Feq  = lbrho  * lbw * (1.0 + lbcu /lbcs2 - lbu2 /(2.0*lbcs2) + (lbcu *lbcu)/(2.0*lbcs2*lbcs2));
                    const double Feqp = lbrhop * lbw * (1.0 + lbcup/lbcs2 - lbu2p/(2.0*lbcs2) + (lbcup*lbcup)/(2.0*lbcs2*lbcs2));

                    const double geq  = Feq *lbcs2  + (lbp  - lbrho *lbcs2)  * lbw;
                    const double geqp = Feqp*lbcs2  + (lbpp - lbrhop*lbcs2) * lbw;

                    const double gbarp = FL.lbPopulations(ip,jp,kp,{nf})(ii,jj,kk);
                    FL.lbPopulations(i,j,k,{nf})(ii,jj,kk) = geq + gbarp - geqp;
                }
            }
        }
    };

    // -------------------------
    // Patch bounds: convert global A,B to local for tangential indices
    // -------------------------
    const int la_i = ia - Grid.OffsetX;
    const int la_j = ja - Grid.OffsetY;
    const int la_k = ka - Grid.OffsetZ;

    const int lb_i = ib - Grid.OffsetX;
    const int lb_j = jb - Grid.OffsetY;
    const int lb_k = kb - Grid.OffsetZ;

    if (axis == 0)
    {
        const int jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        const int jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);
        const int kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        const int kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);

        for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
        {
            const int gx = faceLocal + Grid.OffsetX; // should equal ia
            const int gy = j + Grid.OffsetY;
            const int gz = k + Grid.OffsetZ;

            const dVector3 u = imposed_velocity(gx, gy, gz);
            apply_at_node(faceLocal, j, k, u, 0);
        }
        return;
    }

    if (axis == 1)
    {
        const int iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        const int iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        const int kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        const int kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);

        for (int k = kMin; k <= kMax; ++k)
        for (int i = iMin; i <= iMax; ++i)
        {
            const int gx = i + Grid.OffsetX;
            const int gy = faceLocal + Grid.OffsetY; // should equal ja
            const int gz = k + Grid.OffsetZ;

            const dVector3 u = imposed_velocity(gx, gy, gz);
            apply_at_node(i, faceLocal, k, u, 1);
        }
        return;
    }

    // axis == 2
    {
        const int iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        const int iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        const int jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        const int jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);

        for (int j = jMin; j <= jMax; ++j)
        for (int i = iMin; i <= iMax; ++i)
        {
            const int gx = i + Grid.OffsetX;
            const int gy = j + Grid.OffsetY;
            const int gz = faceLocal + Grid.OffsetZ; // should equal ka

            const dVector3 u = imposed_velocity(gx, gy, gz);
            apply_at_node(i, j, faceLocal, u, 2);
        }
        return;
    }
}

void MixtureFlow::SetOpenBoundary(FlowSolverLBM& FL, ExtrapolateOrder order, const dVector3& A, const dVector3& B,  int n)
{
    if (n != 1 && n != -1) return;

    const int ia = A.getX(), ja = A.getY(), ka = A.getZ();
    const int ib = B.getX(), jb = B.getY(), kb = B.getZ();

    const bool xConst = (ia == ib);
    const bool yConst = (ja == jb);
    const bool zConst = (ka == kb);

    // Must be exactly one axis-aligned face
    if ((int)xConst + (int)yConst + (int)zConst != 1) return;

    const int Nx = Grid.Nx, Ny = Grid.Ny, Nz = Grid.Nz;

    auto clamp = [](int v, int lo, int hi){ return std::max(lo, std::min(hi, v)); };

    auto owns_face = [&](int axis)->bool
    {
        if (axis == 0) return ( (n > 0 && Grid.OffsetX == 0) ||
                                (n < 0 && Grid.OffsetX + Nx == Grid.TotalNx) );
        if (axis == 1) return ( (n > 0 && Grid.OffsetY == 0) ||
                                (n < 0 && Grid.OffsetY + Ny == Grid.TotalNy) );
        return            ( (n > 0 && Grid.OffsetZ == 0) ||
                            (n < 0 && Grid.OffsetZ + Nz == Grid.TotalNz) );
    };

    const int axis = xConst ? 0 : (yConst ? 1 : 2);
    if (!owns_face(axis)) return;

    const int faceLocal =
        (axis == 0) ? (ia - Grid.OffsetX) :
        (axis == 1) ? (ja - Grid.OffsetY) :
                      (ka - Grid.OffsetZ);

    // Must be on local boundary plane
    if (axis == 0)
    {
        if (faceLocal < 0 || faceLocal >= Nx) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Nx-1)) return;
    }
    else if (axis == 1)
    {
        if (faceLocal < 0 || faceLocal >= Ny) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Ny-1)) return;
    }
    else
    {
        if (faceLocal < 0 || faceLocal >= Nz) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Nz-1)) return;
    }

    // Patch bounds in LOCAL coordinates
    const int la_i = ia - Grid.OffsetX, lb_i = ib - Grid.OffsetX;
    const int la_j = ja - Grid.OffsetY, lb_j = jb - Grid.OffsetY;
    const int la_k = ka - Grid.OffsetZ, lb_k = kb - Grid.OffsetZ;

    int iMin=0,iMax=0,jMin=0,jMax=0,kMin=0,kMax=0;

    if (axis == 0)
    {
        jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);
        kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);
    }
    else if (axis == 1)
    {
        iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);
    }
    else
    {
        iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);
    }

    const bool firstOrder = (order == ExtrapolateOrder::FirstOrder);

    // Core: boundary cell b gets values from interior along normal
    auto apply_at_node = [&](int bi, int bj, int bk)
    {
        if (FL.Obstacle(bi,bj,bk)) return;

        int im=bi, jm=bj, km=bk;
        int imm=bi, jmm=bj, kmm=bk;

        if (axis == 0) { im = bi + n;  imm = bi + 2*n; }
        if (axis == 1) { jm = bj + n;  jmm = bj + 2*n; }
        if (axis == 2) { km = bk + n;  kmm = bk + 2*n; }

        // Need at least first interior
        if (im < 0 || im >= Nx || jm < 0 || jm >= Ny || km < 0 || km >= Nz) return;
        if (FL.Obstacle(im,jm,km)) return;

        // For first order, need second interior too
        if (firstOrder)
        {
            if (imm < 0 || imm >= Nx || jmm < 0 || jmm >= Ny || kmm < 0 || kmm >= Nz) return;
            if (FL.Obstacle(imm,jmm,kmm)) return;
        }

        for (size_t c = 0; c < FL.N_Fluid_Comp; ++c)
        {
            // Update ONLY incoming directions (normal component == n)
            if (axis == 0)
            {
                const int ii = n;
                for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                {
                    const double f1 = FL.lbPopulations(im ,jm ,km ,{c})(ii,jj,kk);
                    if (!firstOrder)
                    {
                        FL.lbPopulations(bi,bj,bk,{c})(ii,jj,kk) = f1;
                    }
                    else
                    {
                        const double f2 = FL.lbPopulations(imm,jmm,kmm,{c})(ii,jj,kk);
                        FL.lbPopulations(bi,bj,bk,{c})(ii,jj,kk) = 2.0*f1 - f2;
                    }
                }
            }
            else if (axis == 1)
            {
                const int jj = n;
                for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                for (int kk = -Grid.dNz; kk <= Grid.dNz; ++kk)
                {
                    const double f1 = FL.lbPopulations(im ,jm ,km ,{c})(ii,jj,kk);
                    if (!firstOrder)
                    {
                        FL.lbPopulations(bi,bj,bk,{c})(ii,jj,kk) = f1;
                    }
                    else
                    {
                        const double f2 = FL.lbPopulations(imm,jmm,kmm,{c})(ii,jj,kk);
                        FL.lbPopulations(bi,bj,bk,{c})(ii,jj,kk) = 2.0*f1 - f2;
                    }
                }
            }
            else // axis == 2
            {
                const int kk = n;
                for (int ii = -Grid.dNx; ii <= Grid.dNx; ++ii)
                for (int jj = -Grid.dNy; jj <= Grid.dNy; ++jj)
                {
                    const double f1 = FL.lbPopulations(im ,jm ,km ,{c})(ii,jj,kk);
                    if (!firstOrder)
                    {
                        FL.lbPopulations(bi,bj,bk,{c})(ii,jj,kk) = f1;
                    }
                    else
                    {
                        const double f2 = FL.lbPopulations(imm,jmm,kmm,{c})(ii,jj,kk);
                        FL.lbPopulations(bi,bj,bk,{c})(ii,jj,kk) = 2.0*f1 - f2;
                    }
                }
            }
        }
    };

    // Loop boundary patch
    if (axis == 0)
    {
        for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
            apply_at_node(faceLocal, j, k);
        return;
    }

    if (axis == 1)
    {
        for (int k = kMin; k <= kMax; ++k)
        for (int i = iMin; i <= iMax; ++i)
            apply_at_node(i, faceLocal, k);
        return;
    }

    // axis == 2
    for (int j = jMin; j <= jMax; ++j)
    for (int i = iMin; i <= iMax; ++i)
        apply_at_node(i, j, faceLocal);
}
