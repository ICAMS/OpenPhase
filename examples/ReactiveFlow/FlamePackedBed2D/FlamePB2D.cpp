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

 *   File created :   2023
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung; Reza Namdar 
 *
 */

#include "AdvectionHR.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Temperature.h"
#include "TextOutput.h"
#include "Tools/TimeInfo.h"
#include "Velocities.h"
#include "HeatSources.h"

#include "ReactiveFlows/SpeciesTransport.h"
#include "ReactiveFlows/EnergyTransport.h"
#include "ReactiveFlows/FlowMixture.h"
#include "ReactiveFlows/SolidBody.h"
#include "ReactiveFlows/ThermoChemistry.h"

#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include "cantera/base/Array.h"
#include "cantera/base/Solution.h"
#include <array>

using namespace openphase;
using namespace std;

int main(int argc, char *argv[])
{ 

#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_FUNNELED, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
#endif

    //feenableexcept(FE_DIVBYZERO);  // Devision by zero
    //feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);   // the result was too large
    //feenableexcept(FE_UNDERFLOW);  // the result was too small

    std::string InputFile            = DefaultInputFileName;
    std::string SettingsInputFile    = "SettingsInput.opi";
    if (argc > 1) InputFile = argv[1]; 

    ThermoChemistry                         TC(InputFile);
    std::shared_ptr<Cantera::Solution>      sol    = Cantera::newSolution(TC.ReactionMechanism, TC.PhaseName, TC.TransportName);
    std::shared_ptr<Cantera::ThermoPhase>   gas    = sol-> thermo();
    std::shared_ptr<Cantera::Transport>     transp = sol-> transport();
    std::shared_ptr<Cantera::Kinetics>      kin    = sol-> kinetics();

    TC.Initialize(gas,SettingsInputFile);

    Settings                                OPSettings(SettingsInputFile);
    RunTimeControl                          RTC(OPSettings,InputFile);
    PhaseField                              Phase(OPSettings);
    BoundaryConditions                      BC(OPSettings,InputFile);
    FlowSolverLBM                           FL(OPSettings, RTC.dt, InputFile);
    Velocities                              Vel(OPSettings); 
    TimeInfo                                Timer;
    FlowMixture                             FM(InputFile);
    EnergyTransport                         ET(OPSettings,InputFile);
    SpeciesTransport                        ST(OPSettings,InputFile);
    SolidBody                               SB(InputFile);
    
    if(RTC.Restart)
    {
        std::cout << "Restart data being read!\n";
        Phase.Read(OPSettings,BC, RTC.tStart);
        FL .Read(OPSettings, BC, RTC.tStart);
        RTC.tStart += 1;
        std::cout << "Done\n";
    }
    else 
    {
        double dx = FL.Grid.dx;
        size_t idx0 = Initializations::Single(Phase, 0, BC);
        if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
        {
            cout<<"PhaseField number ["<<idx0<<"]= is created"<<endl;
        }
        ET.IF_ConstTemp.push_back(false);
        ET.SurfaceTemp.push_back(ET.T0);
        ET.IF_ConstFlux.push_back(false);
        ET.SurfaceFlux.push_back(ET.HeatFlux);  

        SB.DistributeRandomSolidBodies(Phase,OPSettings,""); 
        for (size_t i = 0; i < SB.nParticles; ++i) 
        {
            size_t idx1 = Initializations::Sphere(Phase, 1, SB.rand_Circles[i].r/dx,
                                                SB.rand_Circles[i].x/dx+SB.X0DistZone, 0.0, SB.rand_Circles[i].z/dx+SB.Z0DistZone, BC);
            if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
            {
                cout<<"PhaseField number ["<<idx1<<"]= is created"<<endl;
            }
            cout<<"Particle number["<<i<<"]= is placed at position("<<SB.rand_Circles[i].x/dx+SB.X0DistZone<<","<<SB.rand_Circles[i].z/dx+SB.Z0DistZone<<")"<<endl;
            ET.IF_ConstFlux.push_back(false);
            ET.SurfaceFlux.push_back(ET.HeatFlux);
            ET.IF_ConstTemp.push_back(false);
            ET.SurfaceTemp.push_back(ET.TempSolid);
        }

        FM.DetectObstacles(FL, Phase, FM.DI);
        BC.SetX(FL.Obstacle);
        BC.SetY(FL.Obstacle);
        BC.SetZ(FL.Obstacle);

        TC.ExportingData(ST);
        
        //ET.SetInitial(Phase, FL);
        ST.SetInitial(Phase, FL, ET, 
                       TC.BurntMixture, TC.FuelMixture, TC.AirMixture,
                       TC.BurntTemp,    ET.T0,          ET.T0);
        ST.SetBoundaryConditions(ET, BC);
        ST.UpdateFields(ET);
        ST.CalculateSpeciesHeatCapacitiesAndEnthalpies(ET, FL);
        TC.UpdateProperties(ET,ST,FL, gas, transp, kin);
        ST.SetBoundaryConditions(ET, BC);

        FM.Initialize(OPSettings, Phase, FL, Vel, BC);
        FM.LBMLimits(Phase, FL, FM.MaxU, FM.lbnu, RTC.dt);

        if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
        {
            std::cout << "Write initial state to file\n";
        }
        //  Output to file
        {
            Phase.WriteVTK(OPSettings, 0);
            ET.WriteVTKTemperature(OPSettings,FL,0);
            FL.WriteVTK(OPSettings,Phase, 0);
            ST.WriteVTKMassFractions(OPSettings, 0);
            ST.WriteVTKHHR(OPSettings, 0);
            FM.WriteVTKMixtureVelocity(OPSettings,0);
        }
        if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
        {
            std::cout << "Done\n";
        }
    }
    
    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
    {   
        SB.writeData(OPSettings.TextDir,"Grid size (m)", FL.Grid.dx);
        SB.writeData(OPSettings.TextDir,"TimeStep (s)", RTC.dt);
        SB.writeData(OPSettings.TextDir,"Nx", FL.Grid.TotalNx);
        SB.writeData(OPSettings.TextDir,"Ny", FL.Grid.TotalNy);
        SB.writeData(OPSettings.TextDir,"Nz", FL.Grid.TotalNz);
        SB.writeData(OPSettings.TextDir,"dimension of simulation",OPSettings.Grid.Active());
        SB.writeData(OPSettings.TextDir,"LB kinematic viscosity ", FM.InletViscosity * RTC.dt/FL.Grid.dx/FL.Grid.dx );
        SB.writeData(OPSettings.TextDir,"If second order BB is applied", FM.SecOrdBB);
        SB.writeData(OPSettings.TextDir,"Initial temperature (K) = ",ET.T0);
        SB.writeData(OPSettings.TextDir,"Hot wall temperature (K) = ",ET.TempSolid);
        SB.writeData(OPSettings.TextDir,"Reynolds number",FL.U0X*SB.PartDiameter/FM.InletViscosity);
        SB.writeData(OPSettings.TextDir,"Inlet velocity (m/s)",FL.U0X);
        
        for (size_t iPF = 0; iPF < Phase.FieldsProperties.size(); iPF++)
        {
            if(ET.IF_ConstFlux[iPF]==true)
            {
                string dcheck = "If PF ["+std::to_string(iPF)+"] is constant flux";
                string dHF = "Heat Flux of PF ["+std::to_string(iPF)+"] is equal";
                SB.writeData(OPSettings.TextDir,dcheck, ET.IF_ConstFlux[iPF]);
                SB.writeData(OPSettings.TextDir,dHF, ET.SurfaceFlux[iPF]);
            }
        }
        for (size_t iPF = 0; iPF < Phase.FieldsProperties.size(); iPF++)
        {
            if(ET.IF_ConstTemp[iPF]==true)
            {
                string dcheck = "If PF ["+std::to_string(iPF)+"] is temperature constant";
                string dT = "Temperature of PF ["+std::to_string(iPF)+"] is equal";
                SB.writeData(OPSettings.TextDir,dcheck, ET.IF_ConstTemp[iPF]);
                SB.writeData(OPSettings.TextDir,dT, ET.SurfaceTemp[iPF]);
            }
        }
        for (size_t i = 0; i < TC.SpeciesNames.size(); i++)
        {
            SB.writeData(OPSettings.TextDir,"Species " + std::to_string(i), TC.SpeciesNames.at(i));
        }
        SB.writeData(OPSettings.TextDir,"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", "");
    }

    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
    {
        std::cout << "Entering the Time Loop!!!\n";
    }

    for(RTC.tStep = 1; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep()) 
    {
        FM.Collision(FL,Vel,BC);
        FM.Propagation(FL,Phase,BC);  

        if(FM.Updating_Velocity) FM.UpdatingInletVelocity(FL, ST.FuelConsumptionRate, TC.FuelMixture.at(TC.FuelIndex)); 
        FM.SetInletVelocity(FL,BC,Phase,FL.U0X);
        FL.SetPressureOutlet(Vel, Phase,FM.Kp);
        //FL.SetOutFlow(0);

        ST.UpdateGhostPoints(Phase, ET, FL, SB);

        ST.UpdateFields(ET);

        ST.CalculateAdvection(ET, FL, FM, RTC.dt);
        ST.CalculateDiffusion(Phase, ET, FL, RTC.dt); 
        ST.CalculateReaction(ET, FL, OPSettings.Grid.Active(), RTC.dt);
       
        if(ET.Conjugate) ET.CalculateSolidDiffusion(Phase, FL, RTC.dt);
        ST.SetFreeBCNX(Phase, ET);
        
        ST.CalculateSpeciesHeatCapacitiesAndEnthalpies(ET, FL);
        TC.UpdateProperties(ET,ST,FL, gas, transp, kin);
        ST.SetBoundaryConditions(ET, BC);
        FL.CalculateDensityGradient(Phase,Vel);
        FL.CalculateHydrodynamicPressureAndMomentum(Vel);
    	FL.CalculateFluidVelocities(Vel, Phase, BC);
        FM.UpdateMixtureVelocity(Vel);
		FL.ApplyForces(Phase, Vel);
		if(FM.DI) FM.CalculateForceDragbyPF(Phase,FL,Vel, FM.Austar);
        FM.CalculateDivergenceVelocity(Phase, FL, ET, ST, BC, RTC.dt);
        
        double InletMassFlow=0.0;
        if(Phase.Grid.OffsetX==0)
        {
            FM.CalculatingMassFlowRate(FL, Phase, InletMassFlow, 0);
        }
        #ifdef MPI_PARALLEL
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &InletMassFlow, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        double OutletMassFlow=0.0;
        if(Phase.Grid.OffsetX+Phase.Grid.Nx==Phase.Grid.TotalNx)
        {
            FM.CalculatingMassFlowRate(FL, Phase,OutletMassFlow, Phase.Grid.Nx-2);
        }
        #ifdef MPI_PARALLEL
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &OutletMassFlow, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        if (RTC.WriteVTK())
        {
            if(!ET.Conjugate) ET.SetSolidPhaseTemp(Phase,FM.DI);
            ET.WriteVTKTemperature(OPSettings,FL,RTC.tStep);
            FL.WriteVTK(OPSettings,Phase, RTC.tStep);
            ST.WriteVTKMassFractions(OPSettings, RTC.tStep);
            ST.WriteVTKHHR(OPSettings, RTC.tStep);
            FM.WriteVTKMixtureVelocity(OPSettings,RTC.tStep);
        }

        if(RTC.WriteToScreen())
        {
            if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
            {
                cout <<"Mass Flow rate at inlet  is  = "<<InletMassFlow<<" (kg/s)"<<endl; 
                cout <<"Mass Flow rate at outlet  is = "<<OutletMassFlow<<" (kg/s)"<<endl; 
                cout <<"Inlet Velocity= "<<FL.U0X<<" (m/s)"<<endl;
                cout <<"Timestep : "<<RTC.tStep << endl;
                cout <<"Physical time (s): "<<RTC.tStep*RTC.dt <<endl;
                cout <<"====================================================================" << endl;
            }
        }

    } //end time loop

    #ifdef MPI_PARALLEL
        OP_MPI_Finalize();
    #endif

    return 0;
}
