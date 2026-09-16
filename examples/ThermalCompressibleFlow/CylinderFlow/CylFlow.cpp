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

#include "ReactiveFlows/Energy.h"
#include "ReactiveFlows/Transport.h"
#include "ReactiveFlows/MixtureFlow.h"
#include "ReactiveFlows/SolidBody.h"


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
    std::string InputFile = DefaultInputFileName;
    if (argc > 1) InputFile = argv[1];

    Settings                                OPSettings(InputFile);
    RunTimeControl                          RTC(OPSettings,InputFile);
    PhaseField                              Phase(OPSettings);
    BoundaryConditions                      BC(OPSettings,InputFile);
    FlowSolverLBM                           FL(OPSettings, InputFile);
    FL.SetTimeStep(RTC.dt);
    Velocities                              Vel(OPSettings);
    TimeInfo                                Timer;
    MixtureFlow                             MF(OPSettings, InputFile);
    Energy                                  EN(OPSettings, InputFile);
    SolidBody                               SB(OPSettings, InputFile);
    Transport                               TR(OPSettings, InputFile);

    // Full-face index-space points used by the point-based boundary API below.
    const dVector3 InletA{0.0, 0.0, 0.0};
    const dVector3 InletB{0.0, double(FL.Grid.TotalNy-1), double(FL.Grid.TotalNz-1)};
    const dVector3 OutletA{double(FL.Grid.TotalNx-1), 0.0, 0.0};
    const dVector3 OutletB{double(FL.Grid.TotalNx-1), double(FL.Grid.TotalNy-1), double(FL.Grid.TotalNz-1)};

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

        SB.DistributeSolidBodies(1, 0.0, "cylindersData.txt");
        for (size_t i = 0; i < SB.nParticles; ++i)
        {
            size_t idx1 = Initializations::Sphere(Phase, 1, SB.rand_Bodies[i].r/dx,
                                                SB.rand_Bodies[i].x/dx+SB.X0DistZone, 0.0, SB.rand_Bodies[i].z/dx+SB.Z0DistZone, BC);
            EN.ThermalSurfaceCondition.push_back({ThermalCondition::Type::Conjugate, 0.0, idx1, 0.0});
            if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
            {
                cout<<"PhaseField number ["<<idx1<<"]= is created"<<endl;
            }
            if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
            {
                cout<<"Particle number["<<i<<"]= is placed at position("<<SB.rand_Bodies[i].x/dx+SB.X0DistZone<<","<<SB.rand_Bodies[i].z/dx+SB.Z0DistZone<<")"<<endl;
            }
        }

        MF.DetectObstacles(FL, Phase, BC, MF.DI);
        SB.CalculateDistanceField(Phase);
        SB.SetBoundaryConditions(BC);

        EN.SetInitial(Phase, FL, MF, SB);
        EN.SetSolidPhaseTemp(Phase,MF.DI);
        EN.SetBoundaryConditions(BC);
        EN.CalculateProperties(Phase, FL, SB, MF, BC);
        MF.SetInitial(Phase, FL, BC);

        if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
        {
            std::cout << "Write initial state to file\n";
        }
        //  Output to file
        {
            Phase.WriteVTK(OPSettings, 0);
            EN.WriteVTKTemperature(OPSettings,0);
            FL.WriteVTK(OPSettings,Phase, 0);
            MF.WriteVTK(OPSettings,FL, 0);
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
        SB.writeData(OPSettings.TextDir,"LB kinematic viscosity ", MF.InletViscosity * RTC.dt/FL.Grid.dx/FL.Grid.dx );
        SB.writeData(OPSettings.TextDir,"If second order BB is applied", MF.SecOrdBB);
        SB.writeData(OPSettings.TextDir,"Initial temperature (K) = ",EN.T0);
        SB.writeData(OPSettings.TextDir,"Hot wall temperature (K) = ",EN.TempSolid);
        SB.writeData(OPSettings.TextDir,"Reynolds number",FL.U0X*SB.PartDiameter/MF.InletViscosity);
        SB.writeData(OPSettings.TextDir,"Inlet velocity (m/s)",FL.U0X);

        for (size_t is = 0; is < EN.ThermalSurfaceCondition.size(); is++)
        {
            const auto& cond = EN.ThermalSurfaceCondition[is];
            string idxStr = std::to_string(cond.solidIdx);
            if(cond.type==ThermalCondition::Type::ConstantFlux)
            {
                SB.writeData(OPSettings.TextDir,"If PF ["+idxStr+"] is constant flux", true);
                SB.writeData(OPSettings.TextDir,"Heat Flux of PF ["+idxStr+"] is equal", cond.value);
            }
            else if(cond.type==ThermalCondition::Type::ConstantTemp)
            {
                SB.writeData(OPSettings.TextDir,"If PF ["+idxStr+"] is temperature constant", true);
                SB.writeData(OPSettings.TextDir,"Temperature of PF ["+idxStr+"] is equal", cond.value);
            }
        }
        SB.writeData(OPSettings.TextDir,"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", "");
    }

    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
    {
        std::cout << "Entering the Time Loop!!!\n";
    }

    for(RTC.tStep = 1; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        MF.Collision(FL,BC);
        MF.Propagation(FL,Phase,SB,BC);

        MF.SetVelocityBoundary(FL, dVector3{FL.U0X,0.0,0.0}, InletA, InletB, 1, MF.POISEFLOW);
        MF.SetPressureBoundary(FL, FL.Poutlet, PressureOutletMode::NonReflecting, OutletA, OutletB, -1);

        EN.UpdateFields();

        TR.CalculateAdvectionDiffusion(Phase, EN, FL, MF, SB, RTC.dt);

        if(EN.CONJUGATE) TR.CalculateSolidDiffusion(Phase, EN, FL, SB, RTC.dt);
        EN.SetOpenBoundary(OutletA, OutletB, -1);
        EN.SetBoundaryConditions(BC);

        EN.CalculateProperties(Phase, FL, SB, MF, BC);

        MF.CalculateDensityGradient(FL, EN);
        FL.CalculateHydrodynamicPressureAndMomentum(Vel);
    	FL.CalculateFluidVelocities(Vel, Phase, BC);
        MF.CalculateVelocityAndPressure(FL);
		MF.ApplyForces(Phase, FL);
        MF.CalculateDivergenceVelocity(Phase, FL, EN, SB, RTC.dt);

        double InletMassFlow=0.0;
        if(Phase.Grid.OffsetX==0)
        {
            InletMassFlow = MF.CalculateMassFlowRate(FL, InletA, InletB, 1);
        }
        #ifdef MPI_PARALLEL
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &InletMassFlow, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        double OutletMassFlow=0.0;
        if(Phase.Grid.OffsetX+Phase.Grid.Nx==Phase.Grid.TotalNx)
        {
            OutletMassFlow = MF.CalculateMassFlowRate(FL, OutletA, OutletB, -1);
        }
        #ifdef MPI_PARALLEL
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &OutletMassFlow, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        if (RTC.WriteVTK())
        {
            EN.WriteVTKTemperature(OPSettings,RTC.tStep);
            FL.WriteVTK(OPSettings,Phase, RTC.tStep);
            MF.WriteVTK(OPSettings,FL, RTC.tStep);
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
