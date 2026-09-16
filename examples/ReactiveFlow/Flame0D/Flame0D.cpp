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

#include "ReactiveFlows/Species.h"
#include "ReactiveFlows/Energy.h"
#include "ReactiveFlows/Transport.h"
#include "ReactiveFlows/MixtureFlow.h"
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

    Settings                                OPSettings(SettingsInputFile);
    RunTimeControl                          RTC(OPSettings,InputFile);
    PhaseField                              Phase(OPSettings);
    BoundaryConditions                      BC(OPSettings,InputFile);
    FlowSolverLBM                           FL(OPSettings, InputFile);
    FL.SetTimeStep(RTC.dt);
    Velocities                              Vel(OPSettings);
    TimeInfo                                Timer;
    MixtureFlow                             MF(OPSettings, InputFile);
    Energy                                  EN(OPSettings, InputFile);
    Species                                 SP(OPSettings, InputFile);
    SolidBody                               SB(OPSettings, InputFile);
    Transport                               TR(OPSettings, InputFile);
    ThermoChemistry                         TC(InputFile);

    std::shared_ptr<Cantera::Solution>      sol    = Cantera::newSolution(TC.ReactionMechanism, TC.GasPhaseName, TC.GasTransportName);
    std::shared_ptr<Cantera::ThermoPhase>   gas    = sol-> thermo();
    std::shared_ptr<Cantera::Transport>     transp = sol-> transport();
    std::shared_ptr<Cantera::Kinetics>      kin    = sol-> kinetics();

    TC.InitializeGasPhase(*gas, EN.T0, FL.Pth0, SettingsInputFile);

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
        size_t idx0 = Initializations::Single(Phase, 0, BC);
        if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
        {
            cout<<"PhaseField number ["<<idx0<<"]= is created"<<endl;
        }

        MF.DetectObstacles(FL, Phase, BC, MF.DI);

        TC.ExportingGasSpeciesData(SP);
        SP.SetInitial(Phase, FL, EN,
                       TC.BurntGasMixture, TC.GasMixture, TC.AirMixture,
                       TC.BurntGasTemp,    EN.T0,          EN.T0);
        SP.SetBoundaryConditions(BC);
        EN.SetBoundaryConditions(BC);
        SP.UpdateFields();
        EN.UpdateFields();
        TC.CalculateGasKinetics(EN, SP, MF, FL, *gas, *transp, *kin, BC);
        SP.SetBoundaryConditions(BC);
        EN.SetBoundaryConditions(BC);

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
            SP.WriteVTK(OPSettings, 0);
            EN.WriteVTKHHR(OPSettings, 0);
            MF.WriteVTK(OPSettings,FL, 0);
        }
        if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
        {
            std::cout << "Done\n";
        }
    }

    if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
    {
        SB.writeData(OPSettings.TextDir,"TimeStep (s)", RTC.dt);
        SB.writeData(OPSettings.TextDir,"Initial temperature (K) = ", EN.T0);
        SB.writeData(OPSettings.TextDir,"Initial pressure (Pa) = ", FL.Pth0);
        for (size_t i = 0; i < TC.GasSpeciesNames.size(); i++)
        {
            SB.writeData(OPSettings.TextDir,"Species " + std::to_string(i), TC.GasSpeciesNames({i}));
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
        SP.UpdateFields();
        EN.UpdateFields();
        TR.CalculateAdvectionDiffusion(Phase, SP, EN, FL, MF, SB, RTC.dt);
        TR.CalculateReaction(SP, EN, MF, FL, OPSettings.Grid.Active(), TC.FuelConsumptionRate, TC.GasFuelIndex, RTC.dt);
        TC.CalculateGasKinetics(EN, SP, MF, FL, *gas, *transp, *kin, BC);
        SP.SetBoundaryConditions(BC);
        EN.SetBoundaryConditions(BC);
        MF.CalculateDensityGradient(FL, EN);
        FL.CalculateHydrodynamicPressureAndMomentum(Vel);
    	FL.CalculateFluidVelocities(Vel, Phase, BC);
        MF.CalculateVelocityAndPressure(FL);
		MF.ApplyForces(Phase, FL);
        MF.CalculateDivergenceVelocity(Phase, FL, EN, SP, SB, RTC.dt);

        if (RTC.WriteVTK())
        {
            EN.WriteVTKTemperature(OPSettings,RTC.tStep);
            SP.WriteVTK(OPSettings, RTC.tStep);
            EN.WriteVTKHHR(OPSettings, RTC.tStep);
        }

        if(RTC.WriteToScreen())
        {
            if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
            {
                cout <<"Temperature (K) : "<<EN.T(0,0,0)<<endl;
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
