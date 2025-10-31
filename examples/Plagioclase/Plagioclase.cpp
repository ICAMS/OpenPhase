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

 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Julia Kindin; Hesham Salama
 *
 */

#include "Settings.h"
#include "InterfaceProperties.h"
#include "RunTimeControl.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "BoundaryConditions.h" 
#include "Initializations.h"
#include "Tools/TimeInfo.h"

using namespace std;
using namespace openphase;
/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{ 
#ifdef MPI_PARALLEL
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif
    Settings                                OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                          RTC(OPSettings);
    PhaseField                              Phi(OPSettings);
    DoubleObstacle                          DO(OPSettings);
    InterfaceProperties                     IP(OPSettings);
    EquilibriumPartitionDiffusionBinary     DF(OPSettings);
    Composition                             Cx(OPSettings);
    Temperature                             Tx(OPSettings);
    DrivingForce                            dG(OPSettings);
    BoundaryConditions                      BC(OPSettings);
    TimeInfo                                Timer(OPSettings, "Execution Time Statistics");

    if(RTC.Restart)
    {
        cout << "Restart data being read!" << endl;
          
        Phi.Read(OPSettings, BC, RTC.tStart);
        Cx.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC, RTC.tStart);
        cout << "Done reading restart parameters!" << endl;
    }
    else 
    {
        int idx0 = Initializations::Single(Phi, 0, BC);
        //int idx1 = Initializations::Sphere(Phi, 1, OPSettings.iWidth, OPSettings.Dimensions.Nx/2, OPSettings.Dimensions.Ny/2, OPSettings.Dimensions.Nz/2, BC);
        int idx1 = Phi.PlantGrainNucleus(1, OPSettings.Grid.Nx/2, OPSettings.Grid.Ny/2, OPSettings.Grid.Nz/2);
        // Orienting the nucleus as desired
        EulerAngles locAngle({0.0, Pi/6.0, 0.0}, XYZ);
        Phi.FieldsProperties[idx1].Orientation = locAngle.getQuaternion();
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
    }

    cout << "Initialization stage done!" << endl;
    cout << "Entering the Time Loop!!!" << endl;

    //-------------- The Time Loop -------------//

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        dG.Clear();
        IP.Set(Phi, BC);

        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);

        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
        dG.Average(Phi, BC);
        dG.MergePhaseFieldIncrements(Phi, IP);

        Timer.SetTimeStamp("Get PsiDot");

        Timer.SetTimeStamp("Driving Force merge to Psi");
        Phi.NormalizeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Psi Normalize");

        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
        Timer.SetTimeStamp("Solve diffusion");
        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
        
        Timer.SetTimeStamp("Set temperature");
        Phi.MergeIncrements(BC, RTC.dt);
        
        Timer.SetTimeStamp("Merge Phase Fields");

        if (RTC.WriteVTK())
        {            
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteVTK(OPSettings, RTC.tStep);
            Tx.WriteVTK(OPSettings, RTC.tStep);
            //Cx.WriteStatistics(RTC.tStep, RTC.dt);
            IP.WriteVTK(OPSettings, RTC.tStep);
        }
        if (RTC.WriteRawData())
        {    
            Phi.Write(OPSettings, RTC.tStep);
            Cx.Write(OPSettings, RTC.tStep);
            Tx.Write(OPSettings, RTC.tStep);
        }
        Timer.SetTimeStamp("File Output");
        
        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", I_En);
            ConsoleOutput::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            //Timer.PrintWallClockSummary();
        }
    } //end time loop
   
    return 0;
#ifdef MPI_PARALLEL
    }
    MPI_Finalize ();
#endif
}
