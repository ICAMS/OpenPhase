#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
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
    Settings                        OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                      RTC(OPSettings);
    PhaseField                          Phi(OPSettings);
    DoubleObstacle                      DO(OPSettings);
    InterfaceProperties                 IP(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                         Cx(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    cout << "Initialization stage! Done!" << endl;

    if(RTC.Restart)
    {
        cout << "Restart data being read!";

        Phi.Read(OPSettings,BC, RTC.tStart);
        Cx.Read(OPSettings,BC, RTC.tStart);
        Tx.Read(OPSettings,BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        double iRadius = 0.25*OPSettings.Grid.Nx;
        // single lamella
        Initializations::Fractional(Phi, 0, 1, 0.5*iRadius, BC);
        Initializations::Sphere(Phi, 2, iRadius, 0.5*(OPSettings.Grid.Nx-1), 0.5*(OPSettings.Grid.Ny-1), 0, BC);

        // two lamellae
        /*Initializations::Fractional(Phi, 0, 1, iRadius/2.0, BC, OPSettings);
        Initializations::Sphere(Phi, 2, iRadius/2, (OPSettings.Nx)/4, (OPSettings.Ny)/2, 0, BC, OPSettings);
        Initializations::Sphere(Phi, 2, iRadius/2, (OPSettings.Nx)*3/4, (OPSettings.Ny)/2, 0, BC, OPSettings);
        */
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
    }

    cout << "Entering the Time Loop!!!" << endl;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
Timer.SetStart();
        dG.Clear();
Timer.SetTimeStamp("dG.Clear");
        IP.Set(Phi, BC);
Timer.SetTimeStamp("IP.Set");
        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
Timer.SetTimeStamp("DO.CalculatePhaseFieldIncrements");
        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
Timer.SetTimeStamp("DF.GetDrivingForce");
        dG.Average(Phi, BC);
Timer.SetTimeStamp("dG.Average");
        dG.MergePhaseFieldIncrements(Phi, IP);
Timer.SetTimeStamp("dG.MergePhaseFieldIncrements");
        Phi.NormalizeIncrements(BC, RTC.dt);
Timer.SetTimeStamp("Phi.NormalizeIncrements");
        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
Timer.SetTimeStamp("DF.Solve");
        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
Timer.SetTimeStamp("Tx.Set");
        Phi.MergeIncrements(BC, RTC.dt);
Timer.SetTimeStamp("Phi.MergeIncrements");

        /*if (Phi.Fields((OPSettings.Nx)/4, (OPSettings.Ny)/2, (OPSettings.Nz)/2)[0] <= 0.4)
        {// Moving frame
            Phi.MoveFrame(0,0,1, BC);
            Tx.MoveFrame(0,0,1, BC);
            Cx.MoveFrame(0,0,1, BC);
        }*/
        if (RTC.WriteVTK())
        {// Output to file in VTK format
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteStatistics(OPSettings,RTC.tStep, RTC.SimulationTime);
            dG.WriteVTK(OPSettings,Phi, RTC.tStep);
            //Tx.WriteVTK(tStep);
        }

        if (RTC.WriteRawData())
        {// Output to file of raw data
            Phi.Write(OPSettings,RTC.tStep);
            Cx.Write(OPSettings,RTC.tStep);
            Tx.Write(OPSettings,RTC.tStep);
        }
Timer.SetTimeStamp("File output");
        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", I_En);
                        message += ConsoleOutput::GetStandard("Interface energy", DO.Energy(Phi, IP));
            ConsoleOutput::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            Timer.PrintWallClockSummary();
        }
    } //end time loop
    return 0;
}
