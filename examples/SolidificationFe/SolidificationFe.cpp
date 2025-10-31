#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "UserDrivingForce.h"
#include "Temperature.h"
#include "HeatDiffusion.h"
#include "HeatSources.h"
#include "BoundaryConditions.h"
#include "Initializations.h"
#include "Nucleation.h"
#include "Tools/TimeInfo.h"

using namespace std;
using namespace openphase;
/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                            OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                      RTC(OPSettings);
    PhaseField                          Phi(OPSettings);
    DoubleObstacle                      DO(OPSettings);
    UserDrivingForce                    UDF(OPSettings);
    InterfaceProperties                 IP(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    HeatDiffusion                       HD(OPSettings);
    HeatSources                         HS(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    Nucleation                          Nuc(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    if(RTC.Restart)
    {
        cout << "Restart data being read... ";

        Phi.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC, RTC.tStart);
        Nuc.Read(OPSettings, BC, RTC.tStart);
        cout << "Done!" << endl;
    }
    else
    {
        int idx0 = Initializations::Single(Phi, 0, BC);

        Tx.SetInitial(BC);
    }

    cout << "Initialization stage done!" << endl;
    cout << "Entering the Time Loop!!!" << endl;

    //-------------- The Time Loop -------------//

    IP.Set(Phi, Tx, BC);

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        dG.Clear();
        HS.Activate(Phi, Tx, RTC);
        Nuc.GenerateNucleationSites(Phi, Tx);
        Nuc.PlantNuclei(Phi, RTC.tStep);
        Timer.SetTimeStamp("Plant Nuclei");

        IP.Set(Phi, Tx, BC);
        Timer.SetTimeStamp("Calculate Interface Properties");

        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);

        Timer.SetTimeStamp("Curvature Contribution");
        UDF.SetDrivingForce(Phi, dG, Tx);
        Timer.SetTimeStamp("Chemical Driving Force");
        dG.Average(Phi, BC);

        Timer.SetTimeStamp("Driving Force Average");
        Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);
        Timer.SetTimeStamp("Check Nuclei");
        if (RTC.WriteVTK()) dG.WriteVTK(OPSettings, RTC.tStep, 1, 0);
        Timer.SetTimeStamp("Write Driving Force");

        dG.MergePhaseFieldIncrements(Phi, IP);

        Timer.SetTimeStamp("Merge Driving Force");
        Phi.NormalizeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Normalize Phase-field Increments");

        HS.Apply(Phi, Tx, HD);
        HD.SetEffectiveProperties(Phi, Tx);
        HD.SolveImplicit(Phi,BC,Tx, RTC.dt);
        Timer.SetTimeStamp("Solve heat diffusion");

        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge Phase Fields");

        //  Output to file
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Tx.WriteVTK(OPSettings, RTC.tStep);
        }
        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
            Tx.Write(OPSettings, RTC.tStep);
            Nuc.Write(OPSettings, RTC.tStep);
        }
        Timer.SetTimeStamp("File Output");

        //  Output to screen
        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message = ConsoleOutput::GetStandard("Interface energy density", I_En);
            ConsoleOutput::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();
            Tx.PrintStatistics();
            Timer.PrintWallClockSummary();
        }
    } //end time loop
    return 0;
}
