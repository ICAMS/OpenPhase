
#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "Initializations.h"
#include "BoundaryConditions.h"
#include "DrivingForce.h"

using namespace std;
using namespace openphase;

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                        OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                  RTC(OPSettings);
    PhaseField                      Phi(OPSettings);
    DoubleObstacle                  DO(OPSettings);
    InterfaceProperties             IP(OPSettings);
    BoundaryConditions              BC(OPSettings);
    DrivingForce                    DF(OPSettings);

    Initializations::Young3(Phi, 0, 1, 2, 3, BC);

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.TimeStep = RTC.StartTimeStep; RTC.TimeStep <= RTC.MaxTimeStep; RTC.IncrementTimeStep())
    {
        DF.Clear();
        IP.Set(Phi, BC);
        DO.CalculatePhaseFieldIncrements(Phi, IP, DF);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);

        if (RTC.WriteVTK())
        {
           //  Output to file
           Phi.WriteVTK(OPSettings, RTC.TimeStep);
           IP.WriteVTK(OPSettings,RTC.TimeStep);
        }

        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", I_En);
                        message += ConsoleOutput::GetStandard("Interface energy", DO.Energy(Phi, IP));
            ConsoleOutput::WriteTimeStep(RTC, message);
        }
    } //end of time loop
    return 0;
}
