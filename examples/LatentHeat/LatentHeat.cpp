#include "Settings.h"
#include "RunTimeControl.h"
#include "PhaseField.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "Temperature.h"
#include "Initializations.h"
#include "BoundaryConditions.h"
#include "HeatDiffusion.h"
#include "Tools/TimeInfo.h"

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
    Temperature                     Tx(OPSettings);
    BoundaryConditions              BC(OPSettings);
    HeatDiffusion                   HD(OPSettings);
    TimeInfo                        Timer(OPSettings, "Execution Time Statistics");

    cout << "Initialization stage! Done!" << endl;

    RTC.tStart = -1;

    int x1 = OPSettings.Grid.Nx/2;
    int y1 = OPSettings.Grid.Ny/2;
    int z1 = OPSettings.Grid.Nz/2;

    int matrixIdx = Initializations::Single(Phi, 0, BC);
    int sphereIdx = Initializations::Sphere(Phi, 1, 10, x1, y1, z1, BC);
    EulerAngles ph1({0,0,0},XYZ);

    Phi.FieldsProperties[matrixIdx].Orientation = ph1.getQuaternion();
    Phi.FieldsProperties[sphereIdx].Orientation = ph1.getQuaternion();
    Tx.SetInitial(BC);

    cout << "Entering the Time Loop!!!" << endl;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteRawData())
        {   
            Tx.Write(OPSettings, RTC.tStep);
        }
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Tx.WriteVTK(OPSettings, RTC.tStep);
        }
        IP.Set(Phi, BC);
        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetStart();
        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
        HD.SolveImplicit(Phi, BC, Tx, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);
        if (RTC.WriteToScreen())
        {
            std::string message  = "";
            ConsoleOutput::WriteTimeStep(RTC, message);
        }
    } //end time loop
   return 0;
}
