#include "Settings.h"
#include "RunTimeControl.h"
#include "PhaseField.h"
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
    Temperature                     Tx(OPSettings);
    BoundaryConditions              BC(OPSettings);
    HeatDiffusion                   HD(OPSettings);
    TimeInfo                        Timer(OPSettings, "Execution Time Statistics");

    cout << "Initialization stage! Done!" << endl;

    int x1 = OPSettings.Grid.Nx/2;
    int y1 = OPSettings.Grid.Ny/2;
    int z1 = OPSettings.Grid.Nz/2;
	
    int matrixIdx = Initializations::Single(Phi, 0, BC);
    EulerAngles ph1({0,0,0},XYZ);

    Phi.FieldsProperties[matrixIdx].Orientation = ph1.getQuaternion();
    Tx.SetInitial(BC);
    Tx(x1,y1,z1) += 100.0;
    cout << "Entering the Time Loop!!!" << endl;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteRawData())
        {   
            Tx.Write(OPSettings, RTC.tStep);
        }
        if (RTC.WriteVTK())
        {
            Tx.WriteVTK(OPSettings, RTC.tStep);
        }
        Timer.SetStart();
        Tx.Set(BC,Phi, RTC.SimulationTime, RTC.dt);
        HD.SolveImplicit(Phi, BC, Tx, RTC.dt);
    Timer.SetTimeStamp("Set Temperature");
    Timer.SetTimeStamp("File Output");
        if (RTC.WriteToScreen())
        {
            time_t rawtime;
            time(&rawtime);
            //  Output to screen
            cout << "+++++++++++++++++++++++++++++++++++++++++\n"
                 << "Time Step:          " << RTC.tStep << "\n"
                 << "Wall Clock Time:    " << ctime(&rawtime) << "\n";
        }
    } //end time loop
   return 0;
}
