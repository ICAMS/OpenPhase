
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
    Settings                    OPSettings;
    OPSettings.ReadInput();

    RunTimeControl              RTC(OPSettings);
    PhaseField                  Phi(OPSettings);
    DoubleObstacle              DO(OPSettings);
    InterfaceProperties         IP(OPSettings);
    BoundaryConditions          BC(OPSettings);
    DrivingForce                DF(OPSettings);

    Initializations::Young4(Phi, 0, 1, 2, 3, BC);

    int x1 = (OPSettings.Grid.Nx)/2;
    int y1 = (OPSettings.Grid.Ny)/2;
    int z1 = (OPSettings.Grid.Nz)/2;

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        DF.Clear();
        IP.Set(Phi, BC);
        DO.CalculatePhaseFieldIncrements(Phi, IP, DF);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);

        // Keeping the quadruple junction in the center of the simulation box
        double x2 = 0;
        double y2 = 0;
        double z2 = 0;
        double count = 1.0;
        for(int i = 0;i < OPSettings.Grid.Nx; ++i)
        for(int j = 0;j < OPSettings.Grid.Ny; ++j)
        for(int k = 0;k < OPSettings.Grid.Nz; ++k)
        if(Phi.Fields(i,j,k).size() == 4)
        {
            count ++;
            x2 += i;
            y2 += j;
            z2 += k;
        }
        if(count)
        {
            x2 /= count;
            y2 /= count;
            z2 /= count;
        }

        double tempdx = x2 - x1;
        double tempdy = y2 - y1;
        double tempdz = z2 - z1;

        int dx = 0;
        int dy = 0;
        int dz = 0;
        if (tempdx >= 1) dx = 1;
        if (tempdy >= 1) dy = 1;
        if (tempdz >= 1) dz = 1;
        if (tempdx <= -1) dx = -1;
        if (tempdy <= -1) dy = -1;
        if (tempdz <= -1) dz = -1;
        if (abs(dx) + abs(dy) + abs(dz) != 0)
        {
            Phi.MoveFrame(dx, dy, dz, BC);
        }

        if (RTC.WriteVTK())
        {
           //  Output to file
           Phi.WriteVTK(OPSettings, RTC.tStep);
        }

        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = ConsoleOutput::GetStandard("Interface energy density [J/m^3]", I_En);
                        message += ConsoleOutput::GetStandard("Interface energy [J]", DO.Energy(Phi, IP));
            ConsoleOutput::WriteTimeStep(RTC, message);

            //  For control (values at quadruple junction)
            Phi.PrintPointStatistics(x1,y1,z1);
            Phi.PrintPointStatistics(int(x2),int(y2),int(z2));
        }
    } //end of time loop
    return 0;
}
