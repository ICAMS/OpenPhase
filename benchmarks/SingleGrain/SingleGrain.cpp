#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "Initializations.h"
#include "BoundaryConditions.h"
#include "Tools/CSVParser.h"
#include "Tools/TimeInfo.h"
#include "ConsoleOutput.h"
#include "DrivingForce.h"

using namespace std;
using namespace openphase;
int main()
{
    Settings                    OPSettings;
    OPSettings.ReadInput();

    RunTimeControl              RTC(OPSettings);
    PhaseField                  Phi(OPSettings);
    DoubleObstacle              DO(OPSettings);
    InterfaceProperties         IP(OPSettings);
    BoundaryConditions          BC(OPSettings);
    DrivingForce                DF(OPSettings);

    // Initialize phase-fields
    Initializations::Single(Phi, 0, BC);
    int locIndex = Initializations::Sphere(Phi, 1, OPSettings.Grid.Nx*0.4,
    (OPSettings.Grid.Nx)/2.0, (OPSettings.Grid.Ny)/2.0, (OPSettings.Grid.Nz)/2.0, BC);

    //-------------------------------------------------//

    double Vol0 = Phi.FieldsProperties[locIndex].Volume;
    double R0 = exp(log(0.75/Pi*Vol0)/3.0);
    double R02 = R0*R0;

    // Create output file
    std::vector<std::string> headernames {"tstep", "R^2(simulation)", "R^2(analytic)", "Ratio"};
    std::vector<double> dataInit {0, R02, R02, 1};
    CSVParser::WriteHeader("R_2_graph.dat", headernames);
    CSVParser::WriteData("R_2_graph.dat", dataInit);

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        DF.Clear();
        IP.Set(Phi, BC);
        DO.CalculatePhaseFieldIncrements(Phi, IP, DF);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);

        if (RTC.WriteVTK())
        {
            /// Write data in VTK format
            Phi.WriteVTK(OPSettings, RTC.tStep);
            /// Output of the grain radius
            double Vol = Phi.FieldsProperties[locIndex].Volume;
            double Radius = exp(log(0.75/Pi*Vol)/3.0);

            std::vector<double> datat;

            datat.push_back(double(RTC.tStep+1));
            datat.push_back(Radius*Radius);
            double dx = OPSettings.Grid.dx;
            double RadiusTheor = R02 - 4.0*IP.InterfaceMobility(0,1).MaxMobility*
                                           IP.InterfaceEnergy(0,1).MaxEnergy*(RTC.tStep+1)*
                                           RTC.dt/dx/dx;
            if (RadiusTheor >= 0) datat.push_back(RadiusTheor);
            else datat.push_back(0.0);
            if (Radius > 0) {datat.push_back((Radius*Radius)/(R02 - 4.0*IP.InterfaceMobility(0,1).MaxMobility*
                                                                        IP.InterfaceEnergy(0,1).MaxEnergy*(RTC.tStep+1)*
                                                                        RTC.dt/dx/dx));};
            CSVParser::WriteData("R_2_graph.dat", datat);
        }

        if (RTC.WriteToScreen())
        {
            int I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", I_En);
                        message += ConsoleOutput::GetStandard("Interface energy", DO.Energy(Phi, IP));
            ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps, message);
        }
    }
    return 0;
}
