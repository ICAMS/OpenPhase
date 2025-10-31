#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
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
    InterfaceProperties                 IP(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                         Cx(OPSettings);
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
        Cx.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC, RTC.tStart);
        Nuc.Read(OPSettings, BC, RTC.tStart);
        cout << "Done!" << endl;
    }
    else
    {
        RTC.tStart = -1;
        int idx0 = Initializations::Single(Phi, 0, BC);
        Phi.FieldsProperties[idx0].State = AggregateStates::Liquid;

        int idx1 = Phi.PlantGrainNucleus(1, (OPSettings.Grid.Nx)/2,
                                            (OPSettings.Grid.Ny)/2,
                                            (OPSettings.Grid.Nz)/2);
        //int idx1 = Phi.PlantGrainNucleus(1, 0, 0, 0);
        Phi.FieldsProperties[idx1].State = AggregateStates::Solid;

        //Initializations::RandomNuclei(Phi, OPSettings, 1, 5);
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
    }

    cout << "Initialization stage done!" << endl;
    cout << "Entering the Time Loop!!!" << endl;

    //-------------- The Time Loop -------------//

    IP.Set(Phi, Tx, BC);

    RTC.dt = min(DF.ReportMaximumTimeStep(Tx), IP.ReportMaximumTimeStep());
    cout << "Time step: " << RTC.dt << endl;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        dG.Clear();
        Timer.SetTimeStamp("Clear Driving Force");
        HS.Activate(Phi, Tx, RTC);
        Nuc.GenerateNucleationSites(Phi, Tx);
        Nuc.PlantNuclei(Phi, RTC.tStep);
        Timer.SetTimeStamp("Plant Nuclei");

        IP.Set(Phi, Tx, BC);
        DF.CalculateInterfaceMobility(Phi, Tx, BC, IP);
        Timer.SetTimeStamp("Calculate Interface Properties");

        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);

        Timer.SetTimeStamp("Curvature Contribution");

        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);

        Timer.SetTimeStamp("Chemical Driving Force");
        dG.Average(Phi, BC);

        Timer.SetTimeStamp("Driving Force Average");
        Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);
        Timer.SetTimeStamp("Check Nuclei");
        if (RTC.WriteVTK()) dG.WriteVTK(OPSettings, Phi, RTC.tStep);
        Timer.SetTimeStamp("Write Driving Force");

        dG.MergePhaseFieldIncrements(Phi, IP);

        Timer.SetTimeStamp("Merge Driving Force");
        Phi.NormalizeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Normalize Phase-field Increments");
        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);

        Timer.SetTimeStamp("Solve chemical diffusion");
        //Tx.Set(BC, RTC.dt);
        HS.Apply(Phi, Tx, HD);
        HD.SetEffectiveProperties(Phi, Tx);
        HD.SolveImplicit(Phi,BC,Tx, RTC.dt);

        Timer.SetTimeStamp("Solve heat diffusion");

        Phi.MergeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Merge Phase Fields");

        //  Output to file
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings,RTC.tStep);
            Cx.WriteVTK(OPSettings,RTC.tStep);
            Tx.WriteVTK(OPSettings,RTC.tStep);
            Cx.WriteStatistics(OPSettings, RTC.tStep, RTC.SimulationTime);
        }
        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
            Cx.Write(OPSettings, RTC.tStep);
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

            Tx.PrintStatistics();
            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();
            Timer.PrintWallClockSummary();
        }
    } //end time loop
    return 0;
}
