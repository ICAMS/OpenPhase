#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "Nucleation.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "BoundaryConditions.h"
#include "Initializations.h"
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
    BoundaryConditions                  BC(OPSettings);
    Nucleation                          Nuc(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    if(RTC.Restart)
    {
        cout << "Restart data being read!" << endl;
        cout << "Restart time step: " << RTC.tStart << endl;

        Phi.Read(OPSettings, BC, RTC.tStart);
        Cx.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC, RTC.tStart);

        cout << "Done reading restart parameters!" << endl;
    }
    else
    {
        int idx0 = Initializations::Single(Phi, 0, BC);
        Phi.FieldsProperties[idx0].State = AggregateStates::Liquid;

        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
    }

    cout << "Initialization stage done!" << endl;
    cout << "Entering the Time Loop!!!" << endl;

    IP.Set(Phi, BC);

    //RTC.dt = min(DF.ReportMaximumTimeStep(Tx), IP.ReportMaximumTimeStep());
    cout << "Time step: " << RTC.dt << endl;

    //-------------- The Time Loop -------------//
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();

        dG.Clear();
        Timer.SetTimeStamp("Clear driving force");

        Nuc.GenerateNucleationSites(Phi, Tx);
        Timer.SetTimeStamp("Generate Nucleation Sites");

        Nuc.PlantNuclei(Phi, RTC.tStep);
        Timer.SetTimeStamp("Plant Nuclei");

        IP.Set(Phi, BC);
        //DF.CalculateInterfaceMobility(Phi, Tx, BC, IP);
        Timer.SetTimeStamp("Set interface properties");

        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        Timer.SetTimeStamp("Get PsiDot");

        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("Chemical driving force");

        dG.Average(Phi, BC);
        Timer.SetTimeStamp("Driving force average");

        if (RTC.WriteVTK())
        {
            dG.WriteVTK(OPSettings, Phi, RTC.tStep);
        }
        Timer.SetTimeStamp("Write driving force");

        Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);
        Timer.SetTimeStamp("Check Nuclei");

        dG.MergePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Driving force merge to Psi");

        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Normalize Psi");

        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
        Timer.SetTimeStamp("Solve diffusion");

        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);

        Timer.SetTimeStamp("Set temperature");
        Phi.MergeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Merge phase fields");

        //  Write VTK files
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteVTK(OPSettings, RTC.tStep);
            Tx.WriteVTK(OPSettings, RTC.tStep);
            IP.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteStatistics(OPSettings, RTC.tStep, RTC.SimulationTime);
            Nuc.WriteStatistics(OPSettings, RTC.tStep);
        }
        // Write raw data
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
            std::string message = ConsoleOutput::GetStandard("Interface energy density", to_string(I_En));
            ConsoleOutput::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();
            Timer.PrintWallClockSummary();
        }
    } //end time loop

    return 0;
}
