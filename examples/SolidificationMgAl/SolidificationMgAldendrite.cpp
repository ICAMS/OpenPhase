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

    cout << "Initialization stage done!" << endl;

    if(RTC.Restart)
    {
        cout << "Restart data being read!";
        Phi.Read(OPSettings, BC, RTC.tStart);
        Cx.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC, RTC.tStart);
        Nuc.Read(OPSettings, BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        // Single liquid phase starting configuration
        size_t idx1 = Initializations::Single(Phi, 0, BC);
        // Adding the nucleus in the middle of the simulation domain
        size_t idx2 = Phi.PlantGrainNucleus(1, (OPSettings.Grid.Nx)/2,(OPSettings.Grid.Ny)/2,(OPSettings.Grid.Nz)/2);
        // Orienting the nucleus as desired
        EulerAngles locAngle({Pi/2.0, 0.0*Pi/8.0, Pi/8.0}, XYZ);
        //EulerAngles locAngle({5.0*Pi/8.0, 3.0*Pi/8.0, Pi/8.0}, XYZ);
        Phi.FieldsProperties[idx2].Orientation = locAngle.getQuaternion();

        // Liquid channel surrounded by alpha phase
        /*std::vector<size_t> ids = Initializations::TwoWalls(Phi, 0, 1, 10,BC,OPSettings);
        Phi.FieldsProperties[ids[0]].State = Liquid;
        Phi.FieldsProperties[ids[1]].State = Solid;*/

        // Liquid pocket surrounded by alpha phase
        /*size_t idx1 = Initializations::Single(Phi, 1, BC, OPSettings);
        Phi.FieldsProperties[idx1].State = AggregateStates::Solid;
        //size_t idx2 = Initializations::Sphere(Phi, 0, OPSettings.Nz/2 - 4, (OPSettings.Nx)/2,(OPSettings.Ny)/2,(OPSettings.Nz)/2, BC, OPSettings);
        size_t idx2 = Initializations::Ellipsoid(Phi, 0, OPSettings.Nx/2 - 5, OPSettings.Ny/2 - 5, OPSettings.Nz/2 - 5, (OPSettings.Nx)/2, (OPSettings.Ny)/2, (OPSettings.Nz)/2, BC, OPSettings);
        Phi.FieldsProperties[idx2].State = AggregateStates::Liquid;*/
        // Setting initial values for composition and temperature
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
    }

    fstream tt_out("TimeTemperature.dat", ios::out);
    tt_out << "Time    dt    Temperature" << endl;
    tt_out.close();

    cout << "Entering the Time Loop!!!" << endl;

    //-------------- The Time Loop -------------//
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        dG.Clear();
        Nuc.GenerateNucleationSites(Phi, Tx);
        Nuc.PlantNuclei(Phi, RTC.tStep);
        IP.Set(Phi,Tx, BC);
        DF.CalculateInterfaceMobility(Phi,Tx,BC,IP);
        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
        dG.Average(Phi, BC);
        Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);

        if (RTC.WriteVTK())
        {//  Write data in VTK format
            dG.WriteVTK(OPSettings, RTC.tStep, 1, 0);
            dG.WriteVTK(OPSettings, Phi, RTC.tStep);
        }

        dG.MergePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);
        //Tx.Set(BC, Phi, 6.2e8, 1.773e6, 1, dt);
        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);

        if (RTC.WriteRawData())
        {// Write raw data
            Phi.Write(OPSettings, RTC.tStep);
            Cx.Write(OPSettings, RTC.tStep);
            Tx.Write(OPSettings, RTC.tStep);
            Nuc.Write(OPSettings, RTC.tStep);
        }
        if (RTC.WriteVTK())
        {// Write data in VTK format
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Phi.WriteLaplacianVTK(OPSettings, RTC.tStep, 1);

            Cx.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteStatistics(OPSettings, RTC.tStep, RTC.SimulationTime);
        }
        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", I_En);
                        message += ConsoleOutput::GetStandard("Interface energy", DO.Energy(Phi, IP));
            ConsoleOutput::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();

            fstream tt_out("TimeTemperature.dat", ios::out | ios::app);
            tt_out << RTC.SimulationTime <<  "   " << RTC.dt <<  "   " << Tx(1,1,1) << endl;
            tt_out.close();
        }
    } //end time loop
    return 0;
}
