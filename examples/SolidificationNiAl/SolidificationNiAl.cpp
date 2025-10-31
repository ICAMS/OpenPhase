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
#include "Nucleation.h"
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
    InterfaceProperties                 SigmaMu(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                         Cx(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    BoundaryConditions                  BC_C(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");
    Nucleation                          Nuc(OPSettings);

    //BC_C.BCNZ = Fixed;

    if(RTC.Restart)
    {
        cout << "Restart data being read! ";

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

        Initializations::RandomNuclei(Phi, OPSettings, 1, 5);
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        //Tx.SetInitial(BC, Phi, 0);
    }

    cout << "Initialization stage done!" << endl;
    int x1 = 0.0;
    int y1 = 0.0;
    int z1 = 0.0;

    for(unsigned int n = 1; n < Phi.FieldsProperties.size(); n++)
    if(Phi.FieldsProperties[n].Exist)
    {
        x1 = Phi.FieldsProperties[n].Rcm[0];
        y1 = Phi.FieldsProperties[n].Rcm[1];
        z1 = Phi.FieldsProperties[n].Rcm[2];
        break;
    }
    // 1D diffusion domain extension
    /*
    int NzEXT = 4*OPSettings.Nz;
    vector<double> cZ_1D(NzEXT+2);
    vector<double> delta_cZ_1D(NzEXT+2,0.0);
    double cZ_ave = 0.0;
    for(int i = 0; i < OPSettings.Nx; i++)
    for(int j = 0; j < OPSettings.Ny; j++)
    {
        cZ_ave += Cx.Total(i,j,OPSettings.Nz-1)({0});
    }
    cZ_ave /= double(OPSettings.Nx*OPSettings.Ny);
    for(int k = 0; k < NzEXT+2; k++)
    {
        cZ_1D[k] = cZ_ave;
    }*/
    // End of 1D diffusion domain extension

    cout << "Entering the Time Loop!!!" << endl;

    //-------------- The Time Loop -------------//

    ofstream fit_out;
    RTC.dt = DF.ReportMaximumTimeStep(Tx);

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();

        dG.Clear();
        Nuc.GenerateNucleationSites(Phi, Tx);
        Nuc.PlantNuclei(Phi, RTC.tStep);
        Timer.SetTimeStamp("Plant Nuclei");

        SigmaMu.Set(Phi, Tx, BC);
        DF.CalculateInterfaceMobility(Phi, Tx, BC, SigmaMu);
        Timer.SetTimeStamp("Calculate Interface Properties");

        double I_En = 0.0;
        if (RTC.WriteToScreen()) I_En = DO.Energy(Phi, SigmaMu);

        DO.CalculatePhaseFieldIncrements(Phi, SigmaMu, dG);

        Timer.SetTimeStamp("Curvature Contribution");

        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);

        Timer.SetTimeStamp("Chemical Driving Force");
        dG.Average(Phi, BC);

        Timer.SetTimeStamp("Driving Force Average");
        Nuc.CheckNuclei(Phi, SigmaMu, dG, RTC.tStep);
        Timer.SetTimeStamp("Check Nuclei");
        if (RTC.WriteVTK()) dG.WriteVTK(OPSettings,Phi, RTC.tStep);

        dG.MergePhaseFieldIncrements(Phi, SigmaMu);

        Timer.SetTimeStamp("Merge Driving Force");
        Phi.NormalizeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Normalize Phase-field Increments");
        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
        //DF.Solve(Phi, Cx, Tx, BC_C, RTC.dt);
        // 1D diffusion domain extension
        /*{
            int Nx = OPSettings.Nx;
            int Ny = OPSettings.Ny;
            int Nz = OPSettings.Nz;

            double dt = OPSettings.dt;
            double dx_2 = 1.0/(OPSettings.dx*OPSettings.dx);

            cZ_1D[0] = 0.0;

            for(int i = 0; i < Nx; i++)
            for(int j = 0; j < Ny; j++)
            {
                cZ_1D[0] += Cx.Total(i,j,Nz-1)({0});
            }

            cZ_1D[0] /= double(Nx*Ny);

            for(int k = 1; k < NzEXT; k++)
            {
                delta_cZ_1D[k] = DF.DC(Nx-1,Ny-1,Nz-1)({0})*(cZ_1D[k-1] - 2.0*cZ_1D[k] + cZ_1D[k+1])*dx_2;
            }

            for(int k = 1; k < NzEXT; k++)
            {
                cZ_1D[k] += delta_cZ_1D[k]*dt;
            }

            for(int i = -1; i < Nx+1; i++)
            for(int j = -1; j < Ny+1; j++)
            {
                Cx.Total(i,j,OPSettings.Nz)({0}) = cZ_1D[1];
                Cx.Phase(i,j,OPSettings.Nz)({0,0}) = cZ_1D[1];
            }
        }*/
        // End of 1D diffusion domain extension
        Timer.SetTimeStamp("Solve diffusion");
        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
        //Tx.Set(BC, Phi, 6.2e8, 1.773e6, 0, OPSettings.dt);

        Timer.SetTimeStamp("Set temperature");
        Phi.MergeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Merge Phase Fields");

        //  Output to file
        if (RTC.WriteVTK())
        {
            // Write data in VTK format
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteVTK(OPSettings, RTC.tStep);
            Tx.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteStatistics(OPSettings, RTC.tStep, RTC.SimulationTime);
            //Mu.WriteVTK(tStep, 0, 3);
        }
        if (RTC.WriteRawData())
        {
            // Write raw data
            Phi.Write(OPSettings, RTC.tStep);
            Cx.Write(OPSettings, RTC.tStep);
            Tx.Write(OPSettings, RTC.tStep);
            Nuc.Write(OPSettings, RTC.tStep);
        }
        Timer.SetTimeStamp("File Output");
        time_t rawtime;
        time(&rawtime);

        //  Output to screen
        if(RTC.WriteToScreen())
        {
            cout << "==================================\n"
                 << "Time Step        : " << RTC.tStep << "\n"
                 << "Wall Clock Time  : " << ctime(&rawtime)
                 << "----------------------------------\n"
                 << "Interface Energy : " << I_En <<  "\n"
                 << "==================================\n" << endl;

            //  Statistics
            Phi.PrintPointStatistics(x1,y1,z1);
            Cx.PrintPointStatistics(x1,y1,z1);
            Tx.PrintPointStatistics(x1,y1,z1);

            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            //Phi.PrintVolumeFractions(OPSettings.ElementNames);
            Timer.PrintWallClockSummary();
        }
    } //end time loop
    fit_out.close();

    return 0;
}
