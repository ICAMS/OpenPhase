#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "Initializations.h"
#include "Mechanics.h"
#include "BoundaryConditions.h"
#include "Tools/TimeInfo.h"
#include "PhysicalConstants.h"

using namespace std;
using namespace openphase;
/*
void SetDiffusionCoefficientsCx(PhaseField& Phase, Temperature& Tx,
                                EquilibriumPartitionDiffusionBinary& DF,
                                Composition& Elements)                          ///<  Sets concentration dependent diffusion coefficients in each point
{
    size_t Comp = 0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DF.DC,DF.DC.Bcells(),)
    {
        double invRT  = 1.0/(PhysicalConstants::R*Tx(i,j,k));
        double T = Tx(i,j,k);
        for (size_t n = 0; n < Phase.Nphases; ++n)
        {
            if (n == 0)
            {
                double locCx = Elements.MoleFractions(i,j,k)({0,Comp});
                DF.DC(i,j,k)({n}) = 4.45e-7*(1.0 + (locCx/(1.0 - locCx)) *
                                            (1.0 - (locCx/(1.0 - locCx))) * (8339.9/T)) *
                exp(-((1.0/T) - 2.221e-4)*(17767.0 - (locCx/(1.0 - locCx)) * 26436.0));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ScaleDiffusionCoefficientsStress(PhaseField& Phase, Temperature& Tx,
                                  EquilibriumPartitionDiffusionBinary& DF,
                                  ElasticProperties& EP, double A)              ///<  Scales diffusion coefficient according to its stress dependence
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DF.DC,DF.DC.Bcells(),)
    {
        double invRT  = 1.0/(PhysicalConstants::R*Tx(i,j,k));

        for (size_t n = 0; n < Phase.Nphases; ++n)
        {
            if (n == 0)
            {
                DF.DC(i,j,k)({n}) *=

                (1.0 - (A*EP.Stresses(i,j,k).Pressure()/150e6))*
                exp(-30000.0*100.0*EP.Stresses(i,j,k).Pressure()*invRT);

                // sets minimum diffusion coefficient
                DF.DC(i,j,k)({n}) = std::max(3.3e-14, DF.DC(i,j,k)({n}));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}*/

void SetEffectiveShearFreeElasticConstants(PhaseField& Phase, ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.EffectiveElasticConstants,EP.EffectiveElasticConstants.Bcells(),)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            EP.EffectiveElasticConstants(i,j,k).set_to_zero();
            dMatrix6x6 TempCompliances;
            TempCompliances.set_to_zero();
            bool gamma_present = false;
            bool alpha_present = false;
            bool cementite_present = false;
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                TempCompliances += (EP.ElasticConstants(i,j,k).get_value(alpha->index).inverted() * alpha->value);
                if(Phase.FieldsProperties[alpha->index].Phase == 0) gamma_present = true;
                if(Phase.FieldsProperties[alpha->index].Phase == 1) alpha_present = true;
                if(Phase.FieldsProperties[alpha->index].Phase == 2) cementite_present = true;
            }
            if(gamma_present and alpha_present)

            {
                TempCompliances(3, 3) = 8.6e-11;
                TempCompliances(4, 4) = 8.6e-11;
                TempCompliances(5, 5) = 8.6e-11;
                TempCompliances(0, 1) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(0, 2) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(1, 2) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(1, 0) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(2, 0) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(2, 1) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
            }
            if(gamma_present and cementite_present)
            {
                TempCompliances(3, 3) = 10.42e-11;
                TempCompliances(4, 4) = 10.42e-11;
                TempCompliances(5, 5) = 10.42e-11;
                TempCompliances(0, 1) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(0, 2) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(1, 2) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(1, 0) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(2, 0) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(2, 1) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
            }
            EP.EffectiveElasticConstants(i,j,k) = TempCompliances.inverted();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                        OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                      RTC(OPSettings);
    PhaseField                          Phi(OPSettings);
    DoubleObstacle                      DO(OPSettings);
    InterfaceProperties                 IP(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                         Cx(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    ElasticProperties                   EP(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    ElasticitySolverSpectral            ES(OPSettings, BC);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    cout << "Initialization stage! Done!" << endl;

    if(RTC.Restart)
    {
        cout << "Restart data being read!";
        Phi.Read(OPSettings, BC, RTC.tStart);
        Cx.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        Initializations::Single(Phi, 0, BC);
        Initializations::Layer(Phi, 1, {0.,0.,0.0*OPSettings.Grid.Nz/2}, {0.,0.,1.}, OPSettings.Grid.Nx/10, BC);
        Initializations::Sphere(Phi, 2, OPSettings.Grid.Nx/12,
                                        OPSettings.Grid.Nx/2,
                                        OPSettings.Grid.Ny/2,
                                        0.0*OPSettings.Grid.Nz/2, BC);
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        EP.SetEffectiveProperties(Phi, Cx);
    }

    int x1 = (OPSettings.Grid.Nx+1)/2;
    int y1 = (OPSettings.Grid.Ny+1)/2;
    int z1 = 13;

    double tStepOld = RTC.tStart;
    double AusteniteVolumeOld = Phi.FieldsProperties[0].Volume;

    ofstream GrowthVelocity("GrowthVelocity.dat", ios::out);
    GrowthVelocity << "time(s) \t Velocity(m/s)" << endl;
    GrowthVelocity.close();

    //RTC.dt = DF.ReportMaximumTimeStep();
    cout << "Time step: " << RTC.dt << endl;

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        dG.Clear();
        Timer.SetTimeStamp("Clear Driving Force");
        IP.Set(Phi, Tx, BC);
        Timer.SetTimeStamp("Set Interface Properties");
        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        Timer.SetTimeStamp("Get PhiDot");
        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("Chemical Driving Force");
        EP.SetEffectiveProperties(Phi, Cx);
        Timer.SetTimeStamp("Elastic Properties");
        EP.CalculateChemicalPotentialContribution(Phi, DF);
        Timer.SetTimeStamp("Elastic Chemical Potential");
        EP.CalculateDrivingForce(Phi, dG);
        Timer.SetTimeStamp("Elastic Driving Force");
        if(!(RTC.tStep % 10))
        {
            ES.Solve(EP, BC, RTC.dt);
        }
        Timer.SetTimeStamp("Elastic Solver");
        dG.Average(Phi, BC);
        Timer.SetTimeStamp("Driving Force Average");
        dG.MergePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge Velocities");

        //Setting up modified diffusion coefficient and elastic constants

/*1*/  //SetDiffusionCoefficientsCx(Phi, Tx, DF, Cx);
/*2*/   //ScaleDiffusionCoefficientsStress(Phi, Tx, DF, EP, 10);
/*3*/   //SetEffectiveShearFreeElasticConstants(Phi, EP);

        // End setting up modified diffusion coefficient and elastic constants

        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
        Timer.SetTimeStamp("Diffusion Solver");
        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge PhaseFields");

        if (Phi.Fields((OPSettings.Grid.Nx+1)/2, (OPSettings.Grid.Ny+1)/2, (OPSettings.Grid.Nz+1)/3).get_value(0) <= 0.4)
        {
            double deltaV = Phi.FieldsProperties[0].Volume - AusteniteVolumeOld;
            double deltaL = deltaV*OPSettings.Grid.dx/(OPSettings.Grid.Nx*OPSettings.Grid.Ny);

            double vel = -0.5*deltaL/((RTC.tStep - tStepOld)*RTC.dt);

            GrowthVelocity.open("GrowthVelocity.dat", ios::app);
            GrowthVelocity << RTC.tStep*RTC.dt << "\t" << vel << endl;
            GrowthVelocity.close();

            Phi.ConsumePlane(0,0,1, (OPSettings.Grid.Nx+1)/2, (OPSettings.Grid.Ny+1)/2, (OPSettings.Grid.Nz+1)/2, BC);
            Tx.ConsumePlane(0,0,1, (OPSettings.Grid.Nx+1)/2, (OPSettings.Grid.Ny+1)/2, (OPSettings.Grid.Nz+1)/2, BC);
            Cx.ConsumePlane(0,0,1, (OPSettings.Grid.Nx+1)/2, (OPSettings.Grid.Ny+1)/2, (OPSettings.Grid.Nz+1)/2, BC);

            AusteniteVolumeOld = Phi.FieldsProperties[0].Volume;
            tStepOld = RTC.tStep;
        }

        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
        Timer.SetTimeStamp("Set Temperature");
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteVTK(OPSettings, RTC.tStep);
            EP.WriteStressesVTK(OPSettings, RTC.tStep);
            EP.WriteTotalStrainsVTK(OPSettings, RTC.tStep);
            Cx.WriteStatistics(OPSettings, RTC.tStep, RTC.SimulationTime);
        }
        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
            Cx.Write(OPSettings, RTC.tStep);
            Tx.Write(OPSettings, RTC.tStep);
        }
        Timer.SetTimeStamp("File Output");
        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            double E_En = EP.AverageEnergyDensity(Phi);
            double T_En = E_En + I_En;
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", I_En);
                        message += ConsoleOutput::GetStandard("Elastic energy density", E_En);
                        message += ConsoleOutput::GetStandard("Total energy density", T_En);

            ConsoleOutput::WriteTimeStep(RTC, message);

            Phi.PrintPointStatistics(x1,y1,z1);

            //Cx.PrintPointStatistics(x1,y1,z1);
            //Tx.PrintPointStatistics(x1,y1,z1);
            //EP.PrintPointStatistics(x1,y1,z1);
            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            Phi.PrintVolumeFractions();
            Timer.PrintWallClockSummary();
        }
    } //end time loop

   return 0;
}
