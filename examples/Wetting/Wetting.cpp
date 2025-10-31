#include "BoundaryConditions.h"
#include "ConsoleOutput.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "Nucleation.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Temperature.h"
#include "Tools/TimeInfo.h"
#include "UserDrivingForce.h"
#include "Velocities.h"

using namespace openphase;

class LocalLBM : public FlowSolverLBM                                           ///< Modification of FLowSolverLBM which are only used in this example
{
 public:
    LocalLBM(Settings& locSettings, double in_dt): FlowSolverLBM(locSettings, in_dt){};

    void InitializeSingle();                                                    ///<  Initializes a single density
    void InitializeSphere(const double Radius, const std::array<double,3> pos,
            const double* rho = nullptr);
   // void InitializeY0(const double lbDensity);                                      ///<  Initializes a single sphere
};

double CalculatePhaseMass(PhaseField& Phase, size_t pIdx);

int main(int argc, char *argv[])
{

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    Settings OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                      RTC   (OPSettings);

    BoundaryConditions                  BC    (OPSettings);
    //Composition                         Cx    (OPSettings);
    DoubleObstacle                      DO    (OPSettings);
    DrivingForce                        dG    (OPSettings);
    //EquilibriumPartitionDiffusionBinary DF    (OPSettings);

    UserDrivingForce                    UDF(OPSettings);

    InterfaceProperties                 IP    (OPSettings);
    LocalLBM                            LB    (OPSettings, RTC.dt);
    Nucleation                          Nuc   (OPSettings);
    PhaseField                          Phi   (OPSettings);
    Temperature                         Tx    (OPSettings);
    TimeInfo                            Timer (OPSettings, "Execution Time Statistics");
    Velocities                          Vel   (OPSettings);

    // Set initial conditions
    if(RTC.Restart)
    {
        Phi.Read(OPSettings, BC, RTC.tStart);
        LB .Read(OPSettings, BC, RTC.tStart);
    }
    else
    {
        const double Nx = OPSettings.Grid.Nx;
        const double Ny = OPSettings.Grid.Ny;
        const double Nz = OPSettings.Grid.Nz;
        const double iWidth = Phi.Grid.iWidth;
        const size_t FluidPhaseIdx = 0;
        //const size_t CuPhaseIdx = 1;
        const size_t FePhaseIdx = 2;

        // Read initialization input parameters
        std::fstream inpF (DefaultInputFileName, std::ios::in);
        std::stringstream inp;
        inp << inpF.rdbuf();
        inpF.close();
        int moduleLocation = FileInterface::FindModuleLocation(inp, "Initialization");
        const double LiquidRadius = FileInterface::ReadParameterD(inp, moduleLocation, "dLiquidRadius")/OPSettings.Grid.dx;

        // Initialize vapor phase
        LB.InitializeSingle();
        size_t FluidGrainIdx = Initializations::Single(Phi, FluidPhaseIdx, BC);
        Phi.FieldsProperties[FluidGrainIdx].State = AggregateStates::Liquid;

        // Initialize liquid sphere
        const std::array<double,3> sphere_pos = {Nx/2,Ny/5+19,Nz/2};
        LB.InitializeSphere(LiquidRadius, sphere_pos , &LB.LiquidDensity[0]);

        // Initialize solid substrate
        size_t SolidGrainIdx = Initializations::Rectangular(Phi, FePhaseIdx, Nx+iWidth,Ny/5+iWidth,Nz+iWidth, Nx/2,Ny/10,Nz/2, BC);
        Phi.FieldsProperties[SolidGrainIdx].Mobile = false;
        LB.FinalizeInitialiation(Phi, Vel, BC);

        // Initialize Nucleation in liquid
//        for (int i = 0; i < OPSettings.Nx; i=i+5)
//        {
//            int nuc[OPSettings.Nx];
//            nuc[i] = Phi.PlantGrainNucleus(CuPhaseIdx, i, Ny/5+iWidth, Nz/2);
//        }
         //size_t nuc = Phi.PlantGrainNucleus(CuPhaseIdx, Nx/2, Ny/5+iWidth, Nz/2);
//       size_t nuc1 = Phi.PlantGrainNucleus(CuPhaseIdx, Nx/2+3, Ny/5+5, Nz/2);
//       size_t nuc2 = Phi.PlantGrainNucleus(CuPhaseIdx, Nx/2-3, Ny/5+5, Nz/2);
//       size_t nuc3 = Phi.PlantGrainNucleus(CuPhaseIdx, Nx/2+5, Ny/5+5, Nz/2);
//       size_t nuc4 = Phi.PlantGrainNucleus(CuPhaseIdx, Nx/2-5, Ny/5+5, Nz/2);
//         size_t nuc = Initializations::Sphere(Phi, CuPhaseIdx, 0.5, Nx/2, Ny/5+1, Nz/2, BC, OPSettings);
//         Phi.FieldsProperties[nuc].Density = L12Density;
//       Phi.FieldsProperties[nuc1].Density = L12Density;
//       Phi.FieldsProperties[nuc2].Density = L12Density;
//       Phi.FieldsProperties[nuc3].Density = L12Density;
//       Phi.FieldsProperties[nuc4].Density = L12Density;

        //Cx.SetInitialMoleFractions(Phi,1);
        Tx.SetInitial(BC);
    }

    //DF.SetDiffusionCoefficients(Phi, Tx);
    //IP.Set(Phi, Tx);
    //
    double InitialTotalMass = 0.0;

    std::vector<double> Forces;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteVTK())
        {
            //Cx .WriteVTK(RTC.tStep, OPSettings);
            LB .WriteVTK(OPSettings, Phi, RTC.tStep);
            Phi.WriteVTK(OPSettings, RTC.tStep);
            //Tx .WriteVTK(RTC.tStep, OPSettings);
        }

        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
            LB .Write(OPSettings, RTC.tStep);
        }

        if(RTC.WriteToScreen())
        {
            dG.PrintDiagnostics();
            ConsoleOutput::WriteTimeStep(RTC);
            ConsoleOutput::Write("Interface energy", DO.Energy(Phi, IP));
            ConsoleOutput::Write("");
            ConsoleOutput::Write("Liquid Volume", LB.CalculateLiquidVolume(0));
            ConsoleOutput::Write("Vapor Volume", LB.CalculateVaporVolume(0));
            ConsoleOutput::Write("");
            ConsoleOutput::Write("Liquid Mass", LB.CalculateLiquidMass(0));
            ConsoleOutput::Write("Vapor Mass", LB.CalculateVapordMass(0));
            ConsoleOutput::Write("");
            double FluidMass = LB.CalculateFluidMass()[0];
            double SolidMass = Phi.CalculateSolidMass();
            double TotalMass = FluidMass + SolidMass;
            if (RTC.tStep == 0) InitialTotalMass = TotalMass;
            double DeltaMass = TotalMass - InitialTotalMass;
            if (DeltaMass < std::numeric_limits<double>::epsilon()) DeltaMass = 0.0;
            ConsoleOutput::Write("Fluid Mass", FluidMass);
            ConsoleOutput::Write("Solid Mass", SolidMass);
            ConsoleOutput::Write("-");
            ConsoleOutput::Write("Total Mass", TotalMass);
            ConsoleOutput::Write("Initial Total Mass", InitialTotalMass);
            ConsoleOutput::Write("Mass Conservation Error ", DeltaMass);
            ConsoleOutput::Write("");
            LB.EnforceTotalMassConservation(Phi);

            Phi.PrintVolumeFractions();
            //Tx.PrintPointStatistics(OPSettings.Nx/2, 0,0);
            //Tx.PrintPointStatistics(OPSettings.Nx/2, 0,OPSettings.Ny-1);
            std::cout << "Mximum and Minimum Temperatures "<< std::endl;
            std::cout << "Tmax: "<< Tx.Tmax << std::endl;
            std::cout << "Tmin: "<< Tx.Tmin << std::endl;
            Timer.PrintWallClockSummary();
        }

        Timer.SetStart();
        dG.Clear();
        //HS.Activate(Phi, Tx, RTC);

        //Nuc.GenerateNucleationSites(Phi, Tx);
        //Nuc.PlantNuclei(Phi, RTC.tStep);
        //Timer.SetTimeStamp("Plant Nuclei");
        IP.Set(Phi, Tx, BC);
        Timer.SetTimeStamp("Calculate Interface Properties");
        //DF.CalculateInterfaceMobility(Phi, Cx, Tx, BC, IP);
        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        Timer.SetTimeStamp("Curvature Contribution");
        //DF.GetDrivingForce(Phi, Cx, Tx, dG);
        UDF.SetDrivingForce(Phi, dG, Tx);
        Timer.SetTimeStamp("User Driving Force");
        dG.Average(Phi, BC);
        Timer.SetTimeStamp("Driving Force Average");
        //Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);
        //Timer.SetTimeStamp("Check Nuclei");
        if (RTC.WriteVTK()) dG.WriteVTK(OPSettings, Phi, RTC.tStep);
        Timer.SetTimeStamp("Write Driving Force");
        dG.MergePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Merge Driving Force");
        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Normalize Phase-field Increments");
        //DF.Solve(Phi, Cx, Tx, BC, RTC.dt);
        //Timer.SetTimeStamp("Solve chemical diffusion");
        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
        Timer.SetTimeStamp("Set Temperature");
        //HS.Apply(Phi, Tx, HD);
        //HD.SetEffectiveProperties(Phi,Tx);
        //HD.SolveImplicit(Phi,BC,Tx, RTC.dt);
        //Timer.SetTimeStamp("Solve heat diffusion");
        LB.AcountForAndLimitPhaseTransformation(Phi, RTC.tStep);
        Timer.SetTimeStamp("LB.AcountForAndLimitPhaseTransformation");
        LB.Solve(Phi, Vel, BC);
        Timer.SetTimeStamp("LB.Solve");
        LB.EnforceTotalMassConservation(Phi);
        Timer.SetTimeStamp("LB.EnforceTotalMassConservation");
        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge Phase Fields");
    }
    return 0;
}

void LocalLBM::InitializeSingle()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle       (i,j,k) = 0;
        DensityWetting (i,j,k,{0}) = VaporDensity[0];
        MomentumDensity(i,j,k,{0}) = {0.,0.,0.};
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LocalLBM::InitializeSphere(const double Radius, const std::array<double,3> pos, const double* rho)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        const int ii = i - pos[0];
        const int jj = j - pos[1];
        const int kk = k - pos[2];
        const double rr = sqrt(ii*ii + jj*jj + kk*kk);

        DensityWetting(i,j,k,{0}) = DensityProfile(rr-Radius,0);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

/*void LocalLBM::InitializeY0(const double lbDensity)
{
    //for(long int i = -DensityWetting.BcellsX(); i < 0; i++)
    //for(long int j = -DensityWetting.BcellsY(); j < DensityWetting.sizeY() + DensityWetting.BcellsY(); j++)
    //for(long int k = -DensityWetting.BcellsZ(); k < DensityWetting.sizeZ() + DensityWetting.BcellsZ(); k++)

    for(long int i = -DensityWetting.BcellsX(); i < DensityWetting.sizeX() + DensityWetting.BcellsX(); i++)
    for(long int j = -DensityWetting.BcellsY(); j < 0; j++)
    for(long int k = -DensityWetting.BcellsZ(); k < DensityWetting.sizeZ() + DensityWetting.BcellsZ(); k++)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        Obstacle    (i, j, k) = true;
        lbPopulations (i,j,k)({n}) = EquilibriumDistribution(lbDensity, lbWeights, {0.0,0.0,0.});
        for (int ii = -dNx; ii <= dNx; ii++)
        for (int jj = -dNy; jj <= dNy; jj++)
        for (int kk = -dNz; kk <= dNz; kk++)
        {
            const dVector3 vel {ii*dx/dt,jj*dx/dt,kk*dx/dt};
            DensityWetting (i,j,k)({n}) +=     lbPopulations(i,j,k)({n})(ii,jj,kk)*dRho;
            MomentumDensity(i,j,k)({n}) += vel*lbPopulations(i,j,k)({n})(ii,jj,kk)*dRho;
        }
    }
}*/
