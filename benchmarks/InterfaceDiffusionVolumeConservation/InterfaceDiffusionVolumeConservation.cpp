#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"

namespace op = openphase;

int main(int argc, char **argv)
{
    std::string InputFileName = op::DefaultInputFileName;

#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_FUNNELED, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::BoundaryConditions  BC;
    op::DoubleObstacle      DO;
    op::InterfaceDiffusion  ID;
    op::InterfaceProperties IP;
    op::PhaseField          Phase;
    op::RunTimeControl      RTC;
    op::Settings            OPSettings;
    op::TimeInfo            Timer;

    OPSettings.Initialize();
    OPSettings.ReadInput(InputFileName);

    BC.Initialize(OPSettings);
    BC.ReadInput(InputFileName);
    //ES.Setup_MPI(OPSettings, BC);

    DO   .Initialize(OPSettings);
    ID   .Initialize(OPSettings);
    IP   .Initialize(OPSettings);
    Phase.Initialize(OPSettings);
    RTC  .Initialize(OPSettings);
    Timer.Initialize(OPSettings, "Timer");

    RTC  .ReadInput(InputFileName);
    DO   .ReadInput(InputFileName);
    ID   .ReadInput(InputFileName);
    IP   .ReadInput(InputFileName);
    Phase.ReadInput(InputFileName);

    const unsigned int Nx = OPSettings.Grid.Nx;
    const unsigned int Ny = OPSettings.Grid.Ny;
    const unsigned int Nz = OPSettings.Grid.Nz;

    int index[3] = {};
    // Set initial geometry of the phases
    index[0] = op::Initializations::Single(Phase, 0, BC);
    index[1] = op::Initializations::Sphere(Phase, 1, Ny/4, Nx/2, Ny/2, Nz/2, BC);
    index[2] = op::Initializations::SectionalPlane(Phase,  2,
            {Nx/2.,Ny/2.,Nz/2.}, {0,1,0}, BC);


    const double Volume20 = Phase.FieldsProperties[index[1]].Volume;

    double ElapsedTime = 0.0;
    //----------------------------- The Time Loop ----------------------------//
    for(int tStep = RTC.tStart; tStep <= RTC.nSteps; tStep++)
    {
        // Output to files
        if (!(tStep%RTC.tFileWrite))
        {
            Phase.WriteVTK(OPSettings, tStep);
        }

        // Output to screens
        if(!(tStep%RTC.tScreenWrite))
        {
            const double AvInterfaceEnergie = DO.Energy(Phase,IP);
            const double Volume1 = Phase.FieldsProperties[index[0]].Volume;
            const double Volume2 = Phase.FieldsProperties[index[1]].Volume;
            const double Volume3 = Phase.FieldsProperties[index[2]].Volume;

            double Error = std::abs(Volume2 - Volume20)/Volume20;

            ElapsedTime += RTC.dt;
            op::ConsoleOutput::WriteTimeStep(tStep, RTC.nSteps);
            op::ConsoleOutput::Write("Elapsed time     [s]", ElapsedTime);
            op::ConsoleOutput::Write("Interface energy [J]", AvInterfaceEnergie);
            op::ConsoleOutput::Write("Volume phase 1", Volume1);
            op::ConsoleOutput::Write("Volume phase 2", Volume2);
            op::ConsoleOutput::Write("Volume phase 3", Volume3);
            op::ConsoleOutput::WriteBlankLine();
            op::ConsoleOutput::Write("Relative error", Error);
            op::ConsoleOutput::WriteBlankLine();

            const double MaxError = 50*DBL_EPSILON;
            if ((Error > MaxError) or (tStep == RTC.nSteps))
            {
                if (Error < MaxError)
                {
                    op::ConsoleOutput::WriteBlankLine();
                    op::ConsoleOutput::WriteLine("_");
                    op::ConsoleOutput::WriteBlankLine();
                    op::ConsoleOutput::WriteSimple("Benchmark successfully completed");
                    op::ConsoleOutput::WriteLine("_");
                    op::ConsoleOutput::WriteBlankLine();
                }

                // Write simulation results
                std::ofstream resultsSim ("Results.sim");
                if (resultsSim.is_open())
                {
                    resultsSim << 2 << "\n";
                    resultsSim << "RelativeVolumeError " << Error << " "
                        << MaxError << "\n";
                    resultsSim.flush();
                    resultsSim.close();
                }
                break;
            }
        }


        // Actual calculation of the interface diffusion
        Timer.SetStart();
        IP.Set(Phase, BC);
        Timer.SetTimeStamp("Sigma.Set()");
        ID.CalculatePhaseFieldIncrements(Phase, IP);
        Timer.SetTimeStamp("ID.Calculate()");
        Phase.MergeIncrements(BC, RTC.dt, false);
        Timer.SetTimeStamp("Phase.MergeIncrements()");

        if(!(tStep%RTC.tScreenWrite)) Timer.PrintWallClockSummary();
    }
    return 0;

#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif
}
