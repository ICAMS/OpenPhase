#include "BoundaryConditions.h"
#include "Includes.h"
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include <fstream>
#include <iostream>

namespace op = openphase;

void CalculateContactAngels(
        double (& angle)[3],
        op::PhaseField& Phase);

double DoDiagnostics(
        op::PhaseField& Phase,
        op::RunTimeControl& RTC,
        int tStep,
        double MaxError);

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
    op::InterfaceDiffusion  ID;
    op::InterfaceProperties IP;
    op::PhaseField          Phase;
    op::RunTimeControl      RTC;
    op::Settings            OPSettings;

    OPSettings.Initialize();
    OPSettings.ReadInput(InputFileName);

    BC.Initialize(OPSettings);
    BC.ReadInput(InputFileName);
    //ES.Setup_MPI(OPSettings, BC);

    RTC  .Initialize(OPSettings);
    ID   .Initialize(OPSettings);
    IP   .Initialize(OPSettings);
    Phase.Initialize(OPSettings);

    RTC  .ReadInput(InputFileName);
    ID   .ReadInput(InputFileName);
    IP   .ReadInput(InputFileName);
    Phase.ReadInput(InputFileName);


    // Set initial geometry of the phases
    op::Initializations::Young3(Phase, 0, 0, 0, 0, BC);

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
            double MaxError = 5.;
            double Error    = DoDiagnostics(Phase, RTC, tStep, MaxError);
            // Check if the correct contact angel is reached
            if (Error < MaxError)
            {
                op::ConsoleOutput::WriteBlankLine();
                op::ConsoleOutput::WriteLine("_");
                op::ConsoleOutput::WriteBlankLine();
                op::ConsoleOutput::WriteSimple("Benchmark successfully completed");
                op::ConsoleOutput::WriteLine("_");
                op::ConsoleOutput::WriteBlankLine();
                break;
            }
        }

        // Calculation of the interface diffusion
        IP.Set(Phase, BC);
        ID.CalculatePhaseFieldIncrements(Phase, IP);
        Phase.MergeIncrements(BC, RTC.dt, false);
    }
    return 0;

#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif
}

double DoDiagnostics(
        op::PhaseField& Phase,
        op::RunTimeControl& RTC,
        int tStep,
        double MaxError)
{
    // Calculate contact angels between the phases
    double angle[3] = {};
    CalculateContactAngels(angle, Phase);

    // Calculate deviation from equilibrium
    double Error = 0.;
    for (int i = 0; i < 3; i++)
    {
        Error += std::abs(angle[i]-120);
    }
    Error /= 3.0;

    // Write Results to file
    std::ofstream log("TextData/angles.csv", std::ios_base::app | std::ios_base::out);
    if (log.is_open())
    {
        log  << tStep       << ","
             << Error       << ","
             << angle[0]    << ","
             << angle[1]    << ","
             << angle[2]    << "\n";
    }

    op::ConsoleOutput::WriteTimeStep(tStep, RTC.nSteps);
    op::ConsoleOutput::WriteBlankLine();
    op::ConsoleOutput::Write("Theta 1", angle[0]);
    op::ConsoleOutput::Write("Theta 2", angle[1]);
    op::ConsoleOutput::Write("Theta 3", angle[2]);
    op::ConsoleOutput::WriteBlankLine();
    op::ConsoleOutput::Write("Error", Error);

    // Write simulation results
    if ((Error < MaxError) or (tStep = RTC.nSteps))
    {
        std::ofstream resultsSim ("Results.sim");
        if (resultsSim.is_open())
        {
            resultsSim << 1 << "\n";
            resultsSim << "Error " << Error << " " << MaxError;
            resultsSim.flush();
            resultsSim.close();
        }
    }
    return Error;
}

op::dVector3 InterfaceDirection(const op::PhaseField& Phase, const int Tx,
        const int Tz, const int indexA,
        const int indexB, const int BoxSize)
{
    // Search for point on interface-line AB
    double MinDeviation = 8.;
    double Px = 0.;
    double Pz = 0.;

    // Define a box the triple point
    const int eps = BoxSize;
    const int iMin = ((Tx - eps) > 0)        ? std::abs(Tx - eps) : 0;
    const int kMin = ((Tz - eps) > 0)        ? std::abs(Tz - eps) : 0;
    const int iMax = ((Tx + eps) < Phase.Grid.Nx) ? std::abs(Tx + eps) : Phase.Grid.Nx;
    const int kMax = ((Tz + eps) < Phase.Grid.Nz) ? std::abs(Tz + eps) : Phase.Grid.Nz;

    // Search only on the surface of the box for a point!
    for (int i = iMin; i <= iMax; i++)
    for (int k = kMin; k <= kMax; k++)
    if (((i == iMax) or (i == iMin)) or ((k == kMax) or (k == kMin)))
    {
        {
            // Check if point is on interface-line AB
            double locDeviation = double();
            locDeviation += std::pow(Phase.Fields(i,0,k).get_value(indexA) - 0.5, 2);
            locDeviation += std::pow(Phase.Fields(i,0,k).get_value(indexB) - 0.5, 2);
            locDeviation = std::sqrt(locDeviation);

            // If the point is nearer on the interface line save it
            if (locDeviation < MinDeviation)
            {
                Px = i;
                Pz = k;
                MinDeviation = locDeviation;
            }
        }
    }

    // Calculate direction vector
    op::dVector3 direction = op::dVector3();
    direction.setX(Px-Tx);
    direction.setY(0.0);
    direction.setZ(Pz-Tz);
    direction.normalize();
    return direction;
}

void CalculateContactAngels(double (& angle)[3],
        op::PhaseField& Phase)
{
    // Search for triple point
    double Tx = 0.0;
    double Tz = 0.0;
    double MinDeviation = 1.0;

    const int Nx = Phase.Grid.Nx;
    const int Nz = Phase.Grid.Nz;

    // Search for triple point
    for (int i = 0; i <= Nx/2+Nx/8; i++)
    for (int k = 0; k <= Nz/2;      k++)
    if(Phase.Fields(i,0,k).interface())
    {
        double locDeviation = 0.0;
        locDeviation += std::pow(Phase.Fields(i,0,k).get_value(0) - 1./3, 2);
        locDeviation += std::pow(Phase.Fields(i,0,k).get_value(1) - 1./3, 2);
        locDeviation += std::pow(Phase.Fields(i,0,k).get_value(2) - 1./3, 2);
        locDeviation = std::sqrt(locDeviation);

        if (locDeviation < MinDeviation)
        {
            MinDeviation = locDeviation;
            Tx = i;
            Tz = k;
        }
    }

    // Calculate orientation of interface-line in polar coordinates
    const int eps = Phase.Grid.Nx/4 - Phase.Grid.iWidth/2;
    op::dVector3 O01 = InterfaceDirection(Phase, Tx, Tz, 0, 1, eps);
    op::dVector3 O02 = InterfaceDirection(Phase, Tx, Tz, 0, 2, eps);
    op::dVector3 O12 = InterfaceDirection(Phase, Tx, Tz, 1, 2, eps);

    // Calculate contact angles at the triple point
    const double Rad2Grad = 360./(2.*op::Pi);
    angle[0] = Rad2Grad * acos(O01*O02);
    angle[1] = Rad2Grad * acos(O01*O12);
    angle[2] = Rad2Grad * acos(O02*O12);
}
