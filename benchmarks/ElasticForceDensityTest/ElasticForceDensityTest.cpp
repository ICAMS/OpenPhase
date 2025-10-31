/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2020
 *
 *   Authors:    Raphael Schiedung
 *
 *   Copyright (c) 2009-2020 Interdisciplinary Centre for Advanced Materials
 *                 Simulation (ICAMS). Ruhr-Universitaet Bochum. Germany
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "ElasticProperties.h"
#include "ElasticitySolverSpectral.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "Magnetism/LinearMagneticSolver.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"

namespace op = openphase;

int main(int argc, char **argv)
{

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
    //feenableexcept(FE_OVERFLOW);   // the result was too large
    //feenableexcept(FE_UNDERFLOW);  // the result was too small

    // Initialize general settings of OpenPhase
    op::Settings OPSettings(op::DefaultInputFileName);

    op::BoundaryConditions       BC    (OPSettings);
    op::DoubleObstacle           DO    (OPSettings);
    op::ElasticProperties        EP    (OPSettings);
    op::ElasticitySolverSpectral ES    (OPSettings, BC);
    op::InterfaceProperties      IP    (OPSettings);
    op::PhaseField               Phase (OPSettings);
    op::RunTimeControl           RTC   (OPSettings);
    op::TimeInfo                 Time  (OPSettings, "Interface stress benchmark");

    // Start timer
    Time.SetStart();

    /************************** Set initial condition**************************/
    // Set initial conditions
    const int Nx = OPSettings.Grid.Nx;

    //double MassDensity = 1000;
    //double Radius = 0.0;
    if (RTC.Restart == true)
    {
        // Read restart input
        Phase.Read(OPSettings,BC, RTC.tStart);
    }
    else
    {
        op::Initializations::Single(Phase, 0, BC);

        // Initialize cylindrical specimen
        const double ForceDensity = 1e10;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.ForceDensity, 0, )
        {
            double Phix = sin(2.0*op::Pi/Nx*i);
            EP.ForceDensity(i,j,k) += op::dVector3({ForceDensity*Phix,0,0});
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    Phase.Finalize(BC);
    Phase.SetBoundaryConditions(BC);
    Phase.WriteVTK (OPSettings, 0);

    // Calculate elastic properties
    IP.Set(Phase, BC);
    EP.SetEffectiveProperties(Phase, IP);
    Time.SetTimeStamp("Set initial condition");

    EP.WriteElasticStrainsVTK (OPSettings, 0);
    EP.WriteStressesVTK       (OPSettings, 0);
    EP.WriteForceDensityVTK   (OPSettings, 0);

    /************************** Numeric calculation **************************/
    // Ensure that the right boundary conditions are set
    Phase.SetBoundaryConditions(BC);
    op::ConsoleOutput::WriteTimeStep(0, RTC.nSteps);

    /************************** Numeric calculation **************************/

    // Solve for mechanical equilibrium
    op::ConsoleOutput::Write("Solve problem...");
    //ES.Solve(EP, BC, RTC.dt, 1);
    ES.Solve(EP, BC, RTC.dt);
    Time.SetTimeStamp("Solve for mechanical equilibrium");

    op::ConsoleOutput::WriteTimeStep(1, RTC.nSteps);
    Time.SetTimeStamp("Check solution");

    /***************************** Check Solution ****************************/
    op::ConsoleOutput::Write("Write output");
    op::ConsoleOutput::WriteLine();

    double L = EP.EffectiveElasticConstants(0,0,0)(0,1);
    double G = EP.EffectiveElasticConstants(0,0,0)(0,1);
    double error = 0;
    std::ofstream file(OPSettings.TextDir+"/ux.csv");
    file << "i,uxx\n";
    for (int i = 0; i < Nx; i++)
    {
        double uxx_simulation = EP.Displacements(i,0,0)[0];
        double uxx_analytic   = -1.0/(L+2.0*G)*(Nx/2.0*op::Pi)*(Nx/2.0*op::Pi)*sin(2.0*op::Pi/Nx*i);
        file << i << "," << uxx_simulation << "," << uxx_analytic << "\n";
        error += std::abs(uxx_simulation-uxx_analytic);
    }
    error /= Nx;
    op::ConsoleOutput::Write("Deviation form analytic solution: ", error);

    // Write Output for automated testing
    std::ofstream resultsSim ("Results.sim");
    if (resultsSim.is_open())
    {
        resultsSim << 1 << "\n";
        resultsSim << "Error "    << error << " " << error << "\n";
        resultsSim.flush();
        resultsSim.close();
    }

    EP.WriteElasticStrainsVTK (OPSettings, 1);
    EP.WriteStressesVTK       (OPSettings, 1);
    EP.WriteForceDensityVTK   (OPSettings, 1);

    Time.SetTimeStamp("Write to disk");
    Time.PrintWallClockSummary();

#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif

    return 0;
}
