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

/******************************************************************************
 * RESULTS HAVE BEEN PUBLISHED  IN:
 ******************************************************************************
 *  P., Schiedung, R., Steinbach, I. and Kästner, M., 2021. Benchmark for the
 *  coupled magneto-mechanical boundary value problem in magneto-active
 *  elastomers. Materials, 14(9), p.2380.
 ******************************************************************************
 */

#include "BoundaryConditions.h"
#include "ConsoleOutput.h"
#include "Initializations.h"
#include "Magnetism/LinearMagneticSolver.h"
#include "Mechanics.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"

namespace op = openphase;

int main(int argc, char *argv[])
{
    std::string InputFileName = op::DefaultInputFileName;
    if (argc > 1) InputFileName = argv[1];

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);   // the result was too large
    //feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::Settings                 OPSettings (InputFileName);

    op::BoundaryConditions       BC    (OPSettings);
    op::ElasticProperties        EP    (OPSettings);
    op::ElasticitySolverSpectral ES    (OPSettings, BC);
    op::PhaseField               Phi   (OPSettings);
    op::LinearMagneticSolver     MS    (OPSettings);
    op::RunTimeControl           RTC   (OPSettings);
    op::TimeInfo                 Timer (OPSettings, "Execution Time Statistics");

    std::fstream inpF(InputFileName, std::ios::in);
    if (!inpF)
    {
        std::stringstream message;
        message << "File \"" << InputFileName << "\" could not be opened";
        op::ConsoleOutput::WriteExit(message.str(),"MagnetoactiveElastomerLinear","main");
        op::OP_Exit(EXIT_FAILURE);
    };
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    const double dx = OPSettings.Grid.dx;
    const int  moduleLocation = op::FileInterface::FindModuleLocation(inp, "MagnetoactiveElastomerLinear");
    const double R1           = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string( "dR1"))/dx;
    const double R2           = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string( "dR2"))/dx;
    const double R3           = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string( "dR3"))/dx;
    const double distance     = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string("dDis"))/dx;
    const double H0           = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string("dH0"));
    const double Phi0         = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string("dPhi0"));
    const double DPhi         = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string("dDPhi"));
    const double accuracy     = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string("dAcc"));
    const size_t MinItr       = op::FileInterface::ReadParameterD(inp, moduleLocation, std::string("iMinItr"));

    const size_t Nx = OPSettings.Grid.Nx;
    const size_t Ny = OPSettings.Grid.Ny;
    const size_t Nz = OPSettings.Grid.Nz;

    const double x0 = (Nx>1) ? Nx/2.0 : 0;
    const double y0 = (Ny>1) ? Ny/2.0 : 0;
    const double z0 = (Nz>1) ? Nz/2.0 : 0;

    const double x1 = Nx/2.0 - distance/2;
    const double x2 = Nx/2.0 + distance/2;

    op::ConsoleOutput::Write("Set initial conditions");
    if(RTC.Restart)
    {
        EP.Read  (OPSettings, BC, RTC.tStart);
        Phi.Read (OPSettings, BC, RTC.tStart);
    }
    else
    {
        RTC.tStart = 0;
        op::Initializations::Single(Phi, 0, BC);
        op::Initializations::Cylinder(Phi, 1, R3, Nz+20, 2, x0, y0, z0, BC);
        op::Initializations::Sphere(Phi, 2, R1, x1, y0, z0, BC);
        op::Initializations::Sphere(Phi, 2, R2, x2, y0, z0, BC);
    }
    MS.SetEffectiveSusceptibility(Phi, BC);

    BC.SetX(MS.chi);
    BC.SetY(MS.chi);
    BC.SetZ(MS.chi);

    EP.SetEffectiveProperties(Phi);

    op::ConsoleOutput::Write("Start Calculation");
    std::string FileName = OPSettings.TextDir + "TimeLog.csv";
    std::fstream log(FileName, std::ios::trunc | std::ios::out);
    log << std::scientific << std::setprecision(16);
    for(int tStep = RTC.tStart; tStep <= RTC.nSteps; tStep++)
    {
        // Solve for magnetic field and force density
        MS.H0x = H0*std::cos(op::Pi*(Phi0+tStep*DPhi)/180.0);
        MS.H0y = H0*std::sin(op::Pi*(Phi0+tStep*DPhi)/180.0);
        MS.Solve(BC,0.0001);
        MS.CalcForceDensity(EP, BC);

        // Solve for elastic deformation
        size_t Iteration = 0;
        double DeltaDistance    = 0.0;
        double DeltaDistanceOld = 0.0;
        const int j0 = std::round(y0);
        const int k0 = std::round(z0);
        std::stringstream FName;
        FName << OPSettings.TextDir << "DeltaDistance" << tStep << ".csv";
        std::fstream logDist(FName.str(), std::ios::trunc | std::ios::out);
        auto EndSolve = [&tStep, &OPSettings, &EP, &j0, &k0, &x1, &x2, &dx, &Iteration, &DeltaDistance, &DeltaDistanceOld, &logDist, &MinItr, &accuracy]()
        {
            DeltaDistance = 0.0;
            for (int i = std::round(x1); i < std::round(x2); i++)
            {
                DeltaDistance += EP.ElasticStrains(i,j0,k0)[0]*dx;
            }
            double diff = std::abs((DeltaDistanceOld - DeltaDistance)/DeltaDistance);
            std::cout << "Elastic Solver Interation: "<<  Iteration << "," << "\tParticle Displacement: " << DeltaDistance << "," << "\t Relative Particle Displacement Change: " << diff << "\n";
            logDist   << Iteration << "," << DeltaDistance << "," << diff << "\n";
            //EP.WriteElasticStrainsVTK (tStep, OPSettings);

            DeltaDistanceOld = DeltaDistance;
            Iteration++;
            if (diff < accuracy and Iteration > MinItr) return true;
            else return false;
        };
        ES.Solve(EP, BC, RTC.dt, EndSolve);

        if (RTC.WriteToScreen())
        {
            op::dVector3 H0 = MS.H0(0,0,0); // Resulting external H-field
            op::dVector3 H  = MS.H(0,0,0);  // Resulting external H-field
            op::dVector3 B  = MS.B(0,0,0);  // Resulting external B-field
            op::dVector3 M  = MS.M(std::round(x1), j0, k0); // Resulting particle magnetisation

            DeltaDistance = 0.0;
            for (int i = std::round(x1); i < std::round(x2); i++)
            {
                DeltaDistance += EP.ElasticStrains(i,j0,k0)[0]*dx;
            }

            std::array<std::stringstream,2> line;
            line[1] << std::scientific << std::setprecision(16);
            op::ConsoleOutput::WriteTimeStep(tStep, RTC.nSteps);
            op::ConsoleOutput::WriteWithLog(line, tStep, "Phi", Phi0+tStep*DPhi);
            op::ConsoleOutput::WriteWithLog(line, tStep, "Applied external H0", H0.length());
            op::ConsoleOutput::WriteWithLog(line, tStep, "Resulting external H",  H.length());
            op::ConsoleOutput::WriteWithLog(line, tStep, "Resulting external B",  B.length());
            op::ConsoleOutput::WriteWithLog(line, tStep, "Magnetisation of particle", M.length());
            op::ConsoleOutput::WriteWithLog(line, tStep, "DeltaDistance", DeltaDistance);
            op::ConsoleOutput::WriteWithLog(line, tStep, "Elastic Energy", EP.Energy(Phi));
            op::ConsoleOutput::WriteLineToLogfile(log, line, tStep);
        }

        if (RTC.WriteVTK())
        {
            EP.WriteElasticStrainsVTK (OPSettings, tStep);
            //EP.WriteStressesVTK       (tStep, OPSettings);
            //EP.WriteForceDensityVTK   (tStep, OPSettings);
            //ES.WriteVTK               (tStep, OPSettings);
            MS.WriteVTK               (tStep, OPSettings);
            //Phi.WriteDistortedVTK     (tStep, OPSettings, ES, 1000000);
            Phi.WriteVTK              (OPSettings, tStep);
        }

        if (RTC.WriteRawData())
        {
            EP.Write(OPSettings, tStep);
            Phi.Write(OPSettings, tStep);
        }
    } //end time loop
    return 0;
}
