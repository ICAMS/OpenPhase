/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2016
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
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "ElasticProperties.h"
#include "ElasticitySolverSpectral.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"

namespace op = openphase;

int main()
{
    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::Settings OPSettings(op::DefaultInputFileName);

    op::BoundaryConditions       BC    (OPSettings);
    op::DoubleObstacle           DO    (OPSettings);
    op::ElasticProperties        EP    (OPSettings);
    op::ElasticitySolverSpectral ES    (OPSettings);
    op::InterfaceDiffusion       ID    (OPSettings);
    op::InterfaceProperties      IP    (OPSettings);
    op::PhaseField               Phase (OPSettings);
    op::RunTimeControl           RTC   (OPSettings);


    if (RTC.Restart == true)
    {
        Phase.Read(OPSettings, BC, RTC.tStart);
    }
    else
    {
        const double R  = (OPSettings.Grid.Nx+1)/5.0;
        const double x0 = (OPSettings.Grid.Nx+1)/2.0;
        const double y0 = (OPSettings.Grid.Ny+1)/2.0;
        const double z0 = (OPSettings.Grid.Nz+1)/2.0;

        op::Initializations::Single(Phase, 0, BC);
        op::Initializations::Sphere(Phase, 1, R, x0-R, y0, z0, BC);
        op::Initializations::Sphere(Phase, 1, R, x0+R, y0, z0, BC);
    }
    Phase.SetBoundaryConditions(BC);

    // Initialize elastic properties
    op::ConsoleOutput::WriteSimple("Starting simulation...");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteToScreen())
        {
            //Solve for Mechanical Equilibrium
            EP.SetEffectiveProperties(Phase, IP);
            ES.Solve(EP, BC, RTC.dt);

            Phase.WriteVTK            (OPSettings, RTC.tStep);
            EP.WriteStressesVTK       (OPSettings, RTC.tStep);
            EP.WriteElasticStrainsVTK (OPSettings, RTC.tStep);
            EP.WriteTotalStrainsVTK   (OPSettings, RTC.tStep);
        }

        if(RTC.WriteRawData())
        {
            Phase.Write(OPSettings, RTC.tStep);
        }

        if(RTC.WriteVTK())
        {
            op::ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps);
            op::ConsoleOutput::Write("Simulation Time  [s]", RTC.SimulationTime);
            op::ConsoleOutput::Write("Interface Energy [J]", DO.Energy(Phase,IP));
        }

        IP.Set(Phase, BC);
        ID.CalculatePhaseFieldIncrements(Phase, IP);
        Phase.MergeIncrements(BC, RTC.dt, false);

    }
    return 0;
}
