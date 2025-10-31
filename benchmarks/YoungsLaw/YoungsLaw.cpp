/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2022 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2022 OpenPhase Solutions GmbH,
 *                 Universitaetsstrasse 136, D-44799 Bochum, Germany.
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
 */

#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "Initializations.h"
#include "BoundaryConditions.h"

using namespace std;
using namespace openphase;

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                        OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                  RTC(OPSettings);
    PhaseField                      Phi(OPSettings);
    DoubleObstacle                  DO(OPSettings);
    InterfaceProperties             IP(OPSettings);
    BoundaryConditions              BC(OPSettings);

    Initializations::Fractional(Phi, 0, 1, 25, BC);
    Initializations::Sphere(Phi, 2, 23, 50, 0, 50, BC);
    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        IP.Set(Phi,BC);

        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.KeepPhaseFieldsVolume();
        Phi.MergeIncrements(BC, RTC.dt);
        if (RTC.WriteVTK())
        {
           //  Output to file
           Phi.WriteVTK(OPSettings, RTC.tStep);
        }

        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", to_string(I_En));
            ConsoleOutput::WriteTimeStep(RTC, message);
        }
    } //end of time loop
    return 0;
}
