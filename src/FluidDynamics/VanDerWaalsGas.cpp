/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *     
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *   File created :   2018
 *   Main contributors :   Raphael Schiedung
 *
 */

#include <vector>
#include "FluidDynamics/VanDerWaalsGas.h"

namespace openphase
{
VanDerWaalsGas::EquilibriumValues_t VanDerWaalsGas::EquilibriumValues(
        const double Temperature,
        const double CriticalTemperature,
        const double LaplacePressure)
{
    /// NOTE! The following reduced equilibrium values, Temperature,
    /// VaporDensity, and LiquidDensity, of the Van der Waals equation
    /// should not be changed! They may be updated by more accurate values
    /// or an analytical solution of the Van der Waals equation. They are
    /// used as input parameters of the FlowSolverLBM.
    /// TODO Add surface tension
    const std::vector<VanDerWaalsGas::EquilibriumValues_t> EquilibriumValues =
    {{0.01,3.73587E-144,2.99109E+000},
     {0.02,2.79725E-071,2.98212E+000},
     {0.03,4.89974E-047,2.97310E+000},
     {0.04,5.93062E-035,2.96403E+000},
     {0.05,1.00176E-027,2.95489E+000},
     {0.06,6.36164E-023,2.94570E+000},
     {0.07,1.67599E-019,2.93645E+000},
     {0.08,6.04048E-017,2.92714E+000},
     {0.09,5.80403E-015,2.91777E+000},
     {0.10,2.20916E-013,2.90834E+000},
     {0.11,4.29885E-012,2.89884E+000},
     {0.12,5.06427E-011,2.88928E+000},
     {0.13,4.04833E-010,2.87965E+000},
     {0.14,2.39271E-009,2.86996E+000},
     {0.15,1.09051E-008,2.86041E+000},
     {0.16,4.14938E-008,2.85063E+000},
     {0.17,1.34771E-007,2.84010E+000},
     {0.18,3.81679E-007,2.83046E+000},
     {0.19,9.70874E-007,2.82008E+000},
     {0.20,2.22717E-006,2.81057E+000},
     {0.21,4.71698E-006,2.80034E+000},
     {0.22,9.34579E-006,2.78940E+000},
     {0.23,1.73310E-005,2.77932E+000},
     {0.24,3.04878E-005,2.76855E+000},
     {0.25,5.12821E-005,2.75862E+000},
     {0.26,8.26446E-005,2.74801E+000},
     {0.27,1.28205E-004,2.73673E+000},
     {0.28,1.92678E-004,2.72628E+000},
     {0.29,2.80899E-004,2.71518E+000},
     {0.30,3.98406E-004,2.70416E+000},
     {0.31,5.52486E-004,2.69324E+000},
     {0.32,7.51880E-004,2.68168E+000},
     {0.33,1.00100E-003,2.67023E+000},
     {0.34,1.31062E-003,2.65887E+000},
     {0.35,1.68634E-003,2.64760E+000},
     {0.36,2.14133E-003,2.63574E+000},
     {0.37,2.68097E-003,2.62398E+000},
     {0.38,3.31126E-003,2.61233E+000},
     {0.39,4.04858E-003,2.60010E+000},
     {0.40,4.90196E-003,2.58799E+000},
     {0.41,5.88235E-003,2.57533E+000},
     {0.42,6.99301E-003,2.56345E+000},
     {0.43,8.26446E-003,2.55037E+000},
     {0.44,9.61538E-003,2.53807E+000},
     {0.45,1.12108E-002,2.52525E+000},
     {0.46,1.29534E-002,2.51193E+000},
     {0.47,1.48588E-002,2.49875E+000},
     {0.48,1.69492E-002,2.48571E+000},
     {0.49,1.92308E-002,2.47219E+000},
     {0.50,2.17391E-002,2.45821E+000},
     {0.51,2.44499E-002,2.44439E+000},
     {0.52,2.73973E-002,2.43072E+000},
     {0.53,3.05810E-002,2.41663E+000},
     {0.54,3.38983E-002,2.40211E+000},
     {0.55,3.75940E-002,2.38777E+000},
     {0.56,4.14938E-002,2.37304E+000},
     {0.57,4.56621E-002,2.35793E+000},
     {0.58,5.00000E-002,2.34247E+000},
     {0.59,5.46448E-002,2.32721E+000},
     {0.60,5.98802E-002,2.31160E+000},
     {0.61,6.49351E-002,2.29568E+000},
     {0.62,7.04225E-002,2.27946E+000},
     {0.63,7.63359E-002,2.26296E+000},
     {0.64,8.26446E-002,2.24669E+000},
     {0.65,8.92857E-002,2.22965E+000},
     {0.66,9.61538E-002,2.21239E+000},
     {0.67,1.03734E-001,2.19491E+000},
     {0.68,1.11483E-001,2.17723E+000},
     {0.69,1.19474E-001,2.15889E+000},
     {0.70,1.28041E-001,2.14041E+000},
     {0.71,1.36986E-001,2.12179E+000},
     {0.72,1.46199E-001,2.10261E+000},
     {0.73,1.56006E-001,2.08290E+000},
     {0.74,1.66389E-001,2.06271E+000},
     {0.75,1.77305E-001,2.04248E+000},
     {0.76,1.88679E-001,2.02143E+000},
     {0.77,2.00401E-001,2.00000E+000},
     {0.78,2.12766E-001,1.97824E+000},
     {0.79,2.25734E-001,1.95580E+000},
     {0.80,2.39808E-001,1.93274E+000},
     {0.81,2.53807E-001,1.90913E+000},
     {0.82,2.69542E-001,1.88466E+000},
     {0.83,2.84900E-001,1.85977E+000},
     {0.84,3.02115E-001,1.83385E+000},
     {0.85,3.19489E-001,1.80701E+000},
     {0.86,3.37838E-001,1.77936E+000},
     {0.87,3.58423E-001,1.75070E+000},
     {0.88,3.78788E-001,1.72087E+000},
     {0.89,4.01606E-001,1.68976E+000},
     {0.90,4.25532E-001,1.65728E+000},
     {0.91,4.50450E-001,1.62311E+000},
     {0.92,4.78469E-001,1.58680E+000},
     {0.93,5.10204E-001,1.54823E+000},
     {0.94,5.43478E-001,1.50670E+000},
     {0.95,5.78035E-001,1.46177E+000},
     {0.96,6.21118E-001,1.41203E+000},
     {0.97,6.66667E-001,1.35575E+000},
     {0.98,7.24638E-001,1.28949E+000},
     {0.99,8.06452E-001,1.20351E+000},
     {1.00,1.00000E+000,1.00000E+000}};

    VanDerWaalsGas::EquilibriumValues_t PreviousValue = EquilibriumValues[0];
    VanDerWaalsGas::EquilibriumValues_t Interpolated  = {0.0,0.0,0.0};
    for (auto& Value: EquilibriumValues)
    {
        if (Value.Temperature == Temperature)
        {
            return Value;
        }
        else if (Temperature < Value.Temperature and
                 Temperature > PreviousValue.Temperature)
        {
            Interpolated = PreviousValue;

            const double dT = (Value.Temperature-PreviousValue.Temperature);
            const double DT = (Temperature-PreviousValue.Temperature);
            const double A  = DT/dT;

            Interpolated.Temperature   = Temperature;
            Interpolated.VaporDensity  += (Value.VaporDensity - PreviousValue.VaporDensity)*A;
            Interpolated.LiquidDensity += (Value.LiquidDensity - PreviousValue.LiquidDensity)*A;

            return Interpolated;
        }
        PreviousValue = Value;
    }
    return Interpolated;
}
} //namespace openphase
