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

#ifndef VanDerWaalsGas_H
#define VanDerWaalsGas_H

#include <float.h>
#include <stdexcept>
#include <string>

namespace openphase
{
/// VanDerWaals
/// This class calculated the pressure of a Van der Waals gas which is used by
/// the FlowSolverLBM. Furthermore, it contains the equilibrium solution of the
/// VanDerWaals equation, which are needed as input for the FlowSolverLBM. The
/// values have been obtained by a separate simulation.
/// TODO Add simulation code.
struct VanDerWaalsGas
{
    static inline double ReducedPressure(
            const double Density,
            const double Temperature,
            const double CriticalDensity = 1.0,
            const double CriticalTemperature = 1.0)                             ///< Equation of State (Van der Waals - reduced properties)
    {
#ifdef DEBUG
        if ((Density < DBL_EPSILON) || (Density > DBL_MAX))
            throw std::invalid_argument("Invalid density"
                    + std::to_string(Density));
#endif
        const double density     = Density/CriticalDensity;
        const double temperature = Temperature/CriticalTemperature;
        const double pressure    = 8.0*density*temperature/(3.0-density) -
                                   3.0*density*density;
        return pressure;
    }
    struct EquilibriumValues_t
    {
        double Temperature;
        double VaporDensity;
        double LiquidDensity;
    };
    static EquilibriumValues_t EquilibriumValues(
            const double Temperature,
            const double CriticalTemperature = 1.0,
            const double LaplacePressure = 0.0);                                ///< Calculates the equilibrium densities
};
} //namespace openphase
#endif
