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

 *   File created :   2019
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef BENZIGAS_H
#define BENZIGAS_H

namespace openphase
{
struct BenziGas
{
    struct EquilibriumValues_t
    {
        double BenziParameter;
        double lbVaporDensity;
        double lbLiquidDensity;
        double lbSurfaceTension;
    };
    static EquilibriumValues_t EquilibriumValues(const double BennziParameter); ///< Calculates the equilibrium densities
};
} //namespace openphase
#endif
