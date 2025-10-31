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

 *   File created :   2013
 *   Main contributors :   Oleg Shchyglo, Raphael Schiedung
 *
 */

#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

namespace openphase
{
struct PhysicalConstants
{
    static constexpr double R           = 8.31446261815324;                     ///< Universal gas constant [J/mol/K]
    static constexpr double N_Avogadro  = 6.02214076e23;                        ///< Avogadro constant [1/mol]
    static constexpr double N_A         = 6.02214076e23;                        ///< Avogadro constant [1/mol] (short name)
    static constexpr double k_Boltzmann = 1.3806488e-23;                        ///< Boltzmann constant [Joule/Kelvin]
    static constexpr double k_B         = 1.3806488e-23;                        ///< Boltzmann constant [Joule/Kelvin] (short name)
    static constexpr double epsilon_0   = 8.8541878128e-12;                     ///< Electric vacuum permittivity [kg^(-1) m^(-3) s^4 A^2] (CODATA 2018)
    static constexpr double mu_0        = 1.25663706212e-6;                     ///< Magnetic vacuum permeability [kg m s^(-2) A^(-2)] (CODATA 2018 )
    static constexpr double q_e         = 1.602176634e-19;                      ///< Elementary charge [C] (exact) (CODATA 2018 )
    static constexpr double h           = 6.62607015e-34;                       ///< Planck constant
    static constexpr double Faraday     = 96485.33289;                          ///< Fraday constant [C/mol] (exact) (CODATA 2018)
};
}// namespace openphase
#endif
