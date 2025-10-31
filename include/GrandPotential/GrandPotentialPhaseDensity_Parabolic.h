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

 *   File created :   2021
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef GRANDPOTENTIALPHASEDENSITY_PARABOLIC_H
#define GRANDPOTENTIALPHASEDENSITY_PARABOLIC_H

#include "GrandPotential/GrandPotentialPhaseDensity.h"

namespace openphase
{
class Settings;
struct GrandPotentialPhaseDensity_Parabolic: GrandPotentialPhaseDensity
{
    GrandPotentialPhaseDensity_Parabolic(size_t PhaseIdxInp): GrandPotentialPhaseDensity(PhaseIdxInp){};

    virtual void Initialize (Settings& locSettings) override;                   ///<  Initializes global settings
    virtual void ReadInput  (std::stringstream& InputFile, int moduleLocation) override; ///<  Reads input parameters from a file

    double PhasePressure (double Temperature, const Tensor<double,1>& ChemicalPotential) const override
    {
        double locPressure = 0.0;
        for (size_t comp = 0; comp < Ncomp; ++comp)
        {
            const double& mu = ChemicalPotential({comp});
            locPressure += 0.5*mu*mu/EPS[comp]+C0[comp]*mu;
        }
        return locPressure;
    }
    double PhaseConcentration (double Temperature, double ChemicalPotential, size_t comp) const override
    {
        assert(comp < Ncomp);
        return ChemicalPotential/EPS[comp] + C0[comp];
    }
    double PhaseSusceptibility (double Temperature, double ChemicalPotential, size_t comp) const override
    {
        assert(comp < Ncomp);
        return 1.0/EPS[comp];
    }

    static constexpr auto thisclassname = "GrandPotentialPhaseDensities_Parabolic";

 protected:
    std::vector<double> EPS;                                                    ///< Energy coefficient for parabolic energy
    std::vector<double> C0;                                                     ///< Minimum molar concentration of parabolic free energy
};
}// namespace openphase
#endif
