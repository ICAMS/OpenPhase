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

#ifndef GRANDPOTENTIALPHASEDENSITY_IMPLICITPARABOLIC
#define GRANDPOTENTIALPHASEDENSITY_IMPLICITPARABOLIC

#include "GrandPotential/GrandPotentialPhaseDensity_Parabolic.h"
#include "NumericalMethods/RootFindingAlgorithms.h"

namespace openphase
{
class Settings;
struct GrandPotentialPhaseDensity_ImplicitParabolic: GrandPotentialPhaseDensity_Parabolic
{
    GrandPotentialPhaseDensity_ImplicitParabolic(size_t PhaseIdxInp): GrandPotentialPhaseDensity_Parabolic(PhaseIdxInp){};

    void Initialize (Settings& locSettings) override;                           ///<  Initializes global settings
    void ReadInput(std::stringstream& InputFile, int moduleLocation) override;  ///<  Reads input parameters from a file

    double PhasePressure (double Temperature, const Tensor<double,1>& ChemicalPotential) const override
    {
        double locPressure = 0.0;
        for(size_t comp = 0; comp < Ncomp; ++comp)
        {
            // Use shorthands to improve readability
            const double& eps  = EPS[comp];
            const double& zeta = Zeta[comp];
            const double& mu   = ChemicalPotential({comp});
            const double& c0   = C0[comp];
            const double  c    = PhaseConcentration(Temperature, mu, comp);
            locPressure += eps*(c - c0)*(c - c0)/2 + zeta/c + c*mu;
        }
        return locPressure;
    }
    double PhaseConcentration  (double Temperature, double ChemicalPotential, size_t comp) const override
    {
        assert(comp < Ncomp);
        assert(precision < 1);
        assert(precision > 0);

        // Use shorthands to improve readability
        const double& eps  = EPS[comp];
        const double& zeta = Zeta[comp];
        const double& mu   = ChemicalPotential;
        const double& c0   = C0[comp];
        double c    = c0;//DBL_EPSILON;

        auto f  = [&eps, &c0, &mu, &zeta](double c) { return eps*(c - c0) - mu  - zeta/c/c; };
        auto df = [&eps, &zeta]          (double c) { return eps + 2*zeta/c/c/c; };
        try
        {
            RootFindingAlgorithms::PositivNewton(f, df, c, precision, MaxIterations);
        }
        catch (std::runtime_error& ecep)
        {
            ConsoleOutput::WriteWarning(ecep.what(),thisclassname,"PhaseConcentration");
        }
        return c;
    }
    double PhaseSusceptibility (double Temperature, double ChemicalPotential, size_t comp) const override
    {
        assert(comp < Ncomp);
        assert(precision < 1);
        assert(precision > 0);
        const double& eps  = EPS[comp];
        const double& zeta = Zeta[comp];
        const double c     = PhaseConcentration(Temperature, ChemicalPotential, comp);
        return c*c*c/(eps*c*c*c + 2*zeta);
    }

    static constexpr auto thisclassname = "GrandPotentialPhaseDensities_ImplicitParabolic";

 protected:

    double precision;
    size_t MaxIterations;
    std::vector<double> Zeta;                                                   ///< Energy parameter to prevent negative concentration
};
}// namespace openphase
#endif
