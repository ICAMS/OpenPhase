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

#ifndef GRANDPOTENTIALPHASEDENSITY_COMPRESSIBLEFLUID_H
#define GRANDPOTENTIALPHASEDENSITY_COMPRESSIBLEFLUID_H

#include "GrandPotential/GrandPotentialPhaseDensity.h"
#include "GrandPotential/GrandPotentialPhaseDensity_CompressibleFluid.h"

namespace openphase
{
class Settings;
struct GrandPotentialPhaseDensity_CompressibleFluid: GrandPotentialPhaseDensity ///< Equation of state for compressible fluid or solid
{
    GrandPotentialPhaseDensity_CompressibleFluid(size_t PhaseIdxInp): GrandPotentialPhaseDensity(PhaseIdxInp){};

    void Initialize(Settings& locSettings); ///<  Initializes global settings
    void ReadInput(std::stringstream& InputFile, int moduleLocation) override;  ///<  Reads input parameters from a file

    double PhasePressure ([[maybe_unused]] double Temperature, const Tensor<double,1>& ChemicalPotential) const override
    {
        const double& K   = BulkModulus[0];
        const double& Vm0 = ReferenceMolarVolume[0];
        const double& mu  = ChemicalPotential({0});
        const double& mu0 = ReferenceChemicalPotential[0];
        const double& p0  = ReferencePressure[0];

        assert(K*Vm0 + mu0 - mu > 0 && "Invalid input parameters!");

        return p0 - K*std::log((K*Vm0 + mu0 - mu)/(K*Vm0));
    }
    double PhaseConcentration   ([[maybe_unused]] double Temperature, double ChemicalPotential, size_t comp) const override
    {
        const double& K   = BulkModulus[0];
        const double& Vm0 = ReferenceMolarVolume[0];
        const double& mu  = ChemicalPotential;
        const double& mu0 = ReferenceChemicalPotential[0];

        assert(comp == 0);
        assert(K*Vm0 + mu0 - mu > 0 && "Invalid input parameters!");

        return K/(K*Vm0 + mu0 - mu);
    }
    double PhaseSusceptibility ([[maybe_unused]] double Temperature, double ChemicalPotential, size_t comp) const override
    {
        const double& K   = BulkModulus[0];
        const double& Vm0 = ReferenceMolarVolume[0];
        const double& mu  = ChemicalPotential;
        const double& mu0 = ReferenceChemicalPotential[0];

        assert(comp == 0);
        assert(K*Vm0 + mu0 - mu > 0 && "Invalid input parameters!");

        return K/(K*Vm0 + mu0 - mu)/(K*Vm0 + mu0 - mu);
    }

    double PhasePressure (double height, [[maybe_unused]] double Temperature, const Tensor<double,1>& ChemicalPotential) const override
    {
        const double& K   = BulkModulus[0];
        const double& M   = MolarMass[0];
        const double& Vm0 = ReferenceMolarVolume[0];
        const double& g   = GravityAcceleration;
        const double& h   = height;
        const double& mu  = ChemicalPotential({0});
        const double& mu0 = ReferenceChemicalPotential[0];
        const double& p0  = ReferencePressure[0];

        assert(M*g*h + K*Vm0 + mu0 - mu > 0 && "Invalid input parameters!");

        return p0 - K*std::log((M*g*h + K*Vm0 + mu0 - mu)/(K*Vm0));
    }
    double PhaseConcentration (double height, [[maybe_unused]] double Temperature, double ChemicalPotential, size_t comp) const override
    {
        const double& K   = BulkModulus[0];
        const double& M   = MolarMass[0];
        const double& Vm0 = ReferenceMolarVolume[0];
        const double& g   = GravityAcceleration;
        const double& h   = height;
        const double& mu  = ChemicalPotential;
        const double& mu0 = ReferenceChemicalPotential[0];

        assert(comp == 0);
        assert(M*g*h + K*Vm0 + mu0 - mu > 0 && "Invalid input parameters!");

        return K/(M*g*h + K*Vm0 + mu0 - mu);
    }
    double PhaseSusceptibility (double height, [[maybe_unused]] double Temperature, double ChemicalPotential, size_t comp) const override
    {
        const double& K   = BulkModulus[0];
        const double& M   = MolarMass[0];
        const double& Vm0 = ReferenceMolarVolume[0];
        const double& g   = GravityAcceleration;
        const double& h   = height;
        const double& mu  = ChemicalPotential;
        const double& mu0 = ReferenceChemicalPotential[0];

        assert(comp == 0);
        assert(M*g*h + K*Vm0 + mu0 - mu > 0 && "Invalid input parameters!");

        return K/(M*g*h + K*Vm0 + mu0 - mu)/(M*g*h + K*Vm0 + mu0 - mu);
    }

    static constexpr auto thisclassname = "GrandPotential_PhaseDensities_CompressibleFluid";
    static constexpr bool GravityImplemented = true;

 protected:
    std::vector<double> BulkModulus;                                            ///< Bulk modulus of species in phase
    std::vector<double> ReferenceChemicalPotential;                             ///< Chemical potential of component at the reference pressure
    std::vector<double> ReferenceMolarVolume;                                   ///< Molar volume of species in phase at the reference pressure
    std::vector<double> ReferencePressure;                                      ///< Reference pressure of component
    std::vector<double> MolarMass;                                              ///< Molar mass of component
    double GravityAcceleration;                                                 ///< Gravitational acceleration [m/s^2]
};
}// namespace openphase
#endif
