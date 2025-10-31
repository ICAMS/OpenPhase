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

#ifndef GRANDPOTENTIALPHASEDENSITY_H
#define GRANDPOTENTIALPHASEDENSITY_H

#include "Includes.h"
#include "Settings.h"

namespace openphase
{

struct GrandPotentialPhaseDensity                                               ///< Base class of grand potential densities / equations of state
{
    GrandPotentialPhaseDensity(size_t PhaseIdxInp): PhaseIdx(PhaseIdxInp){};
    virtual void Initialize(Settings& locSettings)                              ///< Allocates and initialises memory
    {
        ElementNames = locSettings.ElementNames;
        PhaseNames   = locSettings.PhaseNames;
        Ncomp        = locSettings.Ncomp;
        dx           = locSettings.Grid.dx;
    }
    virtual void ReadInput(std::stringstream& InputFile, int moduleLocation) = 0;///< Read model parameters

    std::function<double(int,int,int)>        PhasePotential;                   ///< Grand potential density
    std::function<double(int,int,int,size_t)> dChemicalPotential;               ///< First derivative of grand potential density with respect to the chemical potential
    std::function<double(int,int,int,size_t)> dChemicalPotential2;              ///< Second derivative of grand potential density with respect to the chemical potential

    virtual double operator() (int i, int j, int k){return PhasePotential(i,j,k);};  ///< Returns grand potential density of phase

    virtual double PhasePressure       (double Temperature, const Tensor<double,1>& ChemicalPotential) const = 0; ///< Negative grand potential density
    virtual double PhaseConcentration  (double Temperature, double ChemicalPotential, size_t comp)     const = 0; ///< First derivative of grand potential density
    virtual double PhaseSusceptibility (double Temperature, double ChemicalPotential, size_t comp)     const = 0; ///< Second derivative of grand potential density

    virtual double PhasePressure       ([[maybe_unused]] double height, double Temperature, const Tensor<double,1>& ChemicalPotential) const {return PhasePressure       (Temperature, ChemicalPotential);}; ///< Negative grand potential density
    virtual double PhaseConcentration  ([[maybe_unused]] double height, double Temperature, double ChemicalPotential, size_t comp)     const {return PhaseConcentration  (Temperature, ChemicalPotential, comp);}; ///< First derivative of grand potential density
    virtual double PhaseSusceptibility ([[maybe_unused]] double height, double Temperature, double ChemicalPotential, size_t comp)     const {return PhaseSusceptibility (Temperature, ChemicalPotential, comp);}; ///< Second derivative of grand potential density

    void Set (const Storage3D<double,0>& Temperature, const Storage3D<double,1>& ChemicalPotential)
    {
        PhasePotential      = [this, &Temperature, &ChemicalPotential](int i,int j, int k             ) { return -PhasePressure       (Temperature(i,j,k), ChemicalPotential(i,j,k));};
        dChemicalPotential  = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseConcentration  (Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
        dChemicalPotential2 = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseSusceptibility (Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
    }

    void Set (const Storage3D<double,0>& Temperature, const Storage3D<double,1>& ChemicalPotential, const int& GravityDirection)
    {
        if (GravityDirection == 0)
        {
            PhasePotential      = [this, &Temperature, &ChemicalPotential](int i,int j, int k             ) { return -PhasePressure       (i*dx, Temperature(i,j,k), ChemicalPotential(i,j,k));};
            dChemicalPotential  = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseConcentration  (i*dx, Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
            dChemicalPotential2 = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseSusceptibility (i*dx, Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
        }
        else if (GravityDirection == 1)
        {
            PhasePotential      = [this, &Temperature, &ChemicalPotential](int i,int j, int k             ) { return -PhasePressure       (j*dx, Temperature(i,j,k), ChemicalPotential(i,j,k));};
            dChemicalPotential  = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseConcentration  (j*dx, Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
            dChemicalPotential2 = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseSusceptibility (j*dx, Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
        }
        else if (GravityDirection == 2)
        {
            PhasePotential      = [this, &Temperature, &ChemicalPotential](int i,int j, int k             ) { return -PhasePressure       (k*dx, Temperature(i,j,k), ChemicalPotential(i,j,k));};
            dChemicalPotential  = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseConcentration  (k*dx, Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
            dChemicalPotential2 = [this, &Temperature, &ChemicalPotential](int i,int j, int k, size_t comp) { return -PhaseSusceptibility (k*dx, Temperature(i,j,k), ChemicalPotential(i,j,k,{comp}),comp);};
        }
    }
    //void Set (const double& Temperature, const Storage3D<double,1>& ChemicalPotential, const ElasticProperties& EP);

    static constexpr auto thisclassname = "GrandPotentialPhaseDensity";
    static constexpr bool GravityImplemented = false;

 protected:

    double dx;
    size_t Ncomp;                                                               ///< Number of chemical components
    size_t PhaseIdx;                                                            ///< Phase index of grand potential density
    std::vector<std::string> ElementNames;                                      ///< Names of corresponding chemical components
    std::vector<std::string> PhaseNames;                                        ///< Names of corresponding thermodynamic phases
};
}// namespace openphase::
#endif
