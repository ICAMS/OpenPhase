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

 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef EQUILIBRIUMDATA_H
#define EQUILIBRIUMDATA_H

#include <cassert>

namespace openphase
{

class EquilibriumData                                                          ///< Structure for storing pair equilibrium data
{
 public:
    double eq_temperature;                                                      ///< Equilibrium temperature at which the data is taken
    double eq_entropy;                                                          ///< Entropy of the phase in equilibrium
    double eq_fraction;                                                         ///< Volume fraction of the phase in equilibrium
    double heat_capacity;                                                       ///< Molar heat capacity of the phase

    Tensor<double,1> part_coefficient;                                          ///< Partitioning coefficient
    Tensor<double,1> nom_composition;                                           ///< Nominal composition for the local equilibrium
    Tensor<double,1> eq_composition;                                            ///< Equilibrium composition of the phase
    Tensor<double,1> eq_slope;                                                  ///< Slope of the solvus line of the phase

    bool valid;                                                                 ///< True if the dataset has valid data
    bool out_of_range;                                                          ///< True if dataset is out of temperature or composition range

    double eq_composition_extra(size_t comp, double locTx) const                ///< Equilibrium composition of a given component extrapolated to local temperature
    {
        double eq_composition_return = eq_composition({comp});
        if(eq_slope({comp}) != 0.0)
        {
            eq_composition_return += (locTx - eq_temperature)/eq_slope({comp});
        }
        return eq_composition_return;
    }

    ~EquilibriumData()
    {

    }
    EquilibriumData():
        eq_temperature(0.0),
        eq_entropy(0.0),
        eq_fraction(0.0),
        heat_capacity(0.0),
        part_coefficient(),
        nom_composition(),
        eq_composition(),
        eq_slope(),
        valid(false),
        out_of_range(true)
    {

    }
    EquilibriumData(size_t Ncomp):
        eq_temperature(0.0),
        eq_entropy(0.0),
        eq_fraction(0.0),
        heat_capacity(0.0),
        part_coefficient({Ncomp}),
        nom_composition({Ncomp}),
        eq_composition({Ncomp}),
        eq_slope({Ncomp}),
        valid(false),
        out_of_range(true)
    {

    }
    EquilibriumData(const EquilibriumData& rhs):
        eq_temperature(rhs.eq_temperature),
        eq_entropy(rhs.eq_entropy),
        eq_fraction(rhs.eq_fraction),
        heat_capacity(rhs.heat_capacity),
        part_coefficient(rhs.part_coefficient),
        nom_composition(rhs.nom_composition),
        eq_composition(rhs.eq_composition),
        eq_slope(rhs.eq_slope),
        valid(rhs.valid),
        out_of_range(rhs.out_of_range)
    {

    }
    EquilibriumData(EquilibriumData&& rhs):
        eq_temperature(rhs.eq_temperature),
        eq_entropy(rhs.eq_entropy),
        eq_fraction(rhs.eq_fraction),
        heat_capacity(rhs.heat_capacity),
        part_coefficient(rhs.part_coefficient),
        nom_composition(rhs.nom_composition),
        eq_composition(rhs.eq_composition),
        eq_slope(rhs.eq_slope),
        valid(rhs.valid),
        out_of_range(rhs.out_of_range)
    {

    }
    EquilibriumData& operator=(const EquilibriumData& rhs)
    {
        eq_temperature   = rhs.eq_temperature;
        eq_entropy       = rhs.eq_entropy;
        eq_fraction      = rhs.eq_fraction;
        heat_capacity    = rhs.heat_capacity;
        part_coefficient = rhs.part_coefficient;
        nom_composition  = rhs.nom_composition;
        eq_composition   = rhs.eq_composition;
        eq_slope         = rhs.eq_slope;
        valid            = rhs.valid;
        out_of_range     = rhs.out_of_range;

        return *this;
    }
    EquilibriumData& operator=(EquilibriumData&& rhs)
    {
        eq_temperature   = rhs.eq_temperature;
        eq_entropy       = rhs.eq_entropy;
        eq_fraction      = rhs.eq_fraction;
        heat_capacity    = rhs.heat_capacity;
        part_coefficient = rhs.part_coefficient;
        nom_composition  = rhs.nom_composition;
        eq_composition   = rhs.eq_composition;
        eq_slope         = rhs.eq_slope;
        valid            = rhs.valid;
        out_of_range     = rhs.out_of_range;

        return *this;
    }
    void Allocate(size_t n_comp)
    {
        part_coefficient.Allocate({n_comp});
        nom_composition.Allocate({n_comp});
        eq_composition.Allocate({n_comp});
        eq_slope.Allocate({n_comp});
    }
};

class PairEquilibriumData                                                      ///< Structure for storing pair equilibrium data
{
 public:
    double eq_entropy_A;                                                        ///< Entropy of the first phase in equilibrium
    double eq_entropy_B;                                                        ///< Entropy of the second phase in equilibrium

    double eq_fraction_A;                                                       ///< Volume fraction of the first phase in equilibrium
    double eq_fraction_B;                                                       ///< Volume fraction of the second phase in equilibrium

    Tensor<double,1> eq_composition_A;                                          ///< Equilibrium composition of the first phase
    Tensor<double,1> eq_composition_B;                                          ///< Equilibrium composition of the second phase

    bool error;                                                                 ///< Indicates if the data is invalid (true) or valid (false)

    ~PairEquilibriumData()
    {

    }
    PairEquilibriumData():
        eq_entropy_A(0.0),
        eq_entropy_B(0.0),
        eq_fraction_A(0.0),
        eq_fraction_B(0.0),
        eq_composition_A(),
        eq_composition_B(),
        error(false)
    {

    }
    PairEquilibriumData(size_t Ncomp):
        eq_entropy_A(0.0),
        eq_entropy_B(0.0),
        eq_fraction_A(0.0),
        eq_fraction_B(0.0),
        eq_composition_A({Ncomp}),
        eq_composition_B({Ncomp}),
        error(false)
    {

    }
    PairEquilibriumData(const PairEquilibriumData& rhs):
        eq_entropy_A(rhs.eq_entropy_A),
        eq_entropy_B(rhs.eq_entropy_B),
        eq_fraction_A(rhs.eq_fraction_A),
        eq_fraction_B(rhs.eq_fraction_B),
        eq_composition_A(rhs.eq_composition_A),
        eq_composition_B(rhs.eq_composition_B),
        error(false)
    {

    }
    PairEquilibriumData(PairEquilibriumData&& rhs):
        eq_entropy_A(rhs.eq_entropy_A),
        eq_entropy_B(rhs.eq_entropy_B),
        eq_fraction_A(rhs.eq_fraction_A),
        eq_fraction_B(rhs.eq_fraction_B),
        eq_composition_A(rhs.eq_composition_A),
        eq_composition_B(rhs.eq_composition_B),
        error(false)
    {

    }
    PairEquilibriumData& operator=(const PairEquilibriumData& rhs)
    {
        eq_entropy_A     = rhs.eq_entropy_A;
        eq_entropy_B     = rhs.eq_entropy_B;
        eq_fraction_A    = rhs.eq_fraction_A;
        eq_fraction_B    = rhs.eq_fraction_B;
        eq_composition_A = rhs.eq_composition_A;
        eq_composition_B = rhs.eq_composition_B;
        error            = rhs.error;

        return *this;
    }
    PairEquilibriumData& operator=(PairEquilibriumData&& rhs)
    {
        eq_entropy_A     = rhs.eq_entropy_A;
        eq_entropy_B     = rhs.eq_entropy_B;
        eq_fraction_A    = rhs.eq_fraction_A;
        eq_fraction_B    = rhs.eq_fraction_B;
        eq_composition_A = rhs.eq_composition_A;
        eq_composition_B = rhs.eq_composition_B;
        error            = rhs.error;

        return *this;
    }

    bool dual_phase() const                                                     ///< Indicates if both phases are stable at the given conditions
    {
        bool temp = false;
        if((eq_fraction_A > DBL_EPSILON) and (eq_fraction_B > DBL_EPSILON))
        {
            temp = true;
        }
        return temp;
    }
};

class SingleEquilibriumData
{
 public:
    double HeatCapacity;                                                        ///< Molar heat capacity
    Tensor<double,2> Dij;                                                       ///< Diffusion coefficients matrix

    SingleEquilibriumData()
    {
        HeatCapacity = 60.0;
    };
    SingleEquilibriumData(size_t Ncomp)
    {
        HeatCapacity = 60.0;
        Dij.Allocate({Ncomp,Ncomp});
        Dij.set_to_zero();
    }
    SingleEquilibriumData(const SingleEquilibriumData& old)
    {
        HeatCapacity = old.HeatCapacity;
        Dij = old.Dij;
    }
    SingleEquilibriumData& operator=(const SingleEquilibriumData& rhs)
    {
        HeatCapacity = rhs.HeatCapacity;
        Dij = rhs.Dij;
        return *this;
    }
};

class SingleThermodynamicData
{
 public:
    double GibbsEnergy;
    Tensor<double,1> PartialDerivative;
    Tensor<double,2> SecondPartialDerivative;

    SingleThermodynamicData()
    {
        GibbsEnergy = 0.0;
    };
    SingleThermodynamicData(size_t Ncons)
    {
        GibbsEnergy = 0.0;
        PartialDerivative.Allocate({Ncons});
        SecondPartialDerivative.Allocate({Ncons,Ncons});
        PartialDerivative.set_to_zero();
        SecondPartialDerivative.set_to_zero();
    }
    SingleThermodynamicData(const SingleThermodynamicData& old)
    {
        GibbsEnergy = old.GibbsEnergy;
        PartialDerivative = old.PartialDerivative;
        SecondPartialDerivative = old.SecondPartialDerivative;
    }
    SingleThermodynamicData& operator=(const SingleThermodynamicData& rhs)
    {
        GibbsEnergy = rhs.GibbsEnergy;
        PartialDerivative = rhs.PartialDerivative;
        SecondPartialDerivative = rhs.SecondPartialDerivative;
        return *this;
    }
};

class SingleKineticData
{
 public:
    Tensor<double,1> AtomicMobility;

    SingleKineticData()
    {

    };
    SingleKineticData(size_t Ncons)
    {
        AtomicMobility.Allocate({Ncons});
        AtomicMobility.set_to_zero();
    }
    SingleKineticData(const SingleKineticData& old)
    {
        AtomicMobility = old.AtomicMobility;
    }
    SingleKineticData& operator=(const SingleKineticData& rhs)
    {
        AtomicMobility = rhs.AtomicMobility;
        return *this;
    }
};

class SingleChemicalPotentialData
{
 public:
    Tensor<double,1> ChemicalPotential;
    double Gibbs;

    SingleChemicalPotentialData()
    {
        Gibbs = 0.0;
    };
    SingleChemicalPotentialData(size_t Ncomp)
    {
        Gibbs = 0.0;
        ChemicalPotential.Allocate({Ncomp});
        ChemicalPotential.set_to_zero();
    }
    SingleChemicalPotentialData(const SingleChemicalPotentialData& old)
    {
        Gibbs = old.Gibbs;
        ChemicalPotential = old.ChemicalPotential;
    }
    SingleChemicalPotentialData& operator=(const SingleChemicalPotentialData& rhs)
    {
        Gibbs = rhs.Gibbs;
        ChemicalPotential = rhs.ChemicalPotential;
        return *this;
    }
};

}// namespace openphase
#endif
