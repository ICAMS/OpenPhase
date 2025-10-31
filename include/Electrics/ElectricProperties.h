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

 *   File created :   2016
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef ECTRICPROPERTIES_H
#define ECTRICPROPERTIES_H

#include "Includes.h"

namespace openphase
{
class BoundaryConditions;
class Composition;
class Settings;

/// Storage class for electric material properties
class OP_EXPORTS ElectricProperties : public OPObject
{
 public:
    ElectricProperties(void){};                                                 ///< Empty constructor (does nothing)
    ElectricProperties(Settings& locSettings,
         const std::string InputFileName = DefaultInputFileName);               ///< Constructs and initializes class
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;                            ///< Allocates memory
    void ReadInput(const std::string InputFileName) override;                   ///< Reads the initial parameters for the interface diffusion
	void ReadInput(std::stringstream& inp) override;                            ///< Reads the initial parameters for the interface diffusion
    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets boundary conditions for stored fields properties

    // Field storages
    Storage3D<dVector3,0> PolarizationDensity;                                  ///< Polarization density [ C*m^(-2) = A*s*m^(-2) ]
    Storage3D<dVector3,0> ElectricField;                                        ///< Electric field [V*m^(-1) = kg*m*s^(-3)*A^(-1) ]
    Storage3D<dVector3,0> WaveVector;
    Storage3D<double,0>   ElectricPotential;                                    ///< Electric potential [v = kg*m^2⋅s^(-3)⋅A^(-1) ]
    Storage3D<double,0>   ChargeDensity;                                        ///< Charge density [C*m^(-3) = A*s*m^(-3)]
    Storage<double>       MolarCharge;                                          ///< Charge of component per mole [C/mol = A*s/mol]

    dVector3 ElectricDisplacementField (long int i, long int j, long int k) const {
        return ElectricField(i,j,k) * PhysicalConstants::epsilon_0 + PolarizationDensity(i,j,k);
    }

    void CalculateChargeDensity(const Composition& Cx);

    void WriteVTK(const Settings& locSettings, const int tStep, const int precision=16) const;

    size_t Nphases;                                                             ///< Number of phases
    size_t Ncomp;                                                               ///< Number of components
    GridParameters Grid;                                                        ///< Simulation grid parameters

 private:
};
}
#endif
