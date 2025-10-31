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

 *   File created :   2017
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef MAGNETICPROPERTIES_H
#define MAGNETICPROPERTIES_H

#include "OPObject.h"

namespace openphase
{
class BoundaryConditions;
class Orientations;
class PhaseField;
class Settings;

/// Storage class for magnetic material properties
class MagneticProperties : public OPObject
{
 public:
    MagneticProperties(void){};                                                 ///< Constructor (does nothing)
    MagneticProperties(Settings& locSettings);                                  ///< Constructor call Initialize and ReadInput

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;                            ///< Initializes storage class
    void ReadInput(const std::string InputFileName) override;                   ///< Reads the initial parameters for the interface diffusion
    void ReadInput(std::stringstream& inp) override;                            ///< Reads elastic properties from the input file
    //double MagneticEnergy(void) const;                                          ///< Calculates and returns the magnetic energy
    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets boundary conditions for stored fields properties
    void SetGrainsProperties(const PhaseField& Phase);                          ///< Sets elastic properties for each grain
    void SetGrainsProperties(const PhaseField& Phase, const Orientations& OR);  ///< Sets elastic properties for each grain according to its orientation
    void WriteMagneticFieldVTK(const int tStep, const Settings& locSettings) const; ///< Writes magnetic field in vtk format
    void WriteMagnetisationVTK(const int tStep, const Settings& locSettings) const; ///< Writes magnetisation field in vtk format

    // Material property storage
    Storage<dMatrix3x3> Permeability;                                           ///< Permeability for each phase
    Storage<dMatrix3x3> PermeabilityPhase;                                      ///< Permeability for each phase, orientation is not considered.
    Storage<dMatrix3x3> RelativePermeability;                                   ///< Relative Permeability for each phase
    Storage<dMatrix3x3> RelativePermeabilityPhase;                              ///< Relative Permeability for each phase, orientation is not considered.

    // Field storages
    Storage3D<dVector3,1> Magnetisation;                                        ///< Magnetisation field of the material
    Storage3D<dVector3,0> MagneticField;                                        ///< Magnetic field of the material

    // System size discretization properties
    unsigned int Sublattices;                                                   ///< Number magnetisation sub lattices`
    unsigned int Nphases;                                                       ///< Number of phases

    GridParameters Grid;                                                        ///< Simulation grid parameters

 private:
};
}
#endif
