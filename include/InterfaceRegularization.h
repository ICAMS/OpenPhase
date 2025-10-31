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

 *   File created :   2023
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#ifndef INTERFACECURVATURE_H
#define INTERFACECURVATURE_H

#include "Includes.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"

namespace openphase
{

class Settings;
class PhaseField;
class BoundaryConditions;
class InterfaceProperties;

/***************************************************************/
class OP_EXPORTS InterfaceRegularization : public OPObject                      ///< The interface regularization module. Provides the storage and manipulation methods.
{
 public:
    InterfaceRegularization() : Averaging(false) {};
    InterfaceRegularization(Settings& locSettings,
                 const std::string InputFileName = DefaultInputFileName);       ///< Initializes the storage and internal variables of the class.
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes the storage and internal variables of the class.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads settings
    void ReadInput(std::stringstream& inp) override;                            ///< Reads settings
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remeshes the storage while keeping the data
    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets the boundary conditions

    void Clear(void);                                                           ///< Deletes curvature in the storage.

    void Average(const PhaseField& Phase, const BoundaryConditions& BC);        ///< Averages curvature across the interface

    void MergePhaseFieldIncrements(PhaseField& Phase, InterfaceProperties& IP); ///< Merges the curvature contribution into the phase field increments.

    void PrintPointStatistics(const int i, const int j, const int k) const;

    void WriteVTK(const Settings& locSettings, const int tStep,
                  const size_t indexA, const size_t indexB,
                  const int precision = 16) const;                              ///< Writes the curvature of the interface between phase-fields with indexA and indexB

    NodeDF at(const double x, const double y, const double z) const;            ///< Arbitrary point access operator for curvature. Uses tri-linear interpolation

    ClearingModes ClearingMode;                                                 ///< Clearing mode for the curvature storage (facilitates VTK output)

    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    int Range;                                                                  ///< Radius of a sphere around a grid point over which the curvature is averaged
    double PhiThreshold;                                                        ///< Outlines the inner part of the interface for curvature averaging.

    Storage3D<NodeDF, 0> Curvature;                                             ///< Curvature storage

    InterfaceRegularization& operator= (const InterfaceRegularization& rhs);    ///< Copy operator for Curvature class

    bool Averaging;                                                             ///< Control parameter for averaging the curvature over the interface.

 protected:
 private:
    void SetWeightsSR(const PhaseField& Phase);                                 ///< Sets averaging weights for curvature across the interface in single resolution
    void CollectAverageSR(const PhaseField& Phase);                             ///< Average curvature across the interface in single resolution
    void DistributeAverageSR(const PhaseField& Phase);                          ///< Average curvature across the interface in single resolution

    void SetWeightsDR(const PhaseField& Phase);                                 ///< Sets averaging weights for curvature across the interface in double resolution
    void CollectAverageDR(const PhaseField& Phase);                             ///< Average curvature across the interface in double resolution
    void DistributeAverageDR(const PhaseField& Phase);                          ///< Average curvature across the interface in double resolution

    void MergePhaseFieldIncrementsSR(PhaseField& Phase, InterfaceProperties& IP);///< Merges the curvature contribution into the phase field increments.
    void MergePhaseFieldIncrementsDR(PhaseField& Phase, InterfaceProperties& IP);///< Merges the curvature contribution into the phase field increments.
};
} // namespace openphase
#endif
