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

#ifndef DRIVINGFORCE_H
#define DRIVINGFORCE_H

#include "Includes.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"

namespace openphase
{

class Settings;
class PhaseField;
class BoundaryConditions;
class InterfaceProperties;
class SemiImplicitCoupling;
class Temperature;
class ThermodynamicPropertiesEQP;
class ThermodynamicPropertiesEQP;

enum class AveragingWeightsModes
{
    Range,
    PhaseFields,
    Counter
};

/***************************************************************/
class OP_EXPORTS DrivingForce : public OPObject                                 ///< The driving force module. Provides the storage and manipulation methods.
{
    friend ThermodynamicPropertiesEQP;
    friend SemiImplicitCoupling;

 public:
    DrivingForce() {};
    DrivingForce(Settings& locSettings,
                 const std::string InputFileName = DefaultInputFileName);       ///< Initializes the storage and internal variables of the driving force class.
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes the storage and internal variables of the driving force class.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads driving force settings
    void ReadInput(std::stringstream& inp) override;                            ///< Reads driving force settings
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remeshes the storage while keeping the data

    void Clear(void);                                                           ///< Deletes driving forces in the storage. Needs to be called at the end/beginning of each time step!

    void Average(const PhaseField& Phase, const BoundaryConditions& BC);        ///< Averages driving forces across the interface

    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets the boundary conditions

    void MergePhaseFieldIncrements(PhaseField& Phase, InterfaceProperties& IP); ///< Merges the driving force into the phase field increments.
    NodeDF CalcUnified(const PhaseField& Phase, const BoundaryConditions& BC);  ///< Reduces raw driving force to a single number
    void   Unify(const PhaseField& Phase, const BoundaryConditions& BC);        ///< Reduces raw driving force to a single number

    double MaxTimeStep(PhaseField& Phase, InterfaceProperties& IP,
                       Settings& OPSettings, double TheorLimit = 1.0,
                       double NumLimit = 1.0e-3);                               ///< Calculates minimum number of phase field iterations for a given time step

    void PrintDiagnostics(void);                                                ///< Prints driving force statistics to screen. Indicates if driving force is overshooting.
    void PrintPointStatistics(const int i, const int j, const int k) const;
  
    void WriteVTK(const Settings& locSettings, const int tStep,
                  const size_t indexA, const size_t indexB,
                  const int precision = 16) const;                              ///< Writes the driving force acting from phase field with indexB onto phase field with indexA
    void WriteVTK(const Settings& locSettings,
                  const PhaseField& Phi,
                  const int tStep,
                  const int precision = 16) const;                              ///< Writes the average driving force acting between all thermodynamic phases, not between individual grains

    NodeDF Force_at(const double x, const double y, const double z) const;      ///< Arbitrary point access operator for driving force. Uses tri-linear interpolation

    double GetDrivingForce(PhaseField& Phi, const int i, const int j, const int k, const size_t alpha, const size_t beta) const;

    Storage3D<NodeDF, 0> Force;                                                 ///< Driving force storage

    DrivingForce& operator= (const DrivingForce& rhs);                          ///< Copy operator for DrivingForce class

    void AverageGlobal(PhaseField& Phase, double time);
    NodeDF AverageDG;
    NodeDF AverageDG1;
    NodeDF AverageDG2;

    AveragingWeightsModes WeightsMode;

 protected:
    Matrix<double> Limit;                                                       ///< Driving force limit for pairs of phases

 private:

    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    int Range;                                                                  ///< Radius of a sphere around a grid point over which the driving force is averaged
    double PhiThreshold;                                                        ///< Outlines the inner part of the interface for driving force averaging.

    size_t OvershootCounter;                                                    ///< Number of driving force overshooting events
    Matrix<double> MAXOvershootPOS;                                             ///< Maximum positive driving force overshoot for each phase pair
    Matrix<double> MAXOvershootNEG;                                             ///< Maximum negative driving force overshoot for each phase pair
    Matrix<double> MAXDrivingForcePOS;                                          ///< Maximum positive driving force value for each phase pair
    Matrix<double> MAXDrivingForceNEG;                                          ///< Maximum negative driving force value for each phase pair

    bool Averaging;                                                             ///< Control parameter for averaging the driving force over the interface.
    bool Unifying;                                                              ///< Control parameter for unification of the driving force over the interface.
    bool Limiting;                                                              ///< If set to true enables the driving force limiting to stabilize interface profile

    double maxPsi;

    void SkipAverage(const PhaseField& Phase);                                  ///< Sets local average driving force to its raw value
    void SetWeights(const PhaseField& Phase);                                   ///< Sets local averaging weights
    void CollectAverage(const PhaseField& Phase);                               ///< First part of Average
    void DistributeAverage(const PhaseField& Phase);                            ///< Second part of Average

    void MergePhaseFieldIncrementsSR(PhaseField& Phase, InterfaceProperties& IP);///< Merges the driving force into the phase field increments.
    void MergePhaseFieldIncrementsDR(PhaseField& Phase, InterfaceProperties& IP);///< Merges the driving force into the phase field increments.

};
} // namespace openphase
#endif
