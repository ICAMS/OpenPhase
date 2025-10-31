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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Johannes Goerler; Raphael Schiedung; Stephan Hubig
 *
 */

#ifndef DOUBLEOBSTACLE_H
#define DOUBLEOBSTACLE_H

#include "Includes.h"

namespace openphase
{

class BoundaryConditions;
class InterfaceProperties;
class PhaseField;
class Settings;
class ElasticProperties;
class DrivingForce;
class InterfaceRegularization;

/*******************************************************************************/
class PotentialCorrections
{
 public:
    double alpha;
    double beta;
    double scal_prod;
    double phi_tl;
    double phi_tu;

    double stencil_weight;

    int d_x;
    int d_y;
    int d_z;
};

class OP_EXPORTS DoubleObstacle : public OPObject
{
 public:
    DoubleObstacle(){};
    DoubleObstacle(Settings& locSettings,
            std::string InputFileName = DefaultInputFileName)                   ///< Initializes module with constructor.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes module.

    void CalculateCurvatureDrivingForce(PhaseField& Phase,
                                        InterfaceProperties& IP,
                                        DrivingForce& dG);                      ///< Calculates generalized curvature contribution to the driving force.

    void CalculatePhaseFieldIncrementsSharp(PhaseField& Phase,
                                            InterfaceProperties& IP);           ///< Calculates interface curvature related driving force using sharp interface formalism of A.Finel et.al

    void CalculatePhaseFieldIncrements(PhaseField& Phase,
                                       InterfaceProperties& IP);                ///< Calculates interface curvature related driving force. Includes triple junction energy.
    void CalculatePhaseFieldIncrements(PhaseField& Phase,
                                       InterfaceProperties& IP,
                                       DrivingForce& dG);                       ///< Calculates interface curvature related driving force. Includes triple junction energy and curvature subtraction

    void CalculatePhaseFieldIncrements(PhaseField& Phase,
                                       InterfaceProperties& IP,
                                       InterfaceRegularization& Kappa);              ///< Calculates interface curvature related driving force. Includes triple junction energy.
    void CalculatePhaseFieldIncrements(PhaseField& Phase,
                                       InterfaceProperties& IP,
                                       DrivingForce& dG,
                                       InterfaceRegularization& Kappa);              ///< Calculates interface curvature related driving force. Includes triple junction energy.

    void StabilizeThinChannels(PhaseField& Phase,
                               InterfaceProperties& IP,
                               double penalty);                                 ///< Stabilizes thin channels between two interfaces with the same phase-field index on both sides of the channel

    double Energy(
            const PhaseField& Phase,
            const InterfaceProperties& IP) const;                               ///< Provides total interface energy in the simulation domain.
    double Energy(
            const PhaseField& Phase,
            const InterfaceProperties& IP,
            const ElasticProperties& EP) const;                                 ///< Provides total interface energy in the simulation domain with accounting for the elastic volume change.
    double AverageEnergyDensity(
            const PhaseField& Phase,
            const InterfaceProperties& IP) const;                               ///< Provides the average interface energy density in the simulation domain.
    double PointEnergy(
            const PhaseField& Phase,
            const InterfaceProperties& IP,
            const int i, const int j, const int k) const;                       ///< Provides the interface energy in a given point (i,j,k).
    double PointEnergy(
            const PhaseField& Phase,
            const InterfaceProperties& IP,
            const ElasticProperties& EP,
            const int i, const int j, const int k) const;                       ///< Provides the interface energy in a given point (i,j,k) with accounting for the elastic volume change.
    void WriteEnergyVTK(const int tStep,
            const Settings& locSettings,
            const PhaseField& Phase,
            const InterfaceProperties& IP) const;                               ///< Writes interface energy in VTK format for a given time step tStep.
    void WriteEnergyVTK(const int tStep,
            const Settings& locSettings,
            const PhaseField& Phase,
            const InterfaceProperties& IP,
            const ElasticProperties& EP) const;                                 ///< Writes interface energy in VTK format for a given time step tStep with accounting for the elastic volume change.

 protected:
 private:
    void CalculatePhaseFieldIncrementsSharpSR(PhaseField& Phase,
                                              InterfaceProperties& IP);         ///< Calculates interface curvature related driving force in single resolution using sharp interface formalism of A.Finel et.al
    void CalculatePhaseFieldIncrementsSharpDR(PhaseField& Phase,
                                              InterfaceProperties& IP);         ///< Calculates interface curvature related driving force in double resolution using sharp interface formalism of A.Finel et.al

    void CalculatePhaseFieldIncrementsSR(PhaseField& Phase,
                                         InterfaceProperties& IP);              ///< Calculates interface curvature related driving force. Includes triple junction energy.
    void CalculatePhaseFieldIncrementsDR(PhaseField& Phase,
                                         InterfaceProperties& IP);              ///< Calculates interface curvature related driving force. Includes triple junction energy.

    void CalculatePhaseFieldIncrementsSR(PhaseField& Phase,
                                         InterfaceProperties& IP,
                                         InterfaceRegularization& Kappa);       ///< Calculates interface curvature related driving force. Includes triple junction energy.
    void CalculatePhaseFieldIncrementsDR(PhaseField& Phase,
                                         InterfaceProperties& IP,
                                         InterfaceRegularization& Kappa);       ///< Calculates interface curvature related driving force. Includes triple junction energy.

    void CalculateFullAnisotropy(const PhaseField& Phase,
                                 InterfaceProperties& IP,
                                 DrivingForce& dG);                             ///< Calculates driving force resulting form anisotropic interface energy.
    };
}
#endif
