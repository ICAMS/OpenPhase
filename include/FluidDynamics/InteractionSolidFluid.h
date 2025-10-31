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

 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung
 *
 */

#ifndef INTERACTIONSOLIDFLUID_H
#define INTERACTIONSOLIDFLUID_H

#include "Globals.h"
#include "GrainsProperties.h"

namespace openphase
{

class BoundaryConditions;
class PhaseField;
class Settings;
class Velocities;

struct OP_EXPORTS InteractionSolidFluid: public OPObject                                   ///< Processes solid liquid interaction
{
    InteractionSolidFluid(Settings& locSettings, std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;
    void ReadInput(std::string InputFileName) override;                         ///< Reads input parameters from a file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input data from the input stream

    static dVector3 LocalGrainVelocity(long i, long j, long k, const Grain& grain,
            const BoundaryConditions& BC, const PhaseField& Phase);             ///< Returns solid velocity of grain at (i,j,k)

    void CalculateSolidVelocities(PhaseField& Phase, Velocities& Vel,
            const BoundaryConditions& BC, const double dt) const;               ///< Calculates velocities of solid particles

 private:
    static void CalculateCenterOfMass(PhaseField& Phase,
            const BoundaryConditions& BC);                                      ///< Calculates the center of mass for each particle (contains OMP and MPI reduction)
    static void CalculateMomentOfInterita(PhaseField& Phase,
            const BoundaryConditions& BC);                                      ///< Calculates the moment of inertia for each particle (contains OMP and MPI reduction)

    // Helper functions
    static void EnforceZeroTotalMomentum(PhaseField& Phase);                    ///< Enforces zero total momentum of rigid bodies
    static void EnforceZeroTotalAngularMomentum(PhaseField& Phase);             ///< Enforces zero total angular momentum of rigid bodies
    static void LimitVelocity(PhaseField& Phase, double VLimit);                ///< Enforces velocity limit of rigid bodies

    void CalculateRigidBodyVelocities(PhaseField& Phase,
            const BoundaryConditions& BC, const double dt) const;               ///< Calculates rigid body velocities (dot(v) = m*F)
    void CalculateRigidBodyVelocitiesWithoutInertia(PhaseField& Phase,
            const BoundaryConditions& BC, const double dt) const;               ///< Calculates rigid body velocities without inertia (v = m*F*tau)
    static double LocalSolidFraction(long i, long j, long k,
            const PhaseField& Phase);                                           ///< Returns local solid fraction
    static void CalculateSolidPhaseVelocities(Velocities& Vel, const PhaseField& Phase,
            const BoundaryConditions& BC);                                      ///< Calculates velocities resulting form rigid body motion (contains OMP and MPI reduction)

    static bool applicable(const Grain& grain)                                  ///< Returns false if grain can not be moved
    {
        if (!grain.Exist) return false;
        if (grain.State != AggregateStates::Solid)  return false;
        if (!grain.Mobile) return false;
        if (grain.Volume <= 1) return false;
        return true;
    }

    double TauTranslation;                                                      ///< Model parameter if no inertial of solid bodies is considered
    double TauRotation;                                                         ///< Model parameter if no inertial of solid bodies is considered
    double VelocityLimit;                                                       ///< Limit of solid body velocity
    bool DoInertia;                                                             ///< If true the inertia of solid bodies is calculated
    bool DoRotations;                                                           ///< If true the rotation of solid bodies is calculated
    bool DoLimitVelocity;                                                       ///< If true the velocity is limited
    bool DoEnforceZeroTotalMomentum;                                            ///< If true the law of momentum conservation is enforced
    bool DoEnforceZeroTotalAngularMomentum;                                     ///< If true the law of angular momentum conservation is enforced
};

} //namespace openphase
#endif
