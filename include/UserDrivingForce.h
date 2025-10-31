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
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef USERDRIVINGFORCE_H
#define USERDRIVINGFORCE_H

#include "Includes.h"
namespace openphase
{

class Settings;
class PhaseField;
class DrivingForce;
class Temperature;
class Composition;

enum class UserDrivingForceModes                                                ///< User driving force modes of operation for each phase pair
{
    None,                                                                       ///< User driving force is inactive
    Value,                                                                      ///< User driving force using constant value
    Formula                                                                     ///< User driving force using constant value and/or temperature dependent formula
};

class OP_EXPORTS UserDrivingForce : public OPObject                             ///< User driving force class
{
 public:
    UserDrivingForce(void){};
    UserDrivingForce(Settings& locSettings,
                     const std::string InputFileName = DefaultInputFileName);

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;
    void ReadInput(std::string InputFileName) override;
    void ReadInput(std::stringstream& inp) override;
    void SetDrivingForce(PhaseField& Phase, DrivingForce& dGab);                ///< Sets the user provided driving force to a specified value
    void SetDrivingForce(PhaseField& Phase, DrivingForce& dGab, Temperature& Tx);///< Sets the user provided driving force considering temperature effect
    void SetDrivingForce(PhaseField& Phase, DrivingForce& dGab, Temperature& Tx, Composition& Cx);///< Sets the user provided driving force considering temperature and composition effect

    size_t Nphases;

    static void SetDrivingForce(PhaseField& Phase, DrivingForce& dGab,
                                int indexA, int indexB, double dGvalue);        ///< Sets specified driving force for a pair of phase fields
 protected:
 private:
    Matrix<UserDrivingForceModes> Mode;
    Matrix<double> Value;
    Matrix<double> Teq;                                                         ///< Equilibrium temperature between pairs of phases
    Matrix<double> Slope;                                                       ///< Slope of the Teq composition dependence
    Matrix<size_t> Component;                                                   ///< Chemical component index for Teq composition dependence
    Matrix<double> LatentHeat;                                                  ///< Latent heat values for the driving force formula
};

}// namespace openphase
#endif
