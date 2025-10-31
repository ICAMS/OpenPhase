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
 *   Main contributors :   Efim Borukhovich; Philipp Engels; Oleg Shchyglo
 *
 */

#ifndef ORIENTATIONS_H
#define ORIENTATIONS_H

#include "Includes.h"
#include "PhaseField.h"

namespace openphase
{

class Settings;
class BoundaryConditions;
class PhaseField;
class Quaternion;
class Crystallography;
class Velocities;
class AdvectionHR;

class OP_EXPORTS Orientations : public OPObject                                 ///< Module which stores and handles orientation data like rotation matrices, Euler angles etc.
{
 public:

    Orientations(){};
    Orientations(Settings& locSettings);

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes the module, allocate the storage, assign internal variables
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remesh and reallocate orientations

    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets boundary conditions
    void WriteVTK(const Settings& locSettings,
                  const int tStep,
                  const int precision=16);                                      ///< Writes rotations to the VTK file
    void WriteTotalVTK(const Settings& locSettings,
                       const PhaseField& Phase,
                       const int tStep,
                       const int precision = 16) const;                         ///< Writes total rotations (including grains orientations) to the VTK file
    bool Write(const Settings& locSettings, const int tStep) const override;    ///< Writes orientations storage to binary file
    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC, const int tStep) override;          ///< Reads orientations from binary file
    void WriteGrainEBSDDataQuaternions(const Settings& locSettings,
                                       const PhaseField& Phase,
                                       const int tStep);                        ///< Writes EBSD file (for MTex). In interface, majority phase is used.
    void WriteMisorientationsVTK(const Settings& locSettings,
                                 const int tStep,
                                 const int precision,
                                 const std::string measure = "deg") const;      ///< Writes misorientations in VTK format
    void WriteRotated100VectorVTK(const Settings& locSettings,
                                  const int tStep,
                                  const int precision = 16);

    void PrintPointStatistics(int x, int y, int z);
    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
                const BoundaryConditions& BC,
                const double dt, const double tStep) override;                  ///< Advects orientations

    GridParameters Grid;                                                        ///< Simulation grid parameters

    int Nphases;

    Storage3D<Quaternion, 0> Quaternions;
    Storage3D<Quaternion, 0> QuaternionsDot;

    Orientations& operator= (const Orientations& rhs);

    dMatrix3x3 getTotalRotation(const PhaseField& Phase,
                                    const int i, const int j, const int k) const///< Returns the local rotation including grain orientation and deformation induced rotation
    {
        return Quaternions(i,j,k).RotationMatrix*Phase.EffectiveOrientation(i,j,k);
    }

    dMatrix3x3 getTotalGrainRotation(const PhaseField& Phase, int i, int j, int k,
                                                                      int alpha) const ///< Returns the total rotation seen by the specified grain including deformation induced rotation
    {
        return Quaternions(i,j,k).RotationMatrix*Phase.FieldsProperties[alpha].Orientation.RotationMatrix;
    }

 protected:
 private:
};
}// namespace openphase
#endif
