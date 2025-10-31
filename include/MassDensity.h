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
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#ifndef MASSDENSITY_H
#define MASSDENSITY_H

#include "Includes.h"

namespace openphase
{
class PhaseField;
class BoundaryConditions;
class Composition;
class Temperature;
class Settings;

class MassDensity : public OPObject                                             ///< Stores the mass densities of phases in each point.
{
 public:
    MassDensity(){};                                                            ///< Empty constructor.
    MassDensity(Settings& locSettings,
            const std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string FileName) override;                        ///< Reads input values from file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input values from file
    void Set(PhaseField& Phase, Composition& Cx, Temperature& Tx);              ///< Sets the density depending on composition and temperature
    void SetInitial(PhaseField& Phase);                                         ///< Sets initial density

    void WriteVTK(const int tStep, const Settings& locSettings,
                  const int precision = 16);                                    ///< Writes density in VTK format
    bool Write(const Settings& locSettings, const int tStep) const override;                                 ///< Writes the raw composition into a file
    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC, const int tStep) override;                                                 ///< Reads the raw composition from a file
    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets the boundary conditions

    void MoveFrame(const int di, const int dj, const int dk,
                   const BoundaryConditions& BC) override;                      ///< Shifts the data on the storage by di, dj and dk in x, y and z directions correspondingly.
    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;                         ///< Remesh the storage while keeping the data
    void PrintPointStatistics(int x, int y, int z);                             ///< Prints to screen density in a given point (x, y, z)

    MassDensity& operator=(const MassDensity& rhs);                             ///< Copy operator for Density class

    size_t Nphases;                                                             ///< Number of thermodynamic phases

    GridParameters Grid;                                                        ///< Simulation grid parameters

    Storage3D<double, 1> Phase;                                                 ///< Phase densities storage
    Storage3D<double, 0> Total;                                                 ///< Total density storage

    Tensor<double, 1> Initial;                                                  ///< Initial density of all phases

 protected:
 private:
};

} // namespace openphase

#endif // DENSITY_H
