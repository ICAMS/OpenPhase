/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 */

#ifndef CHEMICALDIFFUSION_H
#define CHEMICALDIFFUSION_H

#include "Includes.h"

namespace openphase
{

class Settings;
class PhaseField;
class Composition;
class Temperature;
class BoundaryConditions;
class DiffusionProperties;
class ThermodynamicPropertiesEQP;

class ChemicalDiffusion_Impl;

class OP_EXPORTS ChemicalDiffusion: public OPObject                             ///< Chemical diffusion solver
{
 public:
    ~ChemicalDiffusion();
    ChemicalDiffusion();
    ChemicalDiffusion(Settings& locSettings,
                     const std::string InputFileName = DefaultInputFileName);   ///< Initializes storages, sets internal variables.

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from the input stream


    void Solve(PhaseField& Phi, Composition& Cx, Temperature& Tx,
               DiffusionProperties& DP, 
               BoundaryConditions& BC, double dt);                              ///< Solves diffusion equation
    void SolveSF(PhaseField& Phi, Composition& Cx, Temperature& Tx,
               DiffusionProperties& DP, 
               BoundaryConditions& BC, double dt);                              ///< Solves diffusion equation
    void Solve(PhaseField& Phi, Composition& Cx, Temperature& Tx,
               DiffusionProperties& DP, ThermodynamicPropertiesEQP& TP,
               BoundaryConditions& BC, double dt);                              ///< Solves diffusion equation

    void PrintStatistics(const Settings& locSettings, const int tStep);         ///< Prints solver execution statistics
 protected:
 private:
    std::unique_ptr<ChemicalDiffusion_Impl> impl_;

};

} // namespace openphase
#endif
