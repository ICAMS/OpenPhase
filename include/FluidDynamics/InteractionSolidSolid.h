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

 *   File created :   2010
 *   Main contributors :   Oleg Shchyglo; Dmitry Medvedev; Raphael Schiedung
 *
 */

#ifndef INTERACTIONSOLIDSOLID_H
#define INTERACTIONSOLIDSOLID_H

#include "PhaseField.h"

namespace openphase
{

class D3Q27;
class FlowSolverLBM;
class Settings;
class Velocities;

enum class SolidSolidInteractionModel                                           ///< Elasticity homogenization modes and materials models
{
    Standard,                                                                   ///< Standard model developed in ICAMS
    Wang,                                                                       ///< Based on Wang 2006 computer modeling and simulation of solid-state sintering
};

class OP_EXPORTS InteractionSolidSolid: public OPObject                                    ///<  Handles the interaction between solid particles in the fluid environment
{
 public:

    InteractionSolidSolid(Settings& locSettings, std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;
    void ReadInput(std::string InputFileName) override;                         ///< Reads input parameters from a file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input data from the input stream

    void Calculate(PhaseField& Phase,
        const BoundaryConditions& BC,
        const std::function<double(int,int,int)>& MassDensity,
        const double dt) const;                                                 ///< Calculates the solid-solid interaction the entire simulation domain

 private:
    GridParameters Grid;                                                        ///< Simulation grid parameters

    int order;
    int cutoff;
    double strength_wang;                                                       ///< solid-solid interaction strength in [Pa] for wang's model
    double strength;                                                            ///< solid-solid interaction for standard model
    size_t Nphases;

    SolidSolidInteractionModel Model;

    std::vector<std::vector<double>> elastic;
    std::vector<AggregateStates> PhaseAggregateStates;                          ///< Aggregate states of all phases

    void CalculateLocal(const int i, const int j, const int k,
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const double dt) const;                                             ///< Calculates the solid-solid interaction at the point (i,j,k)

    void CalculateLocalWang(const int i, const int j, const int k,
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const std::function<double(int,int,int)>& MassDensity) const;       ///< Calculates the solid-solid interaction at the point (i,j,k) for the Wang model

    static bool applicable(const Grain& grain)
    {
        if (!grain.Exist) return false;
        if (grain.State != AggregateStates::Solid)  return false;
        if (grain.Volume <= std::numeric_limits<double>::epsilon()) return false;
        return true;
    }
};

} //namespace openphase
#endif
