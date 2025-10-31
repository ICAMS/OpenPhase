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
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef SYMMETRYVARIANTS_H
#define SYMMETRYVARIANTS_H

#include "Includes.h"

namespace openphase
{

class Settings;

class OP_EXPORTS SymmetryVariants : OPObject                                    ///< Stores transformation matrices for symmetry variants of thermodynamic phases.
{
 public:
    SymmetryVariants(){};
    SymmetryVariants(Settings& locSettings,
                     const std::string InputFileName = DefaultInputFileName);   ///< Constructor, uses Initialize() and ReadInput()
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads the input from the file InputFileName
    void ReadInput(std::stringstream& inp) override;                            ///< Reads the input from the file InputFileName
    dMatrix3x3& operator()(size_t PhaseIndex, size_t VariantIndex)              ///< Bi-directional access operator
    {
        assert(PhaseIndex < TransformationMatrices.size() && "Access beyond storage range");
        assert(VariantIndex < TransformationMatrices[PhaseIndex].size() && "Access beyond storage range");

        return TransformationMatrices[PhaseIndex][VariantIndex];
    }
    const dMatrix3x3& operator()(size_t PhaseIndex, size_t VariantIndex) const  ///< const access operator
    {
        assert(PhaseIndex < TransformationMatrices.size() && "Access beyond storage range");
        assert(VariantIndex < TransformationMatrices[PhaseIndex].size() && "Access beyond storage range");

        return TransformationMatrices[PhaseIndex][VariantIndex];
    }
    size_t Nvariants(const size_t PhaseIndex) const                             ///< Returns number of symmetry variants of a given phase
    {
        assert(PhaseIndex < TransformationMatrices.size() && "Access beyond storage range");

        return TransformationMatrices[PhaseIndex].size();
    }
    size_t Nphases;                                                             ///< Number of thermodynamic phases
    bool set;                                                                   ///< Indicates if symmetry variants were set
 protected:
     std::vector< std::vector < dMatrix3x3 > > TransformationMatrices;          ///< Transformation matrices storage
 private:
};
} // namespace openphase

#endif
