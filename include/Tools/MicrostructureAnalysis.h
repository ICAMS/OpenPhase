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

 *   File created :   2018
 *   Main contributors :   Oleg Shchyglo; Hesham Salama
 *
 */

#ifndef MICROSTRUCTUREANALYSIS_H
#define MICROSTRUCTUREANALYSIS_H

#include "Includes.h"
#include "VTK.h"
#include "Crystallography.h"
#include "Orientations.h"

namespace openphase
{
class PhaseField;
class SymmetryVariants;
class Crystallography;

class OP_EXPORTS MicrostructureAnalysis
{
 public:
    static void WriteEBSDDataQuaternions(const PhaseField& Phase, const int tStep, double scale = 1.0);
    static void WriteEBSDDataQuaternionsSlice(const PhaseField& Phase, const int tStep, const char Axis, const int Position, double scale = 1.0);

    static void WriteEBSDDataQuaternions(const PhaseField& Phase, const SymmetryVariants& SV, const int tStep, double scale = 1.0);
    static void WriteEBSDDataQuaternionsSlice(const PhaseField& Phase, const SymmetryVariants& SV, const int tStep, const char Axis, const int Position, double scale = 1.0);
    static void WriteEBSDVTK(PhaseField& Phase, Crystallography& CR, const Settings& locSettings, EulerConvention locConvention, const int sd, const int tStep);
    static dVector3 IPFColor(size_t sd, EulerAngles& tempEuler, Crystallography& CR);

    ///< Writes RGB as Vector Field
    template <class T>
    static void WriteVTKRGB(const Storage3D<T, 0>& Field, const Settings& locSettings, const std::string& title, const int& tStep)
    {
        std::vector<VTK::Field_t> ListOfFields;
        ListOfFields.push_back((VTK::Field_t){title, [&Field](int i,int j,int k){return T(Field(i,j,k));}});
        std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, title, tStep, ".vts");

        VTK::Write(Filename, locSettings, ListOfFields);
    };

    static void WriteSymmetryVariantsStatistics(size_t pIndex, const PhaseField& Phase, const SymmetryVariants& SV, const int tStep);
    static double GrainBoundaryStatistics(PhaseField& Phase, std::vector<dVector3> Facets, double DegreeTolerance);
    static void WriteGrainsStatistics(const PhaseField& Phase, const int tStep);

    static void GrainSizeDistribution();
    static void GrainTopologyStatistics();
    static void GrainsSurfaceArea(const PhaseField& Phase);                     ///< Calculates approximate grains surface area.

    static dVector3 FindValuePosition(PhaseField& Phi, size_t index, double value, dVector3 start_position, dVector3 direction, double tolerance = 1.0e-4);
};

}
#endif //MICROSTRUCTUREANALYSIS_H
