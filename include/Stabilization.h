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

 *   File created :   
 *   Main contributors :
 *
 */

#ifndef STABILIZATION_H
#define STABILIZATION_H

#include "Includes.h"
namespace openphase
{

class BoundaryConditions;
class InterfaceProperties;
class PhaseField;
class Settings;
class DrivingForce;
class StabilizationImpl;
/***************************************************************/
class OP_EXPORTS Stabilization
{
 public:
    Stabilization();
    Stabilization(Settings& locSettings);
    ~Stabilization();
    void Initialize(Settings& locSettings);
    void Initialize(Settings& locSettings, LaplacianStencil _DStencil);
    void ReadInput(const std::string FileName);
    void ReadInput(std::stringstream& inp);                    
    void CalculateCurvature(PhaseField& Phase);
    void CalculateCurvature_SPF(PhaseField& Phase);
    void CalculateAverageCurvature(PhaseField& Phase);
    void CalculateAverageCurvatureRegularization(PhaseField& Phase);
    void CalculateAverageCurvatureNormal(PhaseField& Phase);
    void CalculateAverageCurvatureNormalSnap(PhaseField& Phase);
    void CalculateAverageCurvatureSimple(PhaseField& Phase);
    void ComputeVelocity(PhaseField& Phase,
                                InterfaceProperties& IP, DrivingForce& dG);
    void ComputeUpdate(PhaseField& Phase,
                                DrivingForce& dG, double dt);
    double GetMaxDT();
    void Smooth(PhaseField& Phase);
    void WriteVTK(const Settings& locSettings, const int tStep,
                            const size_t indexA, const size_t indexB,
                            const int precision = 16) const;
 protected:
 private:
    std::unique_ptr<StabilizationImpl> impl_;
};
}
#endif
