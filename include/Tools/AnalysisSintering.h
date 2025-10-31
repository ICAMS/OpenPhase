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

 *   File created :   2023
 *   Main contributors : Raphael Schiedung
 *
 */

#ifndef ANALYSISSINTERING_H
#define ANALYSISSINTERING_H

#include<string>
#include<sstream>

#include "OPObject.h"
#include "PhaseField.h"

namespace openphase
{

class DoubleObstacle;
class GrandPotentialDensity;
class GrandPotentialSolver;
class InterfaceProperties;
class Settings;

struct OP_EXPORTS AnalysisSintering: public OPObject
{
    AnalysisSintering(){};
    AnalysisSintering(Settings& locSettings, std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;                            ///< Initializes internal variables and storages
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input data from the user specified input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input data from the user specified string stream

    using Field3D = std::function<double(long,long,long)>;
    double CalculateDensity_BoundingBox      (const Field3D& MassDensity) const;
    static double CalculateDensity_LineIntersection (const PhaseField &Phase, const Field3D& MassDensity);
    static double CalculateDensity_LineIntersectionX(const PhaseField &Phase, const Field3D& MassDensity);
    static double CalculateDensity_LineIntersectionY(const PhaseField &Phase, const Field3D& MassDensity);
    static double CalculateDensity_LineIntersectionZ(const PhaseField &Phase, const Field3D& MassDensity);

    void PrintSolidPhaseFractionInBoundingBox(const PhaseField& Phase, const GrandPotentialDensity& omega, const GrandPotentialSolver& GPS);
    void WriteSolidPhaseFractionInBoundingBox(const PhaseField& Phase, const GrandPotentialDensity& omega, const GrandPotentialSolver& GPS, std::string filename, double time, char separator);

    void DoDiagnostics(Settings& locSettings, PhaseField& Phase,
            DoubleObstacle& DO, InterfaceProperties& IP,
            GrandPotentialDensity& omega, GrandPotentialSolver& GPS,
            RunTimeControl& RTC);

private:
    GridParameters Grid;                                                        ///< Simulation grid parameters

    double ReferenceDensity;                                                    ///< Reference density for relative Density calculation
    bool Do_BoundingBox;                                                        ///< True if bounding box is used for density calculation
    bool Do_LineIntersection;                                                   ///< True if Line Intersection averaged over the results in x,y,z direction is used for density calculation
    bool Do_LineIntersectionX;                                                  ///< True if Line Intersection in x-direction is used for density calculation
    bool Do_LineIntersectionY;                                                  ///< True if Line Intersection in y-direction is used for density calculation
    bool Do_LineIntersectionZ;                                                  ///< True if Line Intersection in z-direction is used for density calculation
    long BoundingBoxSizeX;                                                      ///< Size of bounding box for density calculation in x-direction
    long BoundingBoxSizeY;                                                      ///< Size of bounding box for density calculation in x-direction
    long BoundingBoxSizeZ;                                                      ///< Size of bounding box for density calculation in x-direction
    long locBoundingBox_min_i;                                                  ///< First local x-coordinate for integration
    long locBoundingBox_max_i;                                                  ///< Last local x-coordinate for integration
    long locBoundingBox_min_j;                                                  ///< First local y-coordinate for integration
    long locBoundingBox_max_j;                                                  ///< Last local y-coordinate for integration
    long locBoundingBox_min_k;                                                  ///< First local z-coordinate for integration
    long locBoundingBox_max_k;                                                  ///< Last local z-coordinate for integration
    double BoundingBoxVolume;
    void SetBoundingBoxCoordinates();

    static long First_SolidXCoordinate(const PhaseField &Phase, size_t GasPhaseIdx, long j, long k); ///< Return the first x coordinate which is no gas
    static long Last_SolidXCoordinate (const PhaseField &Phase, size_t GasPhaseIdx, long j, long k); ///< Return the first x coordinate which is no gas
};
}
#endif

