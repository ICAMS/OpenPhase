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

 *   File created :   2015
 *   Main contributors :   Philipp Engels; Hesham Salama
 *
 */

#ifndef CRYSTALLOGRAPHY_H
#define CRYSTALLOGRAPHY_H

#include "Includes.h"

namespace openphase
{

class Settings;

class OP_EXPORTS Crystallography : public OPObject
{
 public:

    LatticeSystems CrystalSym;

    Crystallography(){};
    Crystallography(Settings& locSettings, std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes the module, allocate the storage, assign internal variables
    void ReadInput(std::string InputFileName) override;                         ///< Reads input parameters

    inline static const std::vector<double> SlipSystem_ISO =
    { //  {  n1,  n2,  n3,  d1,  d2,  d3, }
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    };
    
    inline static const std::vector<double> SlipSystem_Cubic =
    { //  {  n1,  n2,  n3,  d1,  d2,  d3, }
        0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
        0.0, 1.0,-1.0, 1.0, 0.0, 0.0,
    
        1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
        -1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    
        1.0, 1.0, 0.0, 0.0, 0.0, 1.0,
        -1.0, 1.0, 0.0, 0.0, 0.0, 1.0
    };
    inline static const std::vector<double> SlipSystem_BCC =
    { //  {  n1,  n2,  n3,  d1,  d2,  d3 }
    
        // {110}<111> slip systems
        0.0, 1.0,-1.0, 1.0, 1.0, 1.0,
        1.0, 0.0,-1.0, 1.0, 1.0, 1.0,
        1.0,-1.0, 0.0, 1.0, 1.0, 1.0,
    
        0.0, 1.0,-1.0,-1.0, 1.0, 1.0,
        1.0, 0.0, 1.0,-1.0, 1.0, 1.0,
        1.0, 1.0, 0.0,-1.0, 1.0, 1.0,
    
        0.0, 1.0, 1.0, 1.0,-1.0, 1.0,
        1.0, 0.0,-1.0, 1.0,-1.0, 1.0,
        1.0, 1.0, 0.0, 1.0,-1.0, 1.0,
    
        0.0, 1.0, 1.0, 1.0, 1.0,-1.0,
        1.0, 0.0, 1.0, 1.0, 1.0,-1.0,
        1.0,-1.0, 0.0, 1.0, 1.0,-1.0,
    
        // {211}<111> slip systems
        -2.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0,-2.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0,-2.0, 1.0, 1.0, 1.0,
    
        2.0, 1.0, 1.0,-1.0, 1.0, 1.0,
        1.0, 2.0,-1.0,-1.0, 1.0, 1.0,
        1.0,-1.0, 2.0,-1.0, 1.0, 1.0,
    
        2.0, 1.0,-1.0, 1.0,-1.0, 1.0,
        1.0, 2.0, 1.0, 1.0,-1.0, 1.0,
        -1.0, 1.0, 2.0, 1.0,-1.0, 1.0,
    
        2.0,-1.0, 1.0, 1.0, 1.0,-1.0,
        -1.0, 2.0, 1.0, 1.0, 1.0,-1.0,
        1.0, 1.0, 2.0, 1.0, 1.0,-1.0,
    
    };
    inline static const std::vector<double> SlipSystem_FCC_XY =
    {
        1.0, 1.0, 0.0, 1.0,-1.0, 0.0,
        1.0, 1.0, 0.0,-1.0, 1.0, 0.0,
        -1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
        -1.0, 1.0, 0.0,-1.0,-1.0, 0.0,
    };
    inline static const std::vector<double> SlipSystem_FCC_XZ =
    {
        1.0, 0.0, 1.0, 1.0, 0.0,-1.0,
        1.0, 0.0, 1.0,-1.0, 0.0, 1.0,
        -1.0, 0.0, 1.0, 1.0, 0.0, 1.0,
        -1.0, 0.0, 1.0,-1.0, 0.0,-1.0,
    };
    inline static const std::vector<double> SlipSystem_FCC_YZ =
    {
        0.0, 1.0, 1.0, 0.0, 1.0,-1.0,
        0.0, 1.0, 1.0, 0.0,-1.0, 1.0,
        0.0,-1.0, 1.0, 0.0, 1.0, 1.0,
        0.0,-1.0, 1.0, 0.0,-1.0,-1.0
    };
    inline static const std::vector<double> SlipSystem_FCC =
    {
        1.0, 1.0, 1.0, 1.0,-1.0, 0.0,
        -1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
        1.0, 1.0,-1.0, 1.0,-1.0, 0.0,
        1.0,-1.0, 1.0, 1.0, 1.0, 0.0,
    
        1.0,-1.0, 1.0, 1.0, 0.0,-1.0,
        -1.0, 1.0, 1.0, 1.0, 0.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 0.0,-1.0,
        1.0, 1.0,-1.0, 1.0, 0.0, 1.0,
    
        -1.0, 1.0, 1.0, 0.0, 1.0,-1.0,
        1.0, 1.0,-1.0, 0.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 0.0, 1.0,-1.0,
        1.0,-1.0, 1.0, 0.0, 1.0, 1.0
    };
    inline static const std::vector<double> SlipSystem_PartialFCC =
    { //  {  n1,  n2,  n3,  d1,  d2,  d3 }
        -2.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0,-2.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0,-2.0, 1.0, 1.0, 1.0,
    
        2.0, 1.0, 1.0,-1.0, 1.0, 1.0,
        -1.0,-2.0, 1.0,-1.0, 1.0, 1.0,
        -1.0, 1.0,-2.0,-1.0, 1.0, 1.0,
    
        -2.0,-1.0, 1.0, 1.0,-1.0, 1.0,
        1.0, 2.0, 1.0, 1.0,-1.0, 1.0,
        1.0,-1.0,-2.0, 1.0,-1.0, 1.0,
    
        -2.0, 1.0,-1.0, 1.0, 1.0,-1.0,
        1.0,-2.0,-1.0, 1.0, 1.0,-1.0,
        1.0, 1.0, 2.0, 1.0, 1.0,-1.0
    };
    inline static const std::vector<double> SlipSystem_HCPBasal =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Basal systems {0001}<11-20>
        0.0, 0.0, 0.0, 1.0, 2.0,-1.0,-1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,-1.0, 2.0,-1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 2.0, 0.0
    };
    inline static const std::vector<double> SlipSystem_HCPPrismatic =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Prismatic systems {10-10}<11-20>
        0.0, 1.0,-1.0, 0.0, 2.0,-1.0,-1.0, 0.0,
        -1.0, 0.0, 1.0, 0.0,-1.0, 2.0,-1.0, 0.0,
        1.0,-1.0, 0.0, 0.0,-1.0,-1.0, 2.0, 0.0
    };
    inline static const std::vector<double> SlipSystem_HCPPyramidal =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Pyramidal systems <a> {10-11}<11-20>
        0.0, 1.0,-1.0, 1.0, 2.0,-1.0,-1.0, 0.0,
        -1.0, 0.0, 1.0, 1.0,-1.0, 2.0,-1.0, 0.0,
        1.0,-1.0, 0.0, 1.0,-1.0,-1.0, 2.0, 0.0,
        -1.0, 1.0, 0.0, 1.0, 1.0, 1.0,-2.0, 0.0,
        0.0,-1.0, 1.0, 1.0,-2.0, 1.0, 1.0, 0.0,
        1.0, 0.0,-1.0, 1.0, 1.0,-2.0, 1.0, 0.0
    };
    inline static const std::vector<double> SlipSystem_HCPPyramidal_1st =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Pyramidal systems (1st) <c+a> {10-11}<11-23>
        -1.0, 1.0, 0.0, 1.0, 2.0, -1.0,-1.0, 3.0,
        -1.0, 1.0, 0.0, 1.0, 1.0, -2.0, 1.0, 3.0,
        1.0, 0.0,-1.0, 1.0,-1.0, -1.0, 2.0, 3.0,
    
        1.0, 0.0,-1.0, 1.0,-2.0, 1.0, 1.0, 3.0,
        0.0,-1.0, 1.0, 1.0,-1.0, 2.0,-1.0, 3.0,
        0.0,-1.0, 1.0, 1.0, 1.0, 1.0,-2.0, 3.0,
    
        1.0,-1.0, 0.0, 1.0,-2.0, 1.0, 1.0, 3.0,
        1.0,-1.0,-0.0, 1.0,-1.0, 2.0,-1.0, 3.0,
        -1.0, 0.0, 1.0, 1.0, 1.0, 1.0,-2.0, 3.0,
    
        -1.0, 0.0, 1.0, 1.0, 2.0,-1.0,-1.0, 3.0,
        0.0, 1.0,-1.0, 1.0, 1.0,-2.0, 1.0, 3.0,
        0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 2.0, 3.0
    };
    inline static const std::vector<double> SlipSystem_HCPPyramidal_2nd =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Pyramidal systems (2nd) <c+a> {11-22}<11-23>
        -2.0, 1.0, 1.0, 2.0, 2.0,-1.0,-1.0, 3.0,
        1.0,-2.0, 1.0, 2.0,-1.0, 2.0,-1.0, 3.0,
        1.0, 1.0,-2.0, 2.0,-1.0,-1.0, 2.0, 3.0,
        2.0,-1.0,-1.0, 2.0,-2.0, 1.0, 1.0, 3.0,
        -1.0, 2.0,-1.0, 2.0, 1.0,-2.0, 1.0, 3.0,
        -1.0,-1.0, 2.0, 2.0, 1.0, 1.0,-2.0, 3.0
    };
    inline static const std::vector<double> TwinSystem_FCC =
    {
        -2.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0,-2.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0,-2.0, 1.0, 1.0, 1.0,
    
        2.0, 1.0, 1.0,-1.0, 1.0, 1.0,
        -1.0,-2.0, 1.0,-1.0, 1.0, 1.0,
        -1.0, 1.0,-2.0,-1.0, 1.0, 1.0,
    
        -2.0,-1.0, 1.0, 1.0,-1.0, 1.0,
        1.0, 2.0, 1.0, 1.0,-1.0, 1.0,
        1.0,-1.0,-2.0, 1.0,-1.0, 1.0,
    
        -2.0, 1.0,-1.0, 1.0, 1.0,-1.0,
        1.0,-2.0,-1.0, 1.0, 1.0,-1.0,
        1.0, 1.0, 2.0, 1.0, 1.0,-1.0
    };
    inline static const std::vector<double> TwinSystem_HCPTensileTwin1 =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Tensile twinning in Co, Mg, Zr, Ti, and Be; compressive twinning in Cd and Zn (T1)
        -1.0, 1.0, 0.0, 2.0, 1.0,-1.0, 0.0, 1.0,
        1.0, 0.0,-1.0, 2.0,-1.0, 0.0, 1.0, 1.0,
        0.0,-1.0, 1.0, 2.0, 0.0, 1.0,-1.0, 1.0,
        1.0,-1.0, 0.0, 2.0,-1.0, 1.0, 0.0, 1.0,
        -1.0, 0.0, 1.0, 2.0, 1.0, 0.0,-1.0, 1.0,
        0.0, 1.0,-1.0, 2.0, 0.0,-1.0, 1.0, 1.0
    };
    inline static const std::vector<double> TwinSystem_HCPTensileTwin2 =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Tensile twinning in Co, Re, and Zr (T2)
        -2.0, 1.0, 1.0, 1.0, 2.0,-1.0,-1.0, 6.0,
        1.0,-2.0, 1.0, 1.0,-1.0, 2.0,-1.0, 6.0,
        1.0, 1.0,-2.0, 1.0,-1.0,-1.0, 2.0, 6.0,
        2.0,-1.0,-1.0, 1.0,-2.0, 1.0, 1.0, 6.0,
        -1.0, 2.0,-1.0, 1.0, 1.0,-2.0, 1.0, 6.0,
        -1.0,-1.0, 2.0, 1.0, 1.0, 1.0,-2.0, 6.0
    };
    inline static const std::vector<double> TwinSystem_HCPCompressionTwin1 =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Compressive twinning in Ti and Zr (C1)
        2.0,-1.0,-1.0, 2.0, 2.0,-1.0,-1.0,-3.0,
        -1.0, 2.0,-1.0, 2.0,-1.0, 2.0,-1.0,-3.0,
        -1.0,-1.0, 2.0, 2.0,-1.0,-1.0, 2.0,-3.0,
        -2.0, 1.0, 1.0, 2.0,-2.0, 1.0, 1.0,-3.0,
        1.0,-2.0, 1.0, 2.0, 1.0,-2.0, 1.0,-3.0,
        1.0, 1.0,-2.0, 2.0, 1.0, 1.0,-2.0,-3.0
    };
    inline static const std::vector<double> TwinSystem_HCPCompressionTwin2 =
    { //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
        //Compressive twinning in Mg (C1)
        -1.0, 1.0, 0.0,-2.0,-1.0, 1.0, 0.0, 1.0,
        1.0, 0.0,-1.0,-2.0, 1.0, 0.0,-1.0, 1.0,
        0.0,-1.0, 1.0,-2.0, 0.0,-1.0, 1.0, 1.0,
        1.0,-1.0, 0.0,-2.0, 1.0,-1.0, 0.0, 1.0,
        -1.0, 0.0, 1.0,-2.0,-1.0, 0.0, 1.0, 1.0,
        0.0, 1.0,-1.0,-2.0, 0.0, 1.0,-1.0, 1.0
    };

    Storage<dMatrix3x3> SymmetriesCubic;                                        ///< Stores cubic symmetry operations
    Storage<dMatrix3x3> CrystalSymmetries;

    size_t nsym = 1;

    double Phi_min = 0.0;
    double Phi_max = 0.0;
    double Theta_min = 0.0;
    double Theta_max = 0.0;

    const double a = std::sqrt(3.0) / 2.0;

 protected:
 private:
};

} // namespace openphase

#endif
