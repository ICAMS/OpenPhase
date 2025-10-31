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
 *   Main contributors :   Matthias Stratmann; Johannes Goerler
 *
 */

#ifndef TEXTOUTPUT_H
#define TEXTOUTPUT_H

#include "Includes.h"

namespace openphase
{
class BoundaryConditions;
class ElasticProperties;
class GrainsProperties;
class InterfaceField;
class PhaseField;
class Settings;
class Temperature;
class Composition;
class Nucleation;
class DrivingForce;
class Orientations;

class TextOutput //: public OPObject                                                               ///< Output of tabulated simulation data
{
 public:
    static void WriteValue(double Value,std::string filename, double time);                        ///< Write double value in column
    static void WriteMultipleValues(std::vector<std::string> Names, std::vector<double> value, 
            std::string filename, double time);                                                    ///< Write vector values and names of value in columns
    static void WritedVectorNValues(std::vector<std::string> Names, dVectorN value, 
            std::string filename, double time);                                                    ///< Write dVectorN values and names of values in columns
    static void LineConcentration(Composition& Cx, PhaseField& Phi,
                                  std::string filename, double timestep,
                                  std::string type, std::string axis, int x, int y,int z);         ///< Write total composition over a straight line in separate files
    static void AveragePlasticStrain(ElasticProperties &EP, std::string filename, 
                                     double timeOrStrain, size_t precision = 5);                   ///< Write average Plastic strain over time
    static void maxElasticRotation(ElasticProperties& EP,
           PhaseField& Phase,Orientations& OR, int tStep, std::string Filename);                   ///< Write max elastic rotations over time
    static NodeAB<double,double> IntegrateDG(PhaseField& Phase, DrivingForce& dG, std::string filename, double RealTime);
    static NodeAB<double,double> IntegrateDG(PhaseField& Phase, DrivingForce& dG, std::string filename, NodeAB<double,double> dgold, double RealTime);
 private:
    static constexpr auto separator = ' ';
};

} // namespace openphase
#endif
