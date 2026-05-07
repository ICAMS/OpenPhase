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
 *  File created :   2024
 *  Main contributors :  Oleg Shchyglo;
 *
 */

#ifndef THERMOKINETICINTERFACE_H
#define THERMOKINETICINTERFACE_H


#include "Includes.h"

namespace openphase
{

class Settings;
class Composition;
class Temperature;
class DiffusionProperties;


class OP_EXPORTS ThermokineticInterface : public OPObject
{
 public:

    ThermokineticInterface(){};                                                 ///< Empty constructor
    ThermokineticInterface(Settings& locSettings,
                const std::string InputFileName = DefaultInputFileName)         ///< Constructor. Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }

    ~ThermokineticInterface(void)
    {
        ConsoleOutput::WriteStandard(thisclassname, "Exited normally");
    };                                                                          ///< Destructor

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from the input stream

    void ReadTDB();                                ///< Read database file
    void ReadTDB(std::stringstream& inp);                                       ///< Read database data from the file stream

    void SetDiffusionData(PhaseField& Phase, Composition& Cx, Temperature& Tx,
                          DiffusionProperties& TP);                             ///< Set equilibrium data for equilibrium partitioning

 private:

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical components

    std::vector<std::string> PhaseNames;                                        ///< Names of thermodynamic phases
    std::vector<std::string> ElementNames;                                      ///< Names of chemical elements

    Tensor<double, 3> DiffusionCoefficients;                                    ///< Diffusion coefficients
    Tensor<double, 3> ActivationEnergies;                                       ///< Activation energies for Arrhenius equation

    std::string DBfile;                                                         ///< Database file name
    std::filesystem::path BaseDir;
};
}
#endif//THERMODYNAMICINTERFACELPD_EXP_H
