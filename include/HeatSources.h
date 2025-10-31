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
 *   Main contributors :   Oleg Shchyglo; Marvin Tegeler; Matthias Stratmann
 *
 */

#ifndef HEATSOURCES_H
#define HEATSOURCES_H

#include "Includes.h"

namespace openphase
{
class Settings;
class HeatDiffusion;
class Temperature;
class PhaseField;
class RunTimeControl;
class BoundaryConditions;

enum class HeatSourceTypes                                                      ///< Types of available heat sources
{
    BCX0,                                                                       ///< Heat source as boundary condition for lower X boundary
    BCXN,                                                                       ///< Heat source as boundary condition for upper X boundary
    BCY0,                                                                       ///< Heat source as boundary condition for lower Y boundary
    BCYN,                                                                       ///< Heat source as boundary condition for upper Y boundary
    BCZ0,                                                                       ///< Heat source as boundary condition for lower Z boundary
    BCZN,                                                                       ///< Heat source as boundary condition for upper Z boundary
    Phase,                                                                      ///< Selected phase will act as a heat source
    Ellipsoidal                                                                 ///< Ellipsoidal shape heat source. Position and size specified by the user
};

enum class HeatSourceEnergyDistributions                                        ///< Energy distributions over the heat source volume
{
    Uniform,                                                                    ///< Uniform energy distribution
    Gaussian                                                                    ///< Gaussian energy distribution
};

enum class HeatSourceModes                                                      ///< Heat source modes
{
    ThreeDimensional,                                                           ///< Volumetric heat source
    TwoDimensional                                                              ///< Area heat source
};

struct HeatSourceStructure
{
    HeatSourceTypes Type;                                                       ///< Heat source type
    EventTriggers TriggerON;                                                    ///< Heat source ON trigger condition selector
    EventTriggers TriggerOFF;                                                   ///< Heat source OFF trigger condition selector
    HeatSourceEnergyDistributions Distribution;                                 ///< Energy distribution over the heat source volume
    HeatSourceModes Mode;                                                       ///< Heat source mode: volumetric or area

    bool Active;                                                                ///< Indicate if the current heat source is active (true) or inactive (false)
    int Repeat;                                                                 ///< How many times to repeat the ON/OFF sequence (-1 for infinite)

    double Value;                                                               ///< Heat flux value [W/m^2] or [W/m^3] depending on the user settings
    double ONtriggerValue;                                                      ///< Trigger parameter value to turn the heat source ON
    double OFFtriggerValue;                                                     ///< Trigger parameter value to turn the heat source OFF
    size_t PhaseIndex;                                                          ///< Phase index for the phase sensitive trigger conditions
    size_t PhaseIndexON;                                                        ///< Phase index for the phase fraction evaluation for ON trigger
    size_t PhaseIndexOFF;                                                       ///< Phase index for the phase fraction evaluation for OFF trigger
    dVector3 Position;                                                          ///< Position of the center of the heat source for ellipsoidal or rectangular shapes [grid coordinates]
    dVector3 Size;                                                              ///< Heat source size for ellipsoidal (three orthogonal radii) and rectangular (three orthogonal side lengths) shape sources [grid coordinates]
    dVector3 Velocity;                                                          ///< Heat source sweep velocity (only for ellipsoidal and rectangular types) [true units]
    void Read(std::ifstream& inp)                                               ///< Reads heat sources data from a given file stream
    {
        inp.read(reinterpret_cast<char*>(&Type           ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&TriggerON      ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&TriggerOFF     ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&Distribution   ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&Mode           ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&Active         ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Repeat         ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&Value          ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&ONtriggerValue ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&OFFtriggerValue), sizeof(double));

        inp.read(reinterpret_cast<char*>(&PhaseIndex     ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&PhaseIndexON   ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&PhaseIndexOFF  ), sizeof(size_t));

        inp.read(reinterpret_cast<char*>(&Position[0]), 3*sizeof(double));
        inp.read(reinterpret_cast<char*>(&Size    [0]), 3*sizeof(double));
        inp.read(reinterpret_cast<char*>(&Velocity[0]), 3*sizeof(double));
    }
    void Write(std::ofstream& outp) const                                       ///< Writes heat sources data to a given file stream
    {
        outp.write(reinterpret_cast<const char*>(&Type           ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&TriggerON      ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&TriggerOFF     ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&Distribution   ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&Mode           ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&Active         ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Repeat         ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&Value          ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&ONtriggerValue ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&OFFtriggerValue), sizeof(double));

        outp.write(reinterpret_cast<const char*>(&PhaseIndex     ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&PhaseIndexON   ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&PhaseIndexOFF  ), sizeof(size_t));

        outp.write(reinterpret_cast<const char*>(&Position[0]), 3*sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Size    [0]), 3*sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Velocity[0]), 3*sizeof(double));
    }
};

class OP_EXPORTS HeatSources : public OPObject                                  ///< Heat sources class
{
 public:
    HeatSources(){};                                                            ///< Constructor
    HeatSources(Settings& locSettings,
                std::string InputFileName = DefaultInputFileName)               ///< Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Allocates memory, initializes the settings);
    void ReadInput(const std::string FileName) override;                        ///< Reads input parameters
    void ReadInput(std::stringstream& FileName) override;                       ///< Reads input parameters

    bool Read(const Settings& locSettings, const int tStep);                    ///< Reads raw data from a file
    bool Write(const Settings& locSettings, const int tStep) const override;    ///< Writes raw data to a file
    bool Read(const std::string FileName);                                      ///< Reads raw data from a file
    bool Write(const std::string FileName) const;                               ///< Writes raw data to a file

    void Apply(PhaseField& Phase, Temperature& Tx, HeatDiffusion& HD);          ///< Applies active heat sources
    void Activate(PhaseField& Phase, Temperature& Tx, RunTimeControl& RTC);     ///< Activates heat sources based on selected trigger conditions

    size_t Add(HeatSourceStructure& myHeatSource)                               ///< Adds custom heat source to the list of sources
    {
        Sources.push_back(myHeatSource);
        return Sources.size() - 1;
    }
    std::vector<HeatSourceStructure> Sources;                                   ///< List of active heat sources
    bool ConsiderHeatSources;                                                   ///< Indicated if there are active heat sources
    size_t Nphases;                                                             ///< Number of phases

 protected:
 private:
    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC,
              const int tStep) override                                         ///< Reads raw data from a file
    {
        // Does nothing and should not be used
        return false;
    }
};

} // namespace openphase
#endif
