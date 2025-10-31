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
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef MECHANICALLOADS_H
#define MECHANICALLOADS_H

#include "Includes.h"

namespace openphase
{
class Settings;
class Temperature;
class PhaseField;
class RunTimeControl;
class ElasticProperties;
class ElasticitySolverSpectral;
class BoundaryConditions;

struct MechanicalLoadStructure
{
    EventTriggers TriggerON;                                                    ///< Load ON trigger condition selector
    EventTriggers TriggerOFF;                                                   ///< Load OFF trigger condition selector

    bool   Active;                                                              ///< Indicate if the current load is active (true) or inactive (false)
    bool   Set;                                                                 ///< Indicates if the load values are set in ElasticProperties (to eliminate repeated load updates)
    int    Repeat;                                                              ///< How many times to repeat the ON/OFF sequence (-1 for infinite)

    int    ONindex;                                                             ///< Stress or strain tensor index for the strain or stress sensitive trigger conditions to turn the load ON
    int    OFFindex;                                                            ///< Stress or strain tensor index for the strain or stress sensitive trigger conditions to turn the load OFF

    double ONtriggerValue;                                                      ///< Trigger parameter value to turn the load ON
    double OFFtriggerValue;                                                     ///< Trigger parameter value to turn the load OFF

    dVector6 AppliedStrainMask;                                                 ///< Applied strain markers for individual strain components
    dVector6 AppliedStressMask;                                                 ///< Applied stress markers for individual stress components
    dVector6 AppliedStrainRateMask;                                             ///< Applied strain rate markers for individual strain components

    vStress AppliedStress;                                                      ///< Applied stress tensor
    vStrain AppliedStrain;                                                      ///< Applied strain tensor
    vStrain AppliedStrainRate;                                                  ///< Applied strain rate tensor

    void Read(std::ifstream& inp)                                               ///< Reads load's data from a given file stream
    {
        inp.read(reinterpret_cast<char*>(&TriggerON      ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&TriggerOFF     ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&Active         ), sizeof(bool));
        inp.read(reinterpret_cast<char*>(&Set            ), sizeof(bool));
        inp.read(reinterpret_cast<char*>(&Repeat         ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&ONindex        ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&OFFindex       ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&ONtriggerValue ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&OFFtriggerValue), sizeof(double));

        inp.read(reinterpret_cast<char*>(AppliedStressMask.data()), 6*sizeof(double));
        inp.read(reinterpret_cast<char*>(AppliedStrainMask.data()), 6*sizeof(double));
        inp.read(reinterpret_cast<char*>(AppliedStrainRateMask.data()), 6*sizeof(double));
        inp.read(reinterpret_cast<char*>(AppliedStress.data()), 6*sizeof(double));
        inp.read(reinterpret_cast<char*>(AppliedStrain.data()), 6*sizeof(double));
        inp.read(reinterpret_cast<char*>(AppliedStrainRate.data()), 6*sizeof(double));
    }
    void Write(std::ofstream& outp) const                                       ///< Writes load's data to a given file stream
    {
        outp.write(reinterpret_cast<const char*>(&TriggerON      ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&TriggerOFF     ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&Active         ), sizeof(bool));
        outp.write(reinterpret_cast<const char*>(&Set            ), sizeof(bool));
        outp.write(reinterpret_cast<const char*>(&Repeat         ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&ONindex        ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&OFFindex       ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&ONtriggerValue ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&OFFtriggerValue), sizeof(double));

        outp.write(reinterpret_cast<const char*>(AppliedStressMask.data()), 6*sizeof(double));
        outp.write(reinterpret_cast<const char*>(AppliedStrainMask.data()), 6*sizeof(double));
        outp.write(reinterpret_cast<const char*>(AppliedStrainRateMask.data()), 6*sizeof(double));
        outp.write(reinterpret_cast<const char*>(AppliedStress.data()), 6*sizeof(double));
        outp.write(reinterpret_cast<const char*>(AppliedStrain.data()), 6*sizeof(double));
        outp.write(reinterpret_cast<const char*>(AppliedStrainRate.data()), 6*sizeof(double));
    }
};

class OP_EXPORTS MechanicalLoads : public OPObject                              ///< Mechanical loads class
{
 public:
    MechanicalLoads(){};                                                        ///< Constructor
    MechanicalLoads(Settings& locSettings,
                std::string InputFileName = DefaultInputFileName)               ///< Initializes storages, sets internal variables
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Allocates memory, initializes the settings
    void ReadInput(const std::string FileName) override;                        ///< Reads input parameters
    void ReadInput(std::stringstream& FileName) override;                       ///< Reads input parameters

    bool Read(const Settings& locSettings, const int tStep);                    ///< Reads raw data from a file
    bool Write(const Settings& locSettings, const int tStep) const override;    ///< Writes raw data to a file
    bool Read(const std::string FileName);                                      ///< Reads raw data from a file
    bool Write(const std::string FileName) const;                               ///< Writes raw data to a file

    void Apply(PhaseField& Phase, ElasticProperties& EP, RunTimeControl& RTC);  ///< Applies active mechanical loads
    void Activate(PhaseField& Phase, ElasticProperties& EP, RunTimeControl& RTC);///< Activates mechanical loads based on selected trigger conditions

    size_t Add(MechanicalLoadStructure& myLoad)                                 ///< Adds custom mechanical load to the list of loads
    {
        Loads.push_back(myLoad);
        return Loads.size() - 1;
    }
    std::vector<MechanicalLoadStructure> Loads;                                 ///< List of active mechanical loads
    bool ConsiderLoads;                                                         ///< Indicates if there are active mechanical loads
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
