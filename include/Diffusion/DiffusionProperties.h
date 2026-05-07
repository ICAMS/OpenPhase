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
 */

#ifndef DIFFUSIONPROPERTIES_H
#define DIFFUSIONPROPERTIES_H

#include "Includes.h"

namespace openphase
{
class Settings;
class PhaseField;
class Composition;
class Temperature;
class BoundaryConditions;
class ChemicalDiffusion;
class ThermokineticInterface;
class ThermodynamicInterfaceLPD;
class ThermodynamicInterfaceOPSTC;
class ThermodynamicPropertiesEQP;

class DiffusionData                                                            ///< Structure for storing pair equilibrium data
{
 public:
    double eq_temperature;                                                      ///< Equilibrium temperature at which the data is taken                                                  ///< Molar heat capacity of the phase
    Tensor<double,1> nom_composition;                                           ///< Nominal composition for the local equilibrium                                            ///< Slope of the solvus line of the phase
    bool valid;                                                                 ///< True if the dataset has valid data
    bool out_of_T_range;                                                        ///< True if dataset is out of temperature range
    bool out_of_C_range;                                                        ///< True if dataset is out of composition range

    ~DiffusionData()
    {

    }
    DiffusionData():
        eq_temperature(0.0),
        nom_composition(),
        valid(false),
        out_of_T_range(false),
        out_of_C_range(false)
    {

    }
    DiffusionData(size_t Ncomp):
        eq_temperature(0.0),
        nom_composition({Ncomp}),
        valid(false),
        out_of_T_range(false),
        out_of_C_range(false)
    {

    }
    DiffusionData(const DiffusionData& rhs):
        eq_temperature(rhs.eq_temperature),
        nom_composition(rhs.nom_composition),
        valid(rhs.valid),
        out_of_T_range(rhs.out_of_T_range),
        out_of_C_range(rhs.out_of_C_range)
    {

    }
    DiffusionData(DiffusionData&& rhs):
        eq_temperature(rhs.eq_temperature),
        nom_composition(rhs.nom_composition),
        valid(rhs.valid),
        out_of_T_range(rhs.out_of_T_range),
        out_of_C_range(rhs.out_of_C_range)
    {

    }
    DiffusionData& operator=(const DiffusionData& rhs)
    {
        eq_temperature  = rhs.eq_temperature;
        nom_composition = rhs.nom_composition;
        valid           = rhs.valid;
        out_of_T_range  = rhs.out_of_T_range;
        out_of_C_range  = rhs.out_of_C_range;

        return *this;
    }
    DiffusionData& operator=(DiffusionData&& rhs)
    {
        eq_temperature  = rhs.eq_temperature;
        nom_composition = rhs.nom_composition;
        valid           = rhs.valid;
        out_of_T_range  = rhs.out_of_T_range;
        out_of_C_range  = rhs.out_of_C_range;
        return *this;
    }
    void Allocate(size_t n_comp)
    {
        nom_composition.Allocate({n_comp});
    }
};

class OP_EXPORTS DiffusionProperties : public OPObject                          ///< Stores and processes diffusion parameters
{
 friend ChemicalDiffusion;
 friend ThermokineticInterface;
 friend ThermodynamicInterfaceLPD;
 friend ThermodynamicInterfaceOPSTC;
 friend ThermodynamicPropertiesEQP;

 public:
    DiffusionProperties(){};
    DiffusionProperties(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Allocates memory and initializes settings
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input from the file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input from the file
    void SetBoundaryConditions(const BoundaryConditions& BC) override;
    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;                         ///< Remeshes the diffusion coefficients storage while keeping the data

    void WriteVTK(const Settings& locSettings, const Composition& Cx, const int tStep);
    void CheckRange(PhaseField& Phi, Composition& Cx, Temperature& Tx);

    void PrintPointStatistics(int i, int j, int k);

    double MaxDiffusionCoefficient;                                             ///< Maximum diffusion coefficient value in the simulation domain

    Storage3D < double, 3 > DiffusionCoefficients;                              ///< Diffusion coefficients storage

    Tensor<double,2> GrainBoundaryFactor;                                       ///< Grain boundary diffusion enhancement factor
    Tensor<DiffusionData,1> GlobalExtrapolationData;                            ///< Global extrapolation data
    Storage3D<DiffusionData,1> LocalExtrapolationData;                          ///< Local extrapolation data

 protected:
    ExtrapolationModes ExtrapolationMode;
    double MaxMoleFractionDeviation;                                            ///< Maximum mole fraction deviation for diffusion data validity
    double MaxTemperatureDeviation;                                             ///< Maximum temperature deviation for diffusion data validity
    Tensor<double,3> ScalingFactor;                                             ///< Diffusion coefficients scaling factor

 private:
    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;
    size_t Ncomp;
    std::vector<std::string> PhaseNames;
    std::vector<std::string> ElementNames;
};

}
#endif//DIFFUSIONPROPERTIES_H
