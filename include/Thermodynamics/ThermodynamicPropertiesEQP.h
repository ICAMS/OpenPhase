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
 *  File created :   2011
 *  Main contributors :   Oleg Shchyglo; Marvin Tegeler; Matthias Stratmann;
 *                         Efim Borukhovich;
 *
 */

#ifndef THERMODYNAMICPROPERTIESEQP_H
#define THERMODYNAMICPROPERTIESEQP_H

#include "Includes.h"
#include "Settings.h"

namespace openphase
{
class BoundaryConditions;
class Composition;
class DiffusionProperties;
class ThermodynamicInterfaceLPD;
class ThermodynamicInterfaceOPSTC;
class ThermodynamicInterfaceOPSCA;
class Nucleation;
class PairEquilibriumData;
class InterfaceProperties;
class DrivingForce;
class Temperature;

enum class DrivingForceModels
{
    Standard,
    LowerSlope,
    Weighted,
    User
};

class OP_EXPORTS ThermodynamicPropertiesEQP : public OPObject               ///< Stores and manipulates thermodynamic properties for equilibirum partitioning model
{
 friend ThermodynamicInterfaceLPD;
 friend ThermodynamicInterfaceOPSTC; 
 friend ThermodynamicInterfaceOPSCA; 
	
 public:

    ThermodynamicPropertiesEQP()
    {

    };
    ThermodynamicPropertiesEQP(Settings& locSettings,
                               const std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from the input stream

    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;                         ///< Remeshes the storages
    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets boundary conditions

    void PrintStatistics(Composition& Cx,
                            bool matrix, bool warnings, double tStep);

    void CheckRange(PhaseField& Phi, Composition& Cx, Temperature& Tx);

    void CalculateDrivingForce(PhaseField& Phi, InterfaceProperties& IP,
                               Composition& Cx, Temperature& Tx,
                               DrivingForce& dGab);                             ///< Calculates thermodynamic driving force

    void PrintPointStatistics(int i, int j, int k);                             ///< Prints thermodynamic properties in a given point (x, y, z)
    void WriteVTK(const Settings& locSettings,
                  const int tStep,
                  const int precision=16) const;                                ///< Writes equilibrium partitioning data in VTK format (.vts file)

    const Tensor<EquilibriumData,2>& ExtrapolationData(int i, int j, int k) const
    {
        if(ExtrapolationMode == ExtrapolationModes::Global)
        {
            return GlobalExtrapolationData;
        }
        return LocalExtrapolationData(i,j,k);
    }
    Tensor<EquilibriumData,2> GlobalExtrapolationData;                          ///< Global extrapolation data
    Storage3D<EquilibriumData,2> LocalExtrapolationData;                        ///< Local extrapolation data
    Storage3D< double, 2> dMu;                                                  ///< External chemical potential contribution

    void ClearChemicalPotential(void);                                          ///< Clears chemical potential storage

 protected:
    ExtrapolationModes ExtrapolationMode;

 private:
    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical components

    std::vector<std::string> PhaseNames;                                        ///< Names of thermodynamic phases
    std::vector<std::string> ElementNames;                                      ///< Names of chemical elements
    
    Tensor<double,2> IntersectionTemperature;
    Tensor<double,2> ReferenceTemperature;

    bool AtStart;
    //bool Projection = true;

    Table<bool> EnableDrivingForce;                                             ///< Switches phase pairs driving force contribution on/off
    DrivingForceModels DrivingForceModel;                                       ///< Driving force calculation model
    Table<bool> DrivingForceMap;                                                ///< Driving force inclusion map for user driving force model

    double MaxCompositionDeviation;                                             ///< Max mole fraction deviation for thermodynamic data validity
    double MaxTemperatureDeviation;                                             ///< Max temperature deviation for thermodynamic data validity
};

}
#endif//THERMODYNAMICPROPERTIESEQP_H
