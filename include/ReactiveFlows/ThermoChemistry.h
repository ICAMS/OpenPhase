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
 *  File created :   2025
 *  Main contributors :   Oleg Shchyglo; Reza Namdar
 *
 */
#ifndef THERMOCHEMISTRY_H_INCLUDED
#define THERMOCHEMISTRY_H_INCLUDED  

#ifdef CANTERA

#include <fstream>
#include <istream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>

#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Velocities.h"
#include "BoundaryConditions.h"
#include "VTK.h"


#include "cantera/base/Array.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/transport/TransportFactory.h"
#include <cantera/core.h>

#include "Species.h"
#include "Energy.h"
#ifdef ONNX
#include "NeuralNetwork.h"
#endif
using namespace std;


namespace openphase {class Species;

class Energy;
#ifdef ONNX
class NeuralNetwork;
#endif
}

namespace openphase {
class OP_EXPORTS ThermoChemistry
{
    public:
    ThermoChemistry(const std::string InputFileName)                                   ///< Constructor
    {
        ReadInput(InputFileName);
    }
    void ReadInput(string InputFile);
    void WriteSpeciesNamesInSettingsInputFile(string InputFile);
    
    void InitializeGasPhase(Cantera::ThermoPhase& gas, const double T0, const double P0, string SettingsInputFile);
#ifdef ONNX
    void CalculateGasKinetics(Energy &EN, Species& SP, MixtureFlow &MF,
                                       FlowSolverLBM& FL, NeuralNetwork& NN);
#endif

    void CalculateGasKinetics(Energy &EN, Species& SP, MixtureFlow &MF, FlowSolverLBM& FL,
						Cantera::ThermoPhase &gas, Cantera::Transport &transp, Cantera::Kinetics &kin, BoundaryConditions& BC);
    void CalculateGasKinetics(PhaseField& Phase, Energy &EN, Species& SP, MixtureFlow &MF, SolidBody& SB, FlowSolverLBM& FL,
						Cantera::ThermoPhase &gas, Cantera::Transport &transp, Cantera::Kinetics &kin, BoundaryConditions& BC);
    Tensor<double, 1> GettingGasFuelMixture(Cantera::ThermoPhase &gas, const std::string FuelName, const std::string Oxidizer, double EquiValRatio);
    Tensor<double, 1> GettingGasMixture(Cantera::ThermoPhase &gas, const double T0, const double P0);
    Tensor<double, 1> GettingAirMixture();
    Tensor<double, 1> GettingBackFlowMixture();
    Tensor<double, 1> GettingBurntGasMixture(Cantera::ThermoPhase& gas, const double T0, const double P0);
    Tensor<double, 1> GettingGasMolecularweight(Cantera::ThermoPhase& gas);
    Tensor<string, 1> GettingGasSpeciesNames(Cantera::ThermoPhase& gas);

    size_t GettingGasFuelIndex(Cantera::ThermoPhase& gas, string FuelName);
    void ExportingGasSpeciesData(Species &SP);
    vector<vector<double>> GettingGasSpeciesPolyConstants(Cantera::ThermoPhase &gas);

	vector<vector<double>>  Coeffs;
    size_t nCoeffs;

    double FuelConsumptionRate=0.0;

    string ReactionMechanism;
    string GasPhaseName;
    string GasTransportName;
    string InitialMixture;

    int nPhases;
    bool STOICH;
    double EquiValRatio;
    string FuelName;
    string Oxidizer;

    double BurntGasTemp;
    size_t GasFuelIndex;

    Tensor<string, 1> GasSpeciesNames;
    Tensor<double, 1> GasMixture;
    Tensor<double, 1> AirMixture;
    Tensor<double, 1> BurntGasMixture;
    Tensor<double, 1> BackFlowMixture;
    Tensor<double, 1> MW;
    size_t GasIdx;
};
}
#endif

#endif
