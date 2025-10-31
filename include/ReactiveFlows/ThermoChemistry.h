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

 *   File created :   2025
 *   Main contributors :   Oleg Shchyglo; Reza Namdar
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

#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include "cantera/base/Array.h"
#include "cantera/base/Solution.h"

#include "SpeciesTransport.h"
#include "EnergyTransport.h"

using namespace std;

namespace openphase {class SpeciesTransport;
class EnergyTransport;}

namespace openphase {
class ThermoChemistry 
{
    public:
    ThermoChemistry(const std::string InputFileName)                                   ///< Constructor
    {
        ReadInput(InputFileName);
    }
    void ReadInput(string InputFile);
    void WriteSpeciesNamesInSettingsInputFile(string InputFile);
    
    void Initialize(shared_ptr<Cantera::ThermoPhase> & gas, string SettingsInputFile);
    void UpdateProperties(EnergyTransport& ET, SpeciesTransport& EST, FlowSolverLBM& FL,
				shared_ptr<Cantera::ThermoPhase> &gas, shared_ptr<Cantera::Transport> &transp, shared_ptr<Cantera::Kinetics> &kin) ;
    void UpdateAllProperties(EnergyTransport& ET, SpeciesTransport& EST, FlowSolverLBM& FL,
				shared_ptr<Cantera::ThermoPhase> &gas, shared_ptr<Cantera::Transport> &transp, shared_ptr<Cantera::Kinetics> &kin) ;
    vector<double> GettingFuelMixture(std::shared_ptr<Cantera::ThermoPhase> & gas, string FuelName, string Oxidizer, double ER);
    vector<double> GettingAirMixture();
    vector<double> GettingBackFlowMixture();
    vector<double> GettingBurntMixture(std::shared_ptr<Cantera::ThermoPhase> & gas, double Temp, double Pressure);
    vector<double> GettingMolecularweight(std::shared_ptr<Cantera::ThermoPhase> & gas);
    vector<vector<double>> GettingSpeciesPolyConstants(std::shared_ptr<Cantera::ThermoPhase> & gas);
    vector<string> GettingSpeciesNames(std::shared_ptr<Cantera::ThermoPhase> & gas);

    size_t GettingFuelIndex(std::shared_ptr<Cantera::ThermoPhase> & gas, string FuelName);
    void ExportingData(SpeciesTransport& ST);

	vector<vector<double>>  Coeffs;
    size_t nCoeffs;

    string ReactionMechanism;
    string PhaseName;
    string TransportName;
    string FuelName;
    string Oxidizer;

    vector<string> SpeciesNames;

    vector<double> FuelMixture;
    vector<double> AirMixture;
    vector<double> BurntMixture;
    vector<double> BackFlowMixture;
    vector<double> MW;

    double ER;
    double BurntTemp;

    int nPhases;

    size_t FuelIndex;
};
}
#endif

#endif
