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

#ifndef SPECIESTRANSPORT_H_INCLUDED
#define SPECIESTRANSPORT_H_INCLUDED

#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Velocities.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <array>
#include <random>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <tuple>

#include "Initializations.h"

#include "EnergyTransport.h"
#include "FlowMixture.h"
#include "SolidBody.h"

using namespace std;

namespace openphase { class FlowMixture;
class EnergyTransport;
class SolidBody; }
namespace openphase  {

class SpeciesTransport 
{
    public:
	SpeciesTransport(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)                                  ///< Constructor
    {
		Initialize(locSettings);
		ReadInput(InputFileName);
    }
	void ReadInput(string InputFile);
	void Initialize(Settings& locSettings);
	void CalculateSpeciesHeatCapacitiesAndEnthalpies(EnergyTransport& ET, FlowSolverLBM& FL);
	void UpdateFields(EnergyTransport& ET);
	void SetBoundaryConditions(EnergyTransport& ET, BoundaryConditions& BC);
	void SetFreeBCNX(PhaseField& Phase, EnergyTransport& ET);
	void CalculateAdvection(EnergyTransport& ET, FlowSolverLBM& FL, FlowMixture& FM, double dt);
	void CalculateDiffusion(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, double dt);
	void CalculateReaction(EnergyTransport& ET, FlowSolverLBM& FL, int nDim, double dt);

	void UpdateGhostPoints(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, SolidBody& SB);
	void CalculateGhostPoints(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, SolidBody& SB);
	void CalculateGhostPointsConjugateHeatTransfer(PhaseField& Phase, EnergyTransport& ET, FlowSolverLBM& FL, SolidBody& SB);

	void SetInitial(PhaseField& Phase, FlowSolverLBM& FL, EnergyTransport& ET,
                    const vector<double>& BM, const vector<double>& FM , const vector<double>& AM, 
                    double TB, double TF, double TA);

	void WriteVTKMassFractions(Settings& locSettings, const int tStep, const int precision = 16);
	void WriteVTKHHR(Settings& locSettings, const int tStep, const int precision = 16);
 
	double CalculateMeanMolarMass(const int i, const int j, const int k, string whichtime);
	double MoleFraction( const int i, const int j, const int k, double MMW, size_t kc);

    GridParameters Grid;                                                        ///< Simulation grid parameters

	bool IF_REACTION;
	size_t nSpecies;

	double X0BurntZone, XNBurntZone;
	double Y0BurntZone, YNBurntZone;
	double Z0BurntZone, ZNBurntZone;

	double X0FuelZone, XNFuelZone;
	double Y0FuelZone, YNFuelZone;
	double Z0FuelZone, ZNFuelZone;

	int BCOrder;
	double TempBurntGas;

	bool AdvectionUpwind; 
	bool AdvectionCentral;
	bool AdvectionVanLeer;
	double CalculatePhiVanLeer(double r);

	bool Species;
    Storage3D<double,0> HRR;                                                   ///< Heat release rate of reaction 
    Storage3D<double,1> Cp_Species;    
    Storage3D<double,1> h_Species;    
	Storage3D<double,1> MassDiff_Species;      
	Storage3D<double,1> W_Species;  
	Storage3D<double,1> MassFractions;  
	Storage3D<double,1> MassFractionsOld;  

	size_t nCoeffs=15;

    Tensor<double, 1> Y0X;                             
    Tensor<double, 1> MolecularWeight;     
	Tensor<double, 2> PolyNomCoeffs;

	double FuelConsumptionRate=0.0;
	size_t FuelIndex;
};
}
#endif