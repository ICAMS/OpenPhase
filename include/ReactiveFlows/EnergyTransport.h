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

#ifndef ENERGYTRANSPORT_H_INCLUDED
#define ENERGYTRANSPORT_H_INCLUDED

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

#include "FlowMixture.h"
#include "SolidBody.h"

using namespace std;

namespace openphase { class FlowMixture;
class SolidBody; }
namespace openphase  {

class EnergyTransport 
{
    public:
	EnergyTransport(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)                                  ///< Constructor
    {
		Initialize(locSettings);
		ReadInput(InputFileName);
    }
	void ReadInput(string InputFile);
	void Initialize(Settings& locSettings);

	void UpdateFields();
	void SetBoundaryConditions(BoundaryConditions& BC);
	void SetFreeBCNX(PhaseField& Phase);
	void SetSolidPhaseTemp(PhaseField& Phase, bool DI);
	void CalculateAdvection(FlowSolverLBM& FL, FlowMixture& FM, double dt);
	void CalculateDiffusion(PhaseField& Phase, FlowSolverLBM& FL, double dt);
	void CalculatePhaseFieldCoupling(PhaseField& Phase, FlowSolverLBM& FL, double hT, double dt);
	void CalculateThermodynamicPressure(FlowSolverLBM& FL, double P0, double T0);
	void UpdateGhostPoints(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB);
	void CalculateGhostPoints(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB);
	void CalculateGhostPointsConjugateHeatTransfer(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB);
	void SetInitial(PhaseField& Phase, FlowSolverLBM& FL);
	void WriteVTKTemperature(Settings& locSettings, FlowSolverLBM& FL, const int tStep, const int precision = 16);
	void IntegrateTemperature(FlowSolverLBM& FL);
	double Lagrange_polynomial4 (double T1, double T2, double T3,  double T4, double T5, double f1, double f2, double f3,  double f4);
	void CalculateMixtureSpecificHeatCapacityAndThermalConductivity(FlowSolverLBM& FL);
	void CalculateSolidDiffusion(PhaseField& Phase, FlowSolverLBM& FL, double dt);

    GridParameters Grid;                                                        ///< Simulation grid parameters
	double T0;
    Storage3D<double,0> Tx;                                                 
    Storage3D<double,0> TxOld;                                        

	vector<double> SurfaceTemp;
	vector<double> SurfaceFlux;
	vector<bool>   IF_ConstTemp;
	vector<bool>   IF_ConstFlux;

	bool IF_ENERGY;

	int BCOrder;

	double TempBurntGas;
	double TempSolid;

	double TempSolidCold;
	double HeatFlux;
	double TempBackFlow;
	double ATstar;
	double Pr;
	double Cp;
	double Mw;

	bool AdvectionUpwind; 
	bool AdvectionCentral;
	bool AdvectionVanLeer;
	double CalculatePhiVanLeer(double r);

    Storage3D<double,0> Cp_Mixture;                                              
    Storage3D<double,0> K_Mixture;    
	double K0; ///thermal conductivity at film temperature                                                 

	bool Conjugate;
	
    Storage3D<double,0> Ts;                                                   ///< Temperature store for conjugated heat tranfer
    Storage3D<double,0> TsOld;                                                ///< Temperature store for conjugated heat tranfer
	Storage3D<double,0> Cp_Solid;                                              
    Storage3D<double,0> K_Solid;   
    Storage3D<double,0> Density_Solid;  
	double Rhos;
	double Ks;
	double Cps;
	double Ts0;

	//Sutherland's constants
	double TMu0;
	double SMu0;
	double Mu0;

};
}
#endif