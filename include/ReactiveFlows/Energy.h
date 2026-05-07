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

#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED

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
#include <cmath>
#include <fstream>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <tuple>
#include "Initializations.h"

#include "MixtureFlow.h"
#include "SolidBody.h"

using namespace std;

namespace openphase { 
class MixtureFlow;
class SolidBody; 
 }
namespace openphase  {

struct ThermalCondition
{
    enum class Type
    {
        ConstantTemp,				///< Constant temperature applied on solid surface		
        ConstantFlux,				///< Constant flux applied on solid surface
        Conjugate    				///< Conjugate heat transfer
    } type;

    double value = 0.0;     // temp or flux
    size_t solidIdx = 0;    // used only for Conjugate
    double visualTemp = 0.0;
};

class OP_EXPORTS Energy 
{
    public:
	Energy(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)                                  ///< Constructor
    {
		ReadInput(InputFileName);
		Initialize(locSettings);
    }
	void Initialize(Settings& locSettings);
	void ReadInput(string InputFile);
	void UpdateFields();
	void SetBoundaryConditions(BoundaryConditions& BC);
	void SetOpenBoundary(const dVector3& A, const dVector3& B, int n);

	void SetSolidPhaseTemp(PhaseField& Phase, bool DI);
	void SetInitialNeighboringFluidPoints(PhaseField& Phase, SolidBody &SB, bool DI);
	void CalculateThermodynamicPressure(FlowSolverLBM& FL, double P0, double T0);

	void SetInitial(PhaseField& Phase, FlowSolverLBM& FL, MixtureFlow &MF, SolidBody &SB);
	double CalculateSutherlandViscosity(double Temp);
	double CalculateIdealGasDensity(double p, double Rm, double Temp);
	void CalculateProperties(PhaseField &Phase, FlowSolverLBM& FL, SolidBody &SB, MixtureFlow &MF, BoundaryConditions &BC);
	void WriteVTKTemperature(Settings& locSettings, const int tStep, const int precision = 16);
	void WriteVTKSolidTemperature(Settings& locSettings, const int tStep, const int precision = 16);

	double Lagrange_polynomial4 (double T1, double T2, double T3,  double T4, double T5, double f1, double f2, double f3,  double f4);

	void WriteVTKHHR(Settings& locSettings, const int tStep, const int precision = 16);

	double SurfaceTemperature(PhaseField& Phase, FlowSolverLBM &FL, SolidBody& SB, 
							  const dVector3& Xb, const dVector3& n, const string time);
	double InterfaceTemperature(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB, const int i, const int j, const int k);
	double CalculatePhiVanLeer(double r);
    AdvectionSchemes AdvScheme;
    GridParameters Grid;                                                        ///< Simulation grid parameters
	double T0;
    Storage3D<double,0> T;                                                 
    Storage3D<double,0> TOld;                                        

	vector<ThermalCondition> ThermalSurfaceCondition;

	bool ENERGY;
	bool REACTION;
	bool CONJUGATE;

	int BCOrder;

	double TempBurntGas;
	double TempSolid;

	double HeatFlux=0.0;

	double TempSolidCold;
	double TempBackFlow;
	double ATstar;

	double Pr;
	double Cp;
	double Mw;
    Storage3D<double,0> CpMix;                                              
    Storage3D<double,0> KMix;    
	double K0; ///thermal conductivity at film temperature                                                 


    Storage3D<double,0> Ts;                                                   ///< Temperature store for conjugated heat tranfer
    Storage3D<double,0> TsOld;                                                ///< Temperature store for conjugated heat tranfer
	Storage3D<double,0> CpSol;                                              
    Storage3D<double,0> KSol;   
    Storage3D<double,0> DensitySol;  

	double Rhos;
	double Ks;
	double Cps;
	double Ts0;

	//Sutherland's constants
	double TMu0;
	double SMu0;
	double Mu0;
	
	double TrackingValue=0.0;

    Storage3D<double,0> HRR;                                                   ///< Heat release rate of reaction 

};
}
#endif
