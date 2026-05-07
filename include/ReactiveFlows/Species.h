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

#ifndef SPECIES_H_INCLUDED
#define SPECIES_H_INCLUDED

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

#include "Energy.h"
#include "MixtureFlow.h"
#include "SolidBody.h"

using namespace std;

namespace openphase { 
class MixtureFlow;
class Energy;
class SolidBody;
 }
namespace openphase  {

class OP_EXPORTS Species
{
    public:
	Species(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)                                  ///< Constructor
    {
		ReadInput(InputFileName);
		Initialize(locSettings);
    }
	void ReadInput(string InputFile);
	void Initialize(Settings& locSettings);
	void CalculateCpSpAndhSp(Energy& EN, FlowSolverLBM& FL);
	void UpdateFields();
	void SetBoundaryConditions(BoundaryConditions& BC);
	void SetOpenBoundary(const dVector3& A, const dVector3& B, int n);
	double Lagrange_polynomial4(double T1, double T2, double T3, double T4, double T5,
                                          double f1, double f2, double f3, double f4);
	Tensor<double, 1> SurfaceSpecies(FlowSolverLBM& FL, SolidBody& SB,
                                     const dVector3& Xb, const dVector3& n, const string time); 
    Tensor<double, 1> InterfaceSpecies(FlowSolverLBM& FL, SolidBody& SB, const int i, const int j, const int k);
	void SetInitial(PhaseField& Phase, FlowSolverLBM& FL, Energy& EN,
                    const Tensor<double, 1>& BM, const Tensor<double, 1>& FM , const Tensor<double, 1>& AM, 
                    double TB, double TF, double TA);

	void WriteVTK(Settings& locSettings, const int tStep, const int precision = 16);
 
	double CalculateMeanMolarMass(const int i, const int j, const int k, string whichtime);
	double MoleFraction( const int i, const int j, const int k, double MMW, size_t kc);

    GridParameters Grid;                                                        ///< Simulation grid parameters

	bool REACTION;
	size_t nSpecies;

	double X0BurntZone, XNBurntZone;
	double Y0BurntZone, YNBurntZone;
	double Z0BurntZone, ZNBurntZone;

	double X0FuelZone, XNFuelZone;
	double Y0FuelZone, YNFuelZone;
	double Z0FuelZone, ZNFuelZone;

	int BCOrder;
	double TempBurntGas;

	bool SPECIES;
    Storage3D<double,1> CpSp;    
    Storage3D<double,1> hSp;    
	Storage3D<double,1> DSp;      
	Storage3D<double,1> WSp;  
	Storage3D<double,1> MassFrac;  
	Storage3D<double,1> MassFracOld;  

	size_t nCoeffs;

    Tensor<double, 1> MwSp;     
	Tensor<double, 2> PolyNomCoeffs;
};
}
#endif
